!====================================================================== 
!    Routines for simple mesh deformations
!    Author:   Prabal Negi
!    Last Modified: 31-Jan-2016
!======================================================================       
      subroutine mvmsh_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'WING_MVMSH'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /WING_MVMSH/ axis_x0,axis_y0,ptch_kred,ptch_amp,
     $         def_st,def_end,ini_aoa,mw_def,ptch_start,msh_rescale


!     default values
      axis_x0         = 0.35
      axis_y0         = 0.034
      ptch_kred       = 0.5
      ptch_amp        = 2.0
      def_st          = 0.1
      def_end         = 1.8 
      ini_aoa         = 0.0 
      mw_def          = 'mv '
      ptch_start      = 0.0  
      msh_rescale     = 0.

!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=WING_MVMSH,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading WING_MVMSH parameters.$')

!     broadcast data
      call bcast(axis_x0,       WDSIZE)
      call bcast(axis_y0,       WDSIZE)
      call bcast(ptch_kred,     WDSIZE)
      call bcast(ptch_amp,      WDSIZE)
      call bcast(def_st,        WDSIZE)
      call bcast(def_end,       WDSIZE)
      call bcast(ini_aoa,       WDSIZE)
      call bcast(mw_def,      3*CSIZE)   
      call bcast(ptch_start,    WDSIZE)
      call bcast(msh_rescale,   WDSIZE)

      return
      end
!-----------------------------------------------------------------------
!     write parameters relaxation term filtering 
      subroutine mvmsh_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'WING_MVMSH'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /WING_MVMSH/ axis_x0,axis_y0,ptch_kred,ptch_amp,
     $         def_st,def_end,ini_aoa,mw_def,ptch_start,msh_rescale

!     read the file
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=WING_MVMSH,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing WING_MVMSH parameters.$')

      return
      end
!----------------------------------------------------------------------

      subroutine my_mvmsh_init

!     Initialize blending function for mesh motion.
!     And set initial angle of attack            
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'               ! param(59)
      include 'WING_MVMSH'

      integer lt
      parameter (lt=lx1*ly1*lz1*lelv)
      
      common /scrns/  h1(lt),h2(lt),rhs(lt),msk(lt),tmp(lt)
      real h1,h2,rhs,msk,tmp
      logical ifaoa           ! if setting initial angle of attack

!     calculate blending function 'basev'
      if (ini_aoa.ne.0) then
        ifaoa=.true.  
!        call my_basev_hlm(basev,h1,h2,rhs,msk,tmp,def_st,mw_def,ifaoa)
        call my_basev_erf(basev,h1,h2,rhs,msk,tmp,
     $       def_st,def_end,mw_def)
        if (abs(msh_rescale).gt.1.0e-12) then
          call gll_rescale(msh_rescale)
          def_st = def_st*(1.+msh_rescale)
          def_end = def_end*(1.+msh_rescale)
          call my_basev_erf(basev,h1,h2,rhs,msk,tmp,
     $       def_st,def_end,mw_def)
        endif        
        call alpha_init2
        ifaoa=.false.  
!        call my_basev_hlm(basev,h1,h2,rhs,msk,tmp,def_st,mw_def,ifaoa)
        call my_basev_erf(basev,h1,h2,rhs,msk,tmp,
     $       def_st,def_end,mw_def)

      else
        call my_basev_erf(basev,h1,h2,rhs,msk,tmp,
     $       def_st,def_end,mw_def)
        if (abs(msh_rescale).gt.1.0e-12) then
          call gll_rescale(msh_rescale)
          def_st = def_st*(1.+msh_rescale)
          def_end = def_end*(1.+msh_rescale)
          call my_basev_erf(basev,h1,h2,rhs,msk,tmp,
     $       def_st,def_end,mw_def)
        endif        
!        call my_basev_hlm(basev,h1,h2,rhs,msk,tmp,def_st,mw_def)
      endif

      param(59) = 1           ! all elements deformed
      call geom_reset(1)

      return
      end subroutine my_mvmsh_init

!---------------------------------------------------------------------- 

      subroutine my_basev_erf(basev,h1,h2,rhs,msk,tmp,x1,x2,surface)

!     Set up base interpolation function.
!
!     This is the expensive part where we find a smooth
!     interpolant between the moving boundary and fixed boundaries.
!
!     Note that it is called only once (or twice).
!
!     Note also that you have other types of motions you might (or
!     might not) need an addition interpolating basis function.


      implicit none

      include 'SIZE_DEF'
      include 'SIZE'

      real basev(lx1,ly1,lz1,lelt),tmp(lx1,ly1,lz1,lelt)
      real h1(1),h2(1),rhs(1),msk(1)
      real h3(lx1*ly1*lz1*lelt)           ! diagnostics

      integer e,f,n,ifld,i

      real x2,x1,y2,y1,rr

      character(3) surface 

      n = nx1*ny1*nz1*nelv

      call rone (h1 ,n)

      call rzero(h3,n)   ! temporary diagnostics
      
c     Modify h1 to make blending function constant near the surface.
c     The idea here is to push the majority of the mesh deformation to the 
c     far field, away from where the fine, boundary-layer-resolving elements
c     are close to the surface.
c

      ifld = 1
      call cheap_dist(h1,ifld,surface)       ! calculate distance from defined "surface"

      y2 = -0.5 
      y1 = 2.0 

      do i=1,n
         
         h1(i) = max(h1(i),x1)
         h1(i) = min(h1(i),x2)
         rr = (y2-y1)/(x2-x1)*h1(i)+(y1*x2-y2*x1)/(x2-x1)
         h3(i) = rr

!     Perhaps more complicated than it needs to be.
!     Works for now.
         if (h1(i).gt.x2) then
               h1(i) = 0.
         else
            h1(i) = (erf(rr)-erf(y2))/(erf(y1)-erf(y2))
         endif
      enddo

      call copy(basev,h1,n)

!      call dsavg(basev)

c     BELOW JUST FOR DIAGNOSTICS
      call outpost(h3,h1,basev,basev,basev,'def')
c      call exitti('quit in my_base_meshv$',nelt)
            
      return
      end subroutine my_basev_erf

!---------------------------------------------------------------------- 

      subroutine my_basev_hlm(basev,h1,h2,rhs,msk,tmp,def_st,surface,
     $              ifaoa)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'

      real basev(lx1,ly1,lz1,lelt),tmp(lx1,ly1,lz1,lelt)
      real h1(1),h2(1),rhs(1),msk(1)
      real h3(lx1*ly1*lz1*lelt)           ! diagnostics
      real m1
      integer i,e,f,n,imsh,ifld,nface
      real tolold,tol

      real rr,arg,delta,z0
      real basev_ofst
      parameter (basev_ofst = 1.0)

      real def_st             ! start of deformation

      character(3) surface
      logical ifaoa 

      n = nx1*ny1*nz1*nelv

      call rone (h1 ,n)  ! H*u = -div (h1 grad u) + h2 u = 0
      call rzero(h2 ,n)  ! h2  = 0
      call rone (msk,n)  ! Mask, for Dirichlet boundary values
      call rzero(tmp,n)  ! RHS for Laplace solver

      call rzero(h3,n)   ! temporary diagnostics. Remove from
                         ! declaration as well
c
c     Modify h1 to make blending function constant near the cylinder.
c     The idea here is to push the majority of the mesh deformation to the 
c     far field, away from where the fine, boundary-layer-resolving elements
c     are close to the cylinder.
c
      ifld = 1
      call cheap_dist(h1,ifld,surface)       ! calculate distance from defined "surface"

!     Has to be manually set in a case dependent fashion.
!     if def_st is too large then the boundary layer formed
!     for mesh movement is in a coarse far field region.
!     Which has issues with convergence.
!     If its 0., then it has convergence issues in sharp corners.
      if (ifaoa) then
        def_st=0.15
      endif  

      delta = 0.5
      do i=1,n
         rr = h1(i)-def_st
         h3(i) = max(rr,0.)
         arg   = -h3(i)/(delta**2)
         h1(i) = 1. + 9.0*exp(arg)
      enddo

      m1 = -1.
      z0 =  0.

      nface = 2*ndim
      do e=1,nelv
        do f=1,nface
!         Set Dirichlet for mesh velocity on all non-interior boundaries
          if (cbc(f,e,1).ne.'E  '.and.(cbc(f,e,1).ne.'P  ')) 
     $      call facev(msk,e,f,z0,nx1,ny1,nz1)

!         Set inhomogeneous Dirichlet data on surface
          if (cbc(f,e,1).eq.surface) then
            call facev(tmp,e,f,m1,nx1,ny1,nz1)
          endif
        enddo
      enddo

      do i=1,n
        if (h3(i).le.0.) then
          msk(i)=0.
          tmp(i,1,1,1)=m1
        endif
      enddo

      call axhelm (rhs,tmp,h1,h2,1,1)

      tolold = param(22)
      tol    = 1.e-12
      param(22) = tol
      imsh   = 1
      ifld = ifield
      ifield = 1
      call hmholtz('mshv',basev,rhs,h1,h2,msk,vmult,imsh,tol,40000,1)

      call dsavg(basev)       ! make continuous at elements

      param(22) = tolold
      ifield = ifld

      call sub2(basev,tmp,n)  ! <--- THIS IS THE GREEN's FUNCTION

      call cmult(basev,1./basev_ofst,n)      
      do i=1,n
        basev(i,1,1,1)=min(basev(i,1,1,1),1.0)
      enddo

c     BELOW JUST FOR DIAGNOSTICS
      call outpost(h3,h1,basev,basev,basev,'def')
c      call exitti('quit in my_base_meshv$',nelt)
      
            
      return
      end subroutine my_basev_hlm

c-----------------------------------------------------------------------
      subroutine my_meshv

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'         ! istep, pi
      include 'INPUT_DEF'
      include 'INPUT'         ! red freq,amplitude of pitching
      include 'GEOM_DEF'
      include 'GEOM'          ! xm1,ym1,zm1
      include 'MVGEOM_DEF'
      include 'MVGEOM'        ! wx,wy,wz
      include 'WING_MVMSH'
     
      integer lt
      parameter (lt=lx1*ly1*lz1*lelv)

      real amp,omega,aa
      real x0,y0              ! pitch axis
      real x,y
      real ucx,ucy,ucz        ! solid body rotation velocities

      integer i,n
      real U0                 ! free stream U
      real diameter           ! cylinder diameter
      real phase_shift        ! if you don't start at mean alpha

!----------------------------------------  

      n = nx1*ny1*nz1*nelv

      pi = 4.0*atan(1.0)

      x0 = axis_x0 
      y0 = axis_y0

      U0 = 1.0
      diameter = 1.0
      if (time.ge.ptch_start) then
        amp = ptch_amp
      else
        amp = 0.0
      endif 

      omega         = ptch_kred*U0/diameter
      phase_shift   = 0.
      aa            = amp*omega*cos(omega*(time-ptch_start)+phase_shift)

      do i=1,n                          ! Translational velocity
        x = xm1(i,1,1,1) - x0
        y = ym1(i,1,1,1) - y0

        ucx =  aa*y
        ucy = -aa*x
        ucz =  0.0

        wx(i,1,1,1) = basev(i)*ucx     ! component.
        wy(i,1,1,1) = basev(i)*ucy
        wz(i,1,1,1) = basev(i)*ucz

        umeshx(i,1,1,1) = wx(i,1,1,1)
        umeshy(i,1,1,1) = wy(i,1,1,1)
        umeshz(i,1,1,1) = wz(i,1,1,1)
      enddo

c     BELOW JUST FOR DIAGNOSTICS
c      call outpost(wx,wy,wz,pr,basev,'   ')
c      if (mod(istep,100).eq.0) then
c      call outpost(wx,basev,wz,pr,wz,'   ')
c      call exitti('quit in my_meshv$',nelt)
c      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine my_meshv_eta

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'          ! xm1,ym1,zm1
      include 'MVGEOM_DEF'
      include 'MVGEOM'        ! wx,wy,wz
      include 'WING_MVMSH'
      include 'FSI'

      real ucx,ucy,ucz        ! mesh velocities
      real dx,dy

      integer i,n

      n = nx1*ny1*nz1*nelv

      do i=1,n                          ! Translational velocity
        if (.not.fsi_ifrot) then
!         Only consider vertical motion
          ucx =  0.
          ucy =  etav
          ucz =  0.0

          wx(i,1,1,1) = basev(i)*ucx
          wy(i,1,1,1) = basev(i)*ucy
          wz(i,1,1,1) = basev(i)*ucz

          umeshx(i,1,1,1) = wx(i,1,1,1)
          umeshy(i,1,1,1) = wy(i,1,1,1)
          umeshz(i,1,1,1) = wz(i,1,1,1)
        else
          dx = xm1(i,1,1,1) - fsi_x0
          dy = ym1(i,1,1,1) - fsi_y0

          ucx =  etav*dy
          ucy = -etav*dx
          ucz =  0.0

          wx(i,1,1,1) = basev(i)*ucx
          wy(i,1,1,1) = basev(i)*ucy
          wz(i,1,1,1) = basev(i)*ucz

          umeshx(i,1,1,1) = wx(i,1,1,1)
          umeshy(i,1,1,1) = wy(i,1,1,1)
          umeshz(i,1,1,1) = wz(i,1,1,1)
        endif
      enddo

      return
      end subroutine my_meshv_eta

c-----------------------------------------------------------------------

      subroutine alpha_init

!     Routine needs basev to have been calculated.
!     Must be called after mvmsh_init           
            
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'WING_MVMSH'
      include 'PARALLEL_DEF'
      include 'PARALLEL'

      integer i,n
      real x0,y0              ! pitch axis
      real x(2)               ! temporary x,y
      real xnew(2)            ! rotated x,y
      real theta              ! initial angle of attack
      real rot_mat(2,2)       ! Rotational matrix
      real pi

      pi    =     3.14159265359

      n = nx1*ny1*nz1*nelv

      x0 = axis_x0 
      y0 = axis_y0

      theta        = ini_aoa*pi/180.

!     Rotational matrix for a clockwise (+ve angle of attack)
!     rotation by an angle theta (in radians)
      rot_mat(1,1) = cos(theta)
      rot_mat(2,1) = -sin(theta)
      rot_mat(1,2) = sin(theta)
      rot_mat(2,2) = cos(theta)
!--------------------

      do i=1,n

         x(1) = xm1(i,1,1,1) - x0
         x(2) = ym1(i,1,1,1) - y0

         call mxm(rot_mat,2,x,2,xnew,1)

         xnew(1) = (xnew(1)-x(1))*basev(i) + x(1)
         xnew(2) = (xnew(2)-x(2))*basev(i) + x(2)

         xm1(i,1,1,1) = xnew(1) + x0
         ym1(i,1,1,1) = xnew(2) + y0

      enddo

      if (nio.eq.0) write(6,*) 'Surface rotated by ',
     $            ini_aoa, ' degrees' 

      return
      end subroutine alpha_init
!----------------------------------------------------------------------
      subroutine alpha_init2

!     Routine needs basev to have been calculated.
!     Must be called after mvmsh_init           
            
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'WING_MVMSH'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'

      integer i,n
      real x0,y0              ! pitch axis
      real x(2)               ! temporary x,y
      real xnew(2)            ! rotated x,y
      real theta              ! initial angle of attack
      real rot_mat(2,2)       ! Rotational matrix
      real pi

      real deltax(lx1,ly1,lz1,lelt),deltay(lx1,ly1,lz1,lelt)
      real dxmask(lx1,ly1,lz1,lelt)
      common /psn_scrtch/ deltax,deltay,dxmask

      integer iface,ifld
      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz

      CHARACTER CB*3
      integer Iel,NFACES,npts
     
      pi    =     3.14159265359

      n = nx1*ny1*nz1*nelv

      x0 = axis_x0 
      y0 = axis_y0

      theta        = ini_aoa*pi/180.

!     Rotational matrix for a clockwise (+ve angle of attack)
!     rotation by an angle theta (in radians)
      rot_mat(1,1) = cos(theta)
      rot_mat(2,1) = -sin(theta)
      rot_mat(1,2) = sin(theta)
      rot_mat(2,2) = cos(theta)
!--------------------

      do i=1,n

         x(1) = xm1(i,1,1,1) - x0
         x(2) = ym1(i,1,1,1) - y0

         call mxm(rot_mat,2,x,2,xnew,1)

         xnew(1) = (xnew(1)-x(1))*basev(i) + x(1)
         xnew(2) = (xnew(2)-x(2))*basev(i) + x(2)

         deltax(i,1,1,1) = xnew(1) + x0 - xm1(i,1,1,1)
         deltay(i,1,1,1) = xnew(2) + y0 - ym1(i,1,1,1)

      enddo

!     I'm trying to use the modified fix_geom to do the rotation.
!     Allows me to rotate the airfoil/object while maintaining inner
!     mapping of GLL points.
!     Hopefully that works.      
!     Creating a mask array
      NFACES=2*NDIM

      n = nx1*ny1*nz1*nelv
      call rzero(dxmask,n)

      ifld = 1
      do  Iel=1,NELV            !do ieg=1,nelgt
        do IFACE = 1,NFACES
          CB = CBC(IFACE,Iel,ifld)
          if (CB.EQ.MW_DEF) then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
!           On the surface boundary I move all points.            
            do iz=KZ1,KZ2
              do iy=KY1,KY2
                do ix=KX1,KX2
                  dxmask(ix,iy,iz,Iel)=1.0
                enddo
              enddo
            enddo
          endif  
          if (CB.EQ.'E  '.or.CB.eq.'P  ') then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
!           On the internal/Parallel surfaces            
!           Only take the deformations on the vertices. 
            dxmask(KX1,KY1,KZ1,Iel)=1.0
            dxmask(KX1,KY1,KZ2,Iel)=1.0
            dxmask(KX1,KY2,KZ1,Iel)=1.0
            dxmask(KX1,KY2,KZ2,Iel)=1.0
            dxmask(KX2,KY1,KZ1,Iel)=1.0
            dxmask(KX2,KY1,KZ2,Iel)=1.0
            dxmask(KX2,KY2,KZ1,Iel)=1.0
            dxmask(KX2,KY2,KZ2,Iel)=1.0
          endif

          if (CB.EQ.'v  '.or.CB.eq.'V  '.or.CB.eq.'o  '
     $       .or.CB.eq.'O  ') then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
!           On the external surfaces ensure that the deformation is zero     
            dxmask(KX1,KY1,KZ1,Iel)=0.0
            dxmask(KX1,KY1,KZ2,Iel)=0.0
            dxmask(KX1,KY2,KZ1,Iel)=0.0
            dxmask(KX1,KY2,KZ2,Iel)=0.0
            dxmask(KX2,KY1,KZ1,Iel)=0.0
            dxmask(KX2,KY1,KZ2,Iel)=0.0
            dxmask(KX2,KY2,KZ1,Iel)=0.0
            dxmask(KX2,KY2,KZ2,Iel)=0.0
          endif
        enddo     ! iface
      enddo       ! nelv

      call dsavg(dxmask)     
      call col2(deltax,dxmask,n)
      call col2(deltay,dxmask,n)
      call rzero(dxmask,n)          ! now using as deltaz array.
      
      call mv_internal_pts(deltax,deltay,dxmask)


      if (nio.eq.0) write(6,*) 'Surface rotated by ',
     $            ini_aoa, ' degrees' 

      return
      end subroutine alpha_init2
!----------------------------------------------------------------------

      subroutine gll_rescale(fact)

!     Routine needs basev to have been calculated.
!     Must be called after mvmsh_init           
            
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'WING_MVMSH'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'

      integer i,n
      real x0,y0              ! pitch axis
      real x(2)               ! temporary x,y
      real xnew(2)            ! rotated x,y
      real theta              ! initial angle of attack
      real rot_mat(2,2)       ! Rotational matrix
      real pi

      real deltax(lx1,ly1,lz1,lelt),deltay(lx1,ly1,lz1,lelt)
      real dxmask(lx1,ly1,lz1,lelt)
      common /psn_scrtch/ deltax,deltay,dxmask

      integer iface,ifld
      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz

      CHARACTER CB*3
      integer Iel,NFACES,npts

      real fact
     
      pi    =     3.14159265359

      n = nx1*ny1*nz1*nelv

      x0 = axis_x0 
      y0 = axis_y0

      theta        = ini_aoa*pi/180.

!     Rotational matrix for a clockwise (+ve angle of attack)
!     rotation by an angle theta (in radians)
      rot_mat(1,1) = cos(theta)
      rot_mat(2,1) = -sin(theta)
      rot_mat(1,2) = sin(theta)
      rot_mat(2,2) = cos(theta)
!--------------------

      do i=1,n

         xnew(1) = xm1(i,1,1,1)*(1. + 
     $         fact*(1.-basev(i)))
         xnew(2) = ym1(i,1,1,1)*
     $          (1+fact*(1.-basev(i)))

         deltax(i,1,1,1) = xnew(1) - xm1(i,1,1,1)
         deltay(i,1,1,1) = xnew(2) - ym1(i,1,1,1)

      enddo

!     I'm trying to use the modified fix_geom to rescale the mesh
!     Hopefully that works.      
!     Creating a mask array
      NFACES=2*NDIM

      n = nx1*ny1*nz1*nelv
      call rzero(dxmask,n)

      ifld = 1
      do  Iel=1,NELV            !do ieg=1,nelgt
        do IFACE = 1,NFACES
          CB = CBC(IFACE,Iel,ifld)
          if (CB.EQ.MW_DEF) then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
!           On the surface boundary I move all points.            
            do iz=KZ1,KZ2
              do iy=KY1,KY2
                do ix=KX1,KX2
                  dxmask(ix,iy,iz,Iel)=1.0
                enddo
              enddo
            enddo
          endif  
          if (CB.EQ.'E  '.or.CB.eq.'P  ') then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
!           On the internal/Parallel surfaces            
!           Only take the deformations on the vertices. 
            dxmask(KX1,KY1,KZ1,Iel)=1.0
            dxmask(KX1,KY1,KZ2,Iel)=1.0
            dxmask(KX1,KY2,KZ1,Iel)=1.0
            dxmask(KX1,KY2,KZ2,Iel)=1.0
            dxmask(KX2,KY1,KZ1,Iel)=1.0
            dxmask(KX2,KY1,KZ2,Iel)=1.0
            dxmask(KX2,KY2,KZ1,Iel)=1.0
            dxmask(KX2,KY2,KZ2,Iel)=1.0
          endif

          if (CB.EQ.'v  '.or.CB.eq.'V  '.or.CB.eq.'o  '
     $       .or.CB.eq.'O  ') then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
!           We need to move external boundaries as well 
            dxmask(KX1,KY1,KZ1,Iel)=1.0
            dxmask(KX1,KY1,KZ2,Iel)=1.0
            dxmask(KX1,KY2,KZ1,Iel)=1.0
            dxmask(KX1,KY2,KZ2,Iel)=1.0
            dxmask(KX2,KY1,KZ1,Iel)=1.0
            dxmask(KX2,KY1,KZ2,Iel)=1.0
            dxmask(KX2,KY2,KZ1,Iel)=1.0
            dxmask(KX2,KY2,KZ2,Iel)=1.0
          endif
        enddo     ! iface
      enddo       ! nelv

      call dsavg(dxmask)     
      call col2(deltax,dxmask,n)
      call col2(deltay,dxmask,n)
      call rzero(dxmask,n)          ! now using as deltaz array.
      
      call mv_internal_pts(deltax,deltay,dxmask)

      if (nio.eq.0) write(6,*) 'Surface rotated by ',
     $            ini_aoa, ' degrees' 

      return
      end subroutine gll_rescale
!----------------------------------------------------------------------

