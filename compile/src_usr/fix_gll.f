!!======================================================================
!!
!!   Author: Prabal Negi
!!   Description: Poisson solver for fixing the mesh
!!   Last Modified: 10-03-2018
!!
!!======================================================================
!---------------------------------------------------------------------- 
      subroutine fix_mygll()

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'

      real dum(3)

      CHARACTER CB*3
      integer Iel,NFACES,npts
      parameter(npts=1000001)
      real xnew,ynew
      real xexact(npts),yexact(npts)

      real deltax(lx1,ly1,lz1,lelt),deltay(lx1,ly1,lz1,lelt)
      real deltaz(lx1,ly1,lz1,lelt)
      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      integer iface,ifld
      integer i,n
      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz

      real psx(lt),psy(lt)

      n=nx1*ny1*nz1*nelt

      call rzero(deltax,lt)
      call rzero(deltay,lt)

      NFACES=2*NDIM

      ifxyo = .true.
      call outpost (vx,vy,vz,pr,t,'gri')
      n = nx1*ny1*nz1*nelv

!     Read GLL points (2D) from file
      open(unit=19,file='naca4412.dat')
      do i=1,npts
        read(19,*) xexact(i),yexact(i)
      enddo
      close (19)

      ifld = 1
      do  Iel=1,NELV            !do ieg=1,nelgt
        do IFACE = 1,NFACES
          CB = CBC(IFACE,Iel,ifld)
          if (CB.EQ.'W  ') then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
            do iz=KZ1,KZ2
              do iy=KY1,KY2
                do  ix=KX1,KX2
                  if (xm1(ix,iy,iz,iel).lt.0.99) then
                    call fix_naca_pts(xnew,ynew,xm1(ix,iy,iz,Iel),
     $              ym1(ix,iy,iz,Iel),xexact,yexact)
                    deltax(ix,iy,iz,Iel)=xnew-xm1(ix,iy,iz,iel)
                    deltay(ix,iy,iz,Iel)=ynew-ym1(ix,iy,iz,iel)
                  else  
                    deltax(ix,iy,iz,Iel)=0.
                    deltay(ix,iy,iz,Iel)=0.
                  endif
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      call my_fixmesh(psx,psy,deltax,deltay)
      call add2(xm1,psx,lt)
      call add2(ym1,psy,lt)
      call fix_geom

      call copy(vx,psx,lt)
      call copy(vy,psy,lt)

      call outpost (vx,vy,vz,pr,t,'gri')

      return
      end

!---------------------------------------------------------------------- 

      subroutine fix_naca_pts(xnew,ynew,xfoil,yfoil,xexact,yexact)

      real*8 xfoil,yfoil,xnew,ynew,dist_min
      integer npts,counter,i
      parameter(npts=1000001)

      real*8 dist

      real*8 xexact(1),yexact(1)

      dist_min = 1.
c     Computed the distance for each point

      do i=1,npts
        dist = sqrt(abs(xexact(i)-xfoil)**2+abs(yexact(i)-yfoil)**2)
        if (dist.lt.dist_min) then
          counter = i
          dist_min = dist
        endif
      enddo

      xnew = xexact(counter)
      ynew = yexact(counter)

      return
      end

c-----------------------------------------------------------------------
      subroutine my_fixmesh(psx,psy,usrfldx,usrfldy)

!     Initialize blending function for mesh motion. 
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'

      integer lt
      parameter (lt=lx1*ly1*lz1*lelv)
      
      common /scrns/  h1(lt),h2(lt),rhs(lt),msk(lt),tmp(lt)
      real h1,h2,rhs,msk,tmp

      real psx(lt),psy(lt)
      real usrfldx(lt),usrfldy(lt)

      call my_poss_soln(psx,psy,h1,h2,rhs,msk,tmp,'W  ',usrfldx,usrfldy)

      return
      end subroutine my_fixmesh

!---------------------------------------------------------------------- 

      subroutine my_poss_soln(psx,psy,h1,h2,rhs,msk,tmp,surface,
     $      usrfldx,usrfldy)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'SOLN_DEF'
      include 'SOLN'          ! vmult

      real tmp(lx1,ly1,lz1,lelt)
      real h1(1),h2(1),rhs(1),msk(1)
      real h3(lx1*ly1*lz1*lelt)           ! diagnostics

      real m1
      integer i,e,f,n,imsh,ifield,ifld

      character(3) surface
      real rr,arg,delta,z0,tol,xavg,tolold

      integer nface

      real psx(1),psy(1)                ! output
      real usrfldx(1),usrfldy(1)        ! BCs

      n = nx1*ny1*nz1*nelv

      call rone (h1 ,n)  ! H*u = -div (h1 grad u) + h2 u = 0
      call rzero(h2 ,n)  ! h2  = 0
      call rone (msk,n)  ! Mask, for Dirichlet boundary values
      call rzero(tmp,n)  ! RHS for Laplace solver

      call rzero(h3,n)   ! temporary diagnostics. Remove from
                         ! declaration as well
c
c     Modify h1 to make blending function constant near the surface.
c     The idea here is to push the majority of the mesh deformation to the 
c     far field, away from where the fine, boundary-layer-resolving elements
c     are close to the cylinder.

      ifld = 1
      call cheap_dist(h1,ifld,surface)       ! calculate distance from defined "surface"
      delta = 0.5
      do i=1,n
        rr = h1(i)
        h3(i) = rr
        arg   = -rr/(delta**2)
        h1(i) = 1. + 9.0*exp(arg)
      enddo

      z0 =  0.

      nface = 2*ndim
      do e=1,nelv
        do f=1,nface
c         Set Dirichlet for mesh velocity on all non-interior boundaries
          if (cbc(f,e,1).ne.'E  '.and.cbc(f,e,1).ne.'P  ') 
     $         call facev(msk,e,f,z0,nx1,ny1,nz1)
        enddo
      enddo

      tol    = 1.e-16
      imsh   = 1
      ifield = 1
      
      tolold = param(22)
      param(22) = tol

!     deltax      
      call copy(tmp,usrfldx,n)
      call chsign(tmp,n)
      call axhelm (rhs,tmp,h1,h2,1,1)
      call hmholtz('mshv',psx,rhs,h1,h2,msk,vmult,imsh,tol,200000,1)
      call sub2(psx,tmp,n)

      call dsavg(psx)       ! This ensures periodic points have exactly the same deltax.
                            ! Should also remove tears from meshes  

!     deltay
      call copy(tmp,usrfldy,n)
      call chsign(tmp,n)
      call axhelm (rhs,tmp,h1,h2,1,1)
      call hmholtz('mshv',psy,rhs,h1,h2,msk,vmult,imsh,tol,200000,1)
      call sub2(psy,tmp,n)

      call dsavg(psy)       ! This ensures periodic points have exactly the same deltay.
                            ! Should also remove tears from meshes  

      param(22) = tolold
           
      return
      end subroutine my_poss_soln

!-----------------------------------------------------------------------
      subroutine fix_naca_gll(series)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'

      real dum(3)

      CHARACTER CB*3
      integer Iel,NFACES,npts
      real xnew,ynew

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      real deltax(lx1,ly1,lz1,lelt),deltay(lx1,ly1,lz1,lelt)
      real deltaz(lx1,ly1,lz1,lelt)
      common /psn_scrtch/ deltax,deltay,deltaz

      integer iface,ifld
      integer i,n
      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz
      real chord,rescale
      real foilx_max,foilx_min
      real glmax,glmin

      integer series

!      series=0012

      n=nx1*ny1*nz1*nelt

      call rzero(deltax,lt)
      call rzero(deltay,lt)
      call rzero(deltaz,lt)

      NFACES=2*NDIM

      ifxyo = .true.
      call outpost (vx,vy,vz,pr,t,'gri')
      n = nx1*ny1*nz1*nelv

      ifld = 1
      do  Iel=1,NELV            !do ieg=1,nelgt
        do IFACE = 1,NFACES
          CB = CBC(IFACE,Iel,ifld)
          if (CB.EQ.'mv ') then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
            do iz=KZ1,KZ2
              do iy=KY1,KY2
                do  ix=KX1,KX2
                  deltax(ix,iy,iz,Iel)=xm1(ix,iy,iz,iel)
                  deltay(ix,iy,iz,Iel)=ym1(ix,iy,iz,iel)
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      foilx_max = glmax(deltax,lt)
      foilx_min = glmin(deltax,lt)
      call cadd(xm1,-foilx_min,lt)         ! leading edge is at x=0 

      chord = foilx_max-foilx_min
      rescale = 1./chord
      call cmult(xm1,rescale,lt)           ! chord length is exactly 1.
      call cmult(ym1,rescale,lt)

      do  Iel=1,NELV            !do ieg=1,nelgt
        do IFACE = 1,NFACES
          CB = CBC(IFACE,Iel,ifld)
          if (CB.EQ.'mv ') then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
            do iz=KZ1,KZ2
              do iy=KY1,KY2
                do  ix=KX1,KX2
                  call fix_naca_bdry(series,xm1(ix,iy,iz,iel),
     $                   ym1(ix,iy,iz,iel),xnew,ynew)
                  deltax(ix,iy,iz,Iel)=xnew-xm1(ix,iy,iz,iel)
                  deltay(ix,iy,iz,Iel)=ynew-ym1(ix,iy,iz,iel)
!                    write(6,'(A3,1x,4(E15.7E2,1x))') 'dx:', 
!     $                  deltax(ix,iy,iz,iel),
!     $                  deltay(ix,iy,iz,iel), 
!     $                  xm1(ix,iy,iz,iel),
!     $                  ym1(ix,iy,iz,iel)                     
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      call outpost (deltax,deltay,deltaz,pr,t,'gri')
      call mv_internal_pts(deltax,deltay,deltaz)
      call outpost (deltax,deltay,deltaz,pr,t,'gri')

      return
      end

!---------------------------------------------------------------------- 
      subroutine fix_naca_bdry(series,xin,yin,xout,yout)

      implicit none

      integer series
      real t,m,p            ! airfoil parameters.
      save t,m,p
      
      logical iffinite_te
      parameter (iffinite_te=.false.)
      
      real a0,a1,a2,a3
      parameter (a0 =  0.2969)
      parameter (a1 = -0.1260)
      parameter (a2 = -0.3516)
      parameter (a3 =  0.2843)

      real a4                 ! depends on iffinite_te
      save a4

      real yt,yc,dycdx
      real xu,yu,xl,yl
      real theta

      real tmp,tmp2
      real xin,yin,xout,yout

      integer icalld
      save icalld
      data icalld /0/
      
      
      if (icalld.eq.0) then
        if (iffinite_te) then
          a4 = -0.1015
        else
          a4 = -0.1036
        endif
        tmp = mod(series,100)+0.
        t   = tmp/100.
        tmp = int(series/1000)+0.
        m   = tmp/100.
        tmp = int(series/100.)+0.
        p   = mod(tmp,10.)/10.

        icalld=icalld+1
      endif

      yt = (t/0.2)*(a0*sqrt(xin)+a1*xin+a2*xin**2+a3*xin**3+a4*xin**4)

      if (p.le.1e-12) then
        xout=xin

!       Upper or lower surface point?
        tmp =(yin-yt)**2
        tmp2=(yin+yt)**2
        if (tmp.le.tmp2) then
          yout=yt
        else
          yout=-yt
        endif    

      else
        if (xin.le.p) then
           yc   = (m/p**2)*(2.*p*xin-xin**2) 
           dycdx= (m/p**2)*(2.*p-2.*xin)
        else
           yc   = (m/(1.-p)**2)*((1-2.*p)+2.*p*xin-xin**2) 
           dycdx= (m/(1.-p)**2)*(2.*p-2.*xin)
        endif    

        theta= atan(dycdx)

        xu = xin-yt*sin(theta)
        yu = yc -yt*cos(theta)

        xl = xin+yt*sin(theta)
        yl = yc -yt*cos(theta)

!       Upper or lower surface? 
        tmp  = sqrt((xu-xin)**2 + (yu-yin)**2) 
        tmp2 = sqrt((xl-xin)**2 + (yl-yin)**2)

        if (tmp.le.tmp2) then
          xout=xu
          yout=yu
        else
          xout=xl
          yout=yl
        endif
      endif



      return
      end subroutine fix_naca_bdry
!---------------------------------------------------------------------- 

      subroutine fix_splitter_gll(th)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'

      real dum(3)

      CHARACTER CB*3
      integer Iel,NFACES,npts
      real xnew,ynew

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      real deltax(lx1,ly1,lz1,lelt),deltay(lx1,ly1,lz1,lelt)
      real deltaz(lx1,ly1,lz1,lelt)
      common /psn_scrtch/ deltax,deltay,deltaz

      integer iface,ifld
      integer i,n
      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz
      real chord,rescale
      real foilx_max,foilx_min
      real glmax,glmin

      real th           ! half-thickness

      n=nx1*ny1*nz1*nelt

      call rzero(deltax,lt)
      call rzero(deltay,lt)
      call rzero(deltaz,lt)

      NFACES=2*NDIM

      ifxyo = .true.
      call outpost (vx,vy,vz,pr,t,'gri')
      n = nx1*ny1*nz1*nelv

      ifld = 1
      do  Iel=1,NELV            !do ieg=1,nelgt
        do IFACE = 1,NFACES
          CB = CBC(IFACE,Iel,ifld)
          if (CB.EQ.'mv ') then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
            do iz=KZ1,KZ2
              do iy=KY1,KY2
                do  ix=KX1,KX2
                  deltax(ix,iy,iz,Iel)=xm1(ix,iy,iz,iel)
                  deltay(ix,iy,iz,Iel)=ym1(ix,iy,iz,iel)
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      do  Iel=1,NELV            !do ieg=1,nelgt
        do IFACE = 1,NFACES
          CB = CBC(IFACE,Iel,ifld)
          if (CB.EQ.'mv ') then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
            do iz=KZ1,KZ2
              do iy=KY1,KY2
                do  ix=KX1,KX2
                  call fix_splitter_bdry(th,xm1(ix,iy,iz,iel),
     $                   ym1(ix,iy,iz,iel),xnew,ynew)
                  deltax(ix,iy,iz,Iel)=xnew-xm1(ix,iy,iz,iel)
                  deltay(ix,iy,iz,Iel)=ynew-ym1(ix,iy,iz,iel)
!                    write(6,'(A3,1x,4(E15.7E2,1x))') 'dx:', 
!     $                  deltax(ix,iy,iz,iel),
!     $                  deltay(ix,iy,iz,iel), 
!     $                  xm1(ix,iy,iz,iel),
!     $                  ym1(ix,iy,iz,iel)                     
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      call outpost (deltax,deltay,deltaz,pr,t,'gri')
      call mv_internal_pts(deltax,deltay,deltaz)
      call outpost (deltax,deltay,deltaz,pr,t,'gri')

      return
      end subroutine fix_splitter_gll

!---------------------------------------------------------------------- 

      subroutine fix_splitter_bdry(th,xin,yin,xout,yout)

      implicit none

      real th               ! half-thickness

      real theta
      real xin,yin,xout,yout

      real rad
      parameter (rad=0.5)     ! cylinder radius

      real rad2
      real xtop,ytop
      real xbot,ybot
      real x0,y0
      real xend

!     Find end of cylinder
      ytop = th
      xtop = sqrt(rad**2 - ytop**2)
      ybot = -th
      xbot = xtop

      xend = xtop+1.0

      rad2 = th   ! rounded end radius is same as half-thickness
 
      if (xin.lt.xtop) then
!       Point on the cylinder
        x0=0.
        y0=0.
        theta = atan2(yin,xin)
        xout = x0+rad*cos(theta) 
        yout = y0+rad*sin(theta)
        return
      elseif(xin.gt.xtop.and.xin.lt.xend) then
!       Point of the flat splitter plate
        if (yin.gt.0) yout = th
        if (yin.lt.0) yout =-th
        xout=xin
      else
!       Point is on the rounded end            
        x0=xend
        y0=0.
        theta = atan2(yin-y0,xin-x0)
        xout = x0+rad2*cos(theta)
        yout = y0+rad2*sin(theta)
      endif  

      return
      end subroutine fix_splitter_bdry
!---------------------------------------------------------------------- 

      subroutine mv_internal_pts(mvx,mvy,mvz) ! fix up geometry irregularities
!     Routine has been modified from fix_geom
!     To move internal points of the domain.

      include 'SIZE'
      include 'TOTAL'

      parameter (lt = lx1*ly1*lz1)
      common /scrns/ xb(lt,lelt),yb(lt,lelt),zb(lt,lelt)
      common /scruz/ tmsk(lt,lelt),tmlt(lt,lelt),w1(lt),w2(lt)
      real mvx(lt,lelt),mvy(lt,lelt),mvz(lt,lelt) ! displaced boundary points
      
      integer e,f
      character*3 cb


      n      = nx1*ny1*nz1*nelt
      nxyz   = nx1*ny1*nz1
      nfaces = 2*ndim
      ifield = 1                   ! velocity field
      if (ifheat) ifield = 2       ! temperature field


      call rone  (tmlt,n)
      call dssum (tmlt,nx1,ny1,nz1)  ! denominator

      call rone  (tmsk,n)
      do e=1,nelfld(ifield)      ! fill mask where bc is periodic
      do f=1,nfaces              ! so we don't translate periodic bcs (z only)
         cb =cbc(f,e,ifield)
         if (cb.eq.'P  ') call facev (tmsk,e,f,0.0,nx1,ny1,nz1)
      enddo
      enddo

      do kpass = 1,ndim+1   ! This doesn't work for 2D, yet.
                            ! Extra pass is just to test convergence

c        call opcopy (xb,yb,zb,xm1,ym1,zm1) ! Must use WHOLE field,
c        call opdssum(xb,yb,zb)             ! not just fluid domain.
         call copy   (xb,xm1,n)
         call copy   (yb,ym1,n)
         call copy   (zb,zm1,n)
           
         call dssum  (xb,nx1,ny1,nz1)
         call dssum  (yb,nx1,ny1,nz1)
         call dssum  (zb,nx1,ny1,nz1)

         xm = 0.
         ym = 0.
         zm = 0.

!         if (nid.eq.0) write(6,*) 'IGEOM kpass:', kpass ! prabal

         do e=1,nelt
            do i=1,nxyz                       ! compute averages of geometry
               s     = 1./tmlt(i,e)
!              Prabal
!              If periodicity is in 'Z' then only zm1 needs special treatment.
!              Likewise for other directions      
               xb(i,e) = s*xb(i,e)
               yb(i,e) = s*yb(i,e)
               zb(i,e) = s*zb(i,e)

!              If there are tears between elemnts   
               xb(i,e) = xb(i,e) - xm1(i,1,1,e)   ! local displacements
               yb(i,e) = yb(i,e) - ym1(i,1,1,e)
               zb(i,e) = zb(i,e) - zm1(i,1,1,e)

               xb(i,e) = xb(i,e)*tmsk(i,e)
               yb(i,e) = yb(i,e)*tmsk(i,e)
               zb(i,e) = zb(i,e)*tmsk(i,e)

!              Prabal                  
!              Add in the boundary displacements in the first pass   
               if (kpass.eq.1) then
                  xb(i,e)=xb(i,e)+mvx(i,e)
                  yb(i,e)=yb(i,e)+mvy(i,e)
                  zb(i,e)=zb(i,e)+mvz(i,e)
               endif   


!               write(6,*) xb(i,e),yb(i,e),zb(i,e) ! prabal

               xm = max(xm,abs(xb(i,e)))
               ym = max(ym,abs(yb(i,e)))
               zm = max(zm,abs(zb(i,e)))
            enddo

            if (kpass.le.ndim) then
               call gh_face_extend(xb(1,e),zgm1,nx1,kpass,w1,w2)
               call gh_face_extend(yb(1,e),zgm1,nx1,kpass,w1,w2)
               call gh_face_extend(zb(1,e),zgm1,nx1,kpass,w1,w2)
            endif

         enddo

         if (kpass.le.ndim) then
            call add2(xm1,xb,n)
            call add2(ym1,yb,n)
            call add2(zm1,zb,n)
         endif
        
         xx = glamax(xb,n)
         yx = glamax(yb,n)
         zx = glamax(zb,n)

         xm = glmax(xm,1)
         ym = glmax(ym,1)
         zm = glmax(zm,1)

         if (nio.eq.0) write(6,1) xm,ym,zm,xx,yx,zx,kpass
    1    format(1p6e12.4,' xyz repair',i2)

      enddo

      param(59) = 1.       ! ifdef = .true.
      call geom_reset(1)   ! reset metrics, etc.
      
      return
      end subroutine mv_internal_pts
c-----------------------------------------------------------------------

      subroutine fix_gll_scale(fact,basevinv)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'

      real dum(3)

      CHARACTER CB*3
      integer Iel,NFACES,npts
      real xnew,ynew

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      real deltax(lx1,ly1,lz1,lelt),deltay(lx1,ly1,lz1,lelt)
      real deltaz(lx1,ly1,lz1,lelt)
      common /psn_scrtch/ deltax,deltay,deltaz

      real basevinv(lx1,ly1,lz1,lelt)

      integer iface,ifld
      integer i,n
      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz
      real chord,rescale
      real foilx_max,foilx_min
      real glmax,glmin

      real fact

!      series=0012

      n=nx1*ny1*nz1*nelt

      call rzero(deltax,lt)
      call rzero(deltay,lt)
      call rzero(deltaz,lt)

      NFACES=2*NDIM

      ifxyo = .true.
      call outpost (vx,vy,vz,pr,t,'gri')
      n = nx1*ny1*nz1*nelv

      ifld = 1
      do  Iel=1,NELV            !do ieg=1,nelgt
        do IFACE = 1,NFACES
          CB = CBC(IFACE,Iel,ifld)
          if (CB.EQ.'mv ') then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
            do iz=KZ1,KZ2
              do iy=KY1,KY2
                do  ix=KX1,KX2
                  deltax(ix,iy,iz,Iel)=xm1(ix,iy,iz,iel)
                  deltay(ix,iy,iz,Iel)=ym1(ix,iy,iz,iel)
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      do  Iel=1,NELV            !do ieg=1,nelgt
        do IFACE = 1,NFACES
          CB = CBC(IFACE,Iel,ifld)
          if (CB.EQ.'mv ') then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $           ,IFACE)
            do iz=KZ1,KZ2
              do iy=KY1,KY2
                do  ix=KX1,KX2
                  xnew = xm1(ix,iy,iz,iel)*
     $                   (1+fact*basevinv(ix,iy,iz,iel))
                  ynew = ym1(ix,iy,iz,iel)*
     $                   (1+fact*basevinv(ix,iy,iz,iel))
                  deltax(ix,iy,iz,Iel)=xnew-xm1(ix,iy,iz,iel)
                  deltay(ix,iy,iz,Iel)=ynew-ym1(ix,iy,iz,iel)
!                    write(6,'(A3,1x,4(E15.7E2,1x))') 'dx:', 
!     $                  deltax(ix,iy,iz,iel),
!     $                  deltay(ix,iy,iz,iel), 
!     $                  xm1(ix,iy,iz,iel),
!     $                  ym1(ix,iy,iz,iel)                     
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      call outpost (deltax,deltay,deltaz,pr,t,'gri')
      call mv_internal_pts(deltax,deltay,deltaz)
      call outpost (deltax,deltay,deltaz,pr,t,'gri')

      return
      end subroutine fix_gll_scale

!---------------------------------------------------------------------- 




