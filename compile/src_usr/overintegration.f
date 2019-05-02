!---------------------------------------------------------------------- 
!     Author: Prabal Negi
!     Description: New routine for overintegration
!---------------------------------------------------------------------- 

      subroutine set_convect_prabal(cr,cs,ct,ux,uy,uz)
C
C     Put vxd,vyd,vzd into rst form on fine mesh
C
C     For rst form, see eq. (4.8.5) in Deville, Fischer, Mund (2002).
C
      include 'SIZE'
      include 'TOTAL'

      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

      real cr(ltd,1),cs(ltd,1),ct(ltd,1)
      real ux(lxy,1),uy(lxy,1),uz(lxy,1)

      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)

      integer e

      call set_dealias_rx

      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd

      ic = 1    ! pointer to vector field C

      do e=1,nelv 

c        Map coarse velocity to fine mesh (C-->F)

         call intp_rstd(cr,ux(1,e),nx1,nxd,if3d,0) ! 0 --> forward
         call intp_rstd(cs,uy(1,e),nx1,nxd,if3d,0) ! 0 --> forward
         if (if3d) call intp_rstd(ct,uz(1,e),nx1,nxd,if3d,0) ! 0 --> forward

      enddo

      return
      end


!---------------------------------------------------------------------- 

      subroutine Calc_convection(cdu,fld)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'INPUT_DEF'
      include 'INPUT'
 
      integer lt
      parameter (lt=lx1*ly1*lz1)
      integer ltd
      parameter (ltd=lxd*lyd*lzd)
      real cxd(ltd,lelv),cyd(ltd,lelv),czd(ltd,lelv)  ! convection on fine grid
      real fxd(ltd,lelv),fyd(ltd,lelv),fzd(ltd,lelv)  ! field derivate on fine grid
      real wkd(ltd,lelv)
     
      real cfld(lt,lelv,3)
      real fld(lt,lelv)
      real dfld(lt,lelv)
      real cdu(lt,lelv)
      real wk(lt,lelv)

      integer ntotd
      integer e

      integer IMESH

      REAL JACM1D(LXD,LYD,LZD,LELV)
      REAL W3MD(LXD,LYD,LZD)
      REAL OVERINTG_MAT(LX1*LY1*LZ1,LXD*LYD*LZD)
      common /OVERINT/ JACM1D,W3MD,OVERINTG_MAT


      if (nid.eq.0) then
        write(6,*) 'Alternate overintegration'
      endif

      ntotd = ltd*nelv
      IMESH = 1               ! velocity mesh

!     Put convecting field in over-integration mesh
      call opcopy(cfld(1,1,1),cfld(1,1,2),cfld(1,1,3),vx,vy,vz)
      call mapw   (cxd ,nxd,cfld(1,1,1) ,nx1,1)
      call mapw   (cyd ,nxd,cfld(1,1,2) ,nx1,1)
      if (if3d) call mapw   (czd ,nxd,cfld(1,1,3) ,nx1,1)

!     Initialize output array
      call rzero(wkd,ltd*nelv)

!     Cx.dfdx
      CALL DUDXYZ (dfld,fld,RXM1,SXM1,TXM1,JACM1,IMESH,1)
      call mapw   (fxd,nxd,dfld,nx1,1)
      CALL COL3   (wk,cxd,fxd,ntotd)
      call add2   (wkd,wk,ntotd)

!     Cy.dfdy
      CALL DUDXYZ (dfld,fld,RYM1,SYM1,TYM1,JACM1,IMESH,2)
      call mapw   (fyd,nxd,dfld,nx1,1)
      CALL COL3   (wk,cyd,fyd,ntotd)
      call add2   (wkd,wk,ntotd)

      if (if3d) then
!       Cz.dfdz
        CALL DUDXYZ (dfld,fld,RZM1,SZM1,TZM1,JACM1,IMESH,3)
        call mapw   (fzd,nxd,dfld,nx1,1)
        CALL COL3   (wk,czd,fzd,ntotd)
        call add2   (wkd,wk,ntotd)
      endif

      call initialize_overintegration

      call col2(wkd,JACM1D,ltd*nelv)

      do e=1,nelv
        call mxm(OVERINTG_MAT,lx1*ly1*lz1,wkd(1,e),ltd,cdu(1,e),1)
      enddo


      return

      end subroutine Calc_convection

!----------------------------------------------------------------------

      subroutine Calc_convection2(cdu,fld)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'INPUT_DEF'
      include 'INPUT'
 
      integer lt
      parameter (lt=lx1*ly1*lz1)
      integer ltd
      parameter (ltd=lxd*lyd*lzd)
      real cxd(ltd,lelv),cyd(ltd,lelv),czd(ltd,lelv)  ! convection on fine grid
      real fxd(ltd,lelv),fyd(ltd,lelv),fzd(ltd,lelv)  ! field derivate on fine grid
      real fx(lt,lelv),fy(lt,lelv),fz(lt,lelv)  ! field derivate on fine grid

      real wkd(ltd,lelv)
     
      real cfld(lt,lelv,3)
      real fld(lt,lelv)             !
      real dfld(lt,lelv)
      real cdu(lt,lelv)
      real wk(lt,lelv)

      integer ntotd
      integer e

      integer IMESH

      REAL JACM1D(LXD,LYD,LZD,LELV)
      REAL W3MD(LXD,LYD,LZD)
      REAL OVERINTG_MAT(LX1*LY1*LZ1,LXD*LYD*LZD)
      common /OVERINT/ JACM1D,W3MD,OVERINTG_MAT


      ntotd = ltd*nelv
      IMESH = 1               ! velocity mesh

      call opcopy(cfld(1,1,1),cfld(1,1,2),cfld(1,1,3),vx,vy,vz)

      do e=1,nelv

        call intp_rstd(cxd(1,e),cfld(1,e,1),nx1,nxd,if3d,0) ! 0 --> interpolate to fine grid
        call intp_rstd(cyd(1,e),cfld(1,e,2),nx1,nxd,if3d,0) ! 0 --> interpolate to fine grid
        if (if3d) then
          call intp_rstd(czd(1,e),cfld(1,e,ndim),nx1,nxd,if3d,0) ! 0 --> interpolate to fine grid
        endif

!        do i=1,nrstd
!           tu(i)=c(i,e,j)*ju(i)   ! C_j*T
!        enddo

!        call intp_rstd(wkd(1,e),fld,nx1,nxd,if3d,0) ! 0 --> forward
!        call grad_rst(fxd(1,e),fyd(1,e),fzd(1,e),wkd,nxd,if3d)    ! Get gradients on Gauss mesh 

!        call grad_rstd(fxd(1,e),fyd(1,e),fzd(1,e),
!     $                         fld(1,e),lx1,lxd,if3d,wkd)  ! gradients on the GL mesh


!       gradients of the field.
        call gradm1(fx,fy,fz,fld)
!       interpolate gradients to fine grid.
        call intp_rstd(fxd(1,e),fx(1,e),nx1,nxd,if3d,0) ! 0 --> interpolate to fine grid
        call intp_rstd(fyd(1,e),fy(1,e),nx1,nxd,if3d,0) ! 0 --> interpolate to fine grid
        if (if3d) then
          call intp_rstd(fzd(1,e),fz(1,e),nx1,nxd,if3d,0) ! 0 --> interpolate to fine grid
        endif

      enddo

      call col2(fxd,cxd,ltd*nelv)
      call col2(fyd,cyd,ltd*nelv)
      call add2(fxd,fyd,ltd*nelv) 
      if (if3d) then
        call col2(fzd,czd,ltd*nelv) 
        call add2(fxd,fzd,ltd*nelv)
      endif

      call initialize_overintegration

      call col2(fxd,JACM1D,ltd*nelv)

      do e=1,nelv
        call mxm(OVERINTG_MAT,lx1*ly1*lz1,fxd(1,e),ltd,cdu(1,e),1)
      enddo


      return

      end subroutine Calc_convection2

!----------------------------------------------------------------------

      subroutine initialize_overintegration

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'DEALIAS_DEF'
      include 'DEALIAS'
      include 'INPUT_DEF'
      include 'INPUT'      

      REAL JACM1D(LXD,LYD,LZD,LELV)
      REAL W3MD(LXD,LYD,LZD)
      REAL OVERINTG_MAT(LX1*LY1*LZ1,LXD*LYD*LZD)
      common /OVERINT/ JACM1D,W3MD,OVERINTG_MAT

      real wxmd(lxd),wymd(lxd),wzmd(lxd)
      real zgmd(lxd,3)
      real wk(lx1*ly1*lz1*lxd*lyd*lzd)
      real legtogauss(lxd,lx1)
      real trialfcn(lx1*ly1*lz1)
      real trialfcn_GL(lxd*lyd*lzd)

      integer ix,iy,iz
      integer e

      integer lt
      parameter (lt=lx1*ly1*lz1)

      integer ltd
      parameter (ltd=lxd*lyd*lzd)

      integer icalld
      save icalld
      data icalld /0/


!     Jacobian is re-calculated since elements can be deformed.
!     For example for moving mesh simulations.

!     Interpolate Jacobian to over-integration mesh.
!     Not entirely sure if it should be interpolated
!     Or a new mapping in higher order should be calculated.
!     Fow low order mappings it should not matter.
      call mapw   (JACM1D,nxd,JACM1,nx1,1)

!     The rest of the variables are on reference element.
!     They only need to be calculated once.
      if (icalld.gt.0) return

      icalld=icalld+1

C     Generate over-integration points and weights
      if (ndim.eq.2) then
        call zwgl (zgmd(1,1),wxmd,nxd)
        call zwgl (zgmd(1,2),wymd,nyd)
        zgmd(nzd,3) = 0.
        wzmd(nzd)   = 1.
        do 100 iy=1,nyd
          do 100 ix=1,nxd
            W3MD(IX,IY,1)=wxmd(ix)*wymd(iy)
  100   continue

      else
        call zwgl (zgmd(1,1),wxmd,nxd)
        call zwgl (zgmd(1,2),wymd,nyd)
        call zwgl (zgmd(1,3),wzmd,nzd)
        do 200 iy=1,nyd
          do 200 ix=1,nxd
            do 200 iz=1,nzd
              W3MD(IX,IY,IZ)=wxmd(ix)*wymd(iy)*wzmd(iz)
  200   continue
      endif     


!     Build an over-integration matrix.

      call rzero(trialfcn,lt)
      e=0
      do 300 iz=1,lz1
      do 300 iy=1,ly1
      do 300 ix=1,lx1
        if (e>0) then
          trialfcn(e)=0
        endif
        e=e+1
        trialfcn(e)=1
        call intp_rstd(trialfcn_GL,trialfcn,nx1,nxd,if3d,0)  ! interpolate trial function
        call col2(trialfcn_GL,W3MD,ltd)                      ! multiply by weights
        call copy(wk((e-1)*ltd+1),trialfcn_GL,ltd)
  300 continue

      call transpose(OVERINTG_MAT,lt,wk,ltd) 

      return
      end subroutine initialize_overintegration

!----------------------------------------------------------------------


c-----------------------------------------------------------------------

      subroutine build_leg_lxd(legtogauss)

!     Build legendre interpolation to over-integration grid

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'

      integer lm, lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)

!     local variables
      integer i, j, k, n, nx, kj
!     Legendre polynomial
      real plegx(lm)
      real z
      real pht(lm2)

      real xm1d(lxd)
      real wts_xm1d(lxd)
      real legtogauss(lxd,lx1)
      real wk(lxd*lxd)

      real legtrnsf(lx1*lx1),invlegtrnsf(lx1*lx1)

      call zwgl  (xm1d,wts_xm1d,lxd)

      nx = LX1
      kj = 0
      n  = nx-1
      do j=1,lxd
        z = xm1d(j)
        call legendre_poly(plegx,z,n)
        do k=1,lx1
           kj = kj+1
           pht(kj) = plegx(k)
        enddo         
      enddo

      call transpose (wk,lxd,pht,lx1)
!      call copy(legtogauss,pht,lxd*lx1)

      call build_leg_transform(legtrnsf,invlegtrnsf)
      call mxm(wk,lxd,legtrnsf,lx1,legtogauss,lx1)

      return
      end subroutine build_leg_lxd

!----------------------------------------------------------------------       
      subroutine build_leg_transform(legtrnsf,invlegtrnsf)

!     Initialise spectral coefficient mapping
!     Modified by Prabal for relaxation term implementation.

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'WZ'
      include 'WZ_DEF'

      integer lm, lm2
      parameter (lm=lx1)
      parameter (lm2=lm*lm)

!     local variables
      integer i, j, k, n, nx, kj
!     Legendre polynomial
      real plegx(lm)
      real z
      real pht(lm2)
      real phi(lm2)
      real wk(lm2)
      real invlegtrnsf(lx1,lx1)
      real legtrnsf(lx1,lx1)
      real indr(lm),ipiv(lm),indc(lm),rmult(lm)
      integer ierr

      nx = LX1
      kj = 0
      n  = nx-1
      do j=1,lx1
        z = ZGM1(j,1)
        call legendre_poly(plegx,z,n)
        do k=1,lx1
          kj = kj+1
          pht(kj) = plegx(k)
        enddo         
      enddo

      call transpose (invlegtrnsf,lx1,pht,lx1)
      call copy(legtrnsf,invlegtrnsf,lx1*lx1)

      call gaujordf(legtrnsf,lx1,lx1,indr,indc,ipiv,ierr,rmult)  ! gauss jordan inverse

      return
      end subroutine build_leg_transform

!---------------------------------------------------------------------- 

 
