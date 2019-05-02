!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!!
!!   Author: Prabal Negi
!!   Email : negi@mech.kth.se 
!!   Description: High pass filtering for stabilization of SEM
!!   Last Modified: 26/01/2017 
!!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!!----------------------------------------------------------------------  
!     Main interface for high-pass filter
      subroutine MAKE_HPF

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'TSTEP'         ! ifield
      include 'HPF'

      integer nxyz
      parameter (nxyz=lx1*ly1*lz1)
      integer n

      integer lm,lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)
      real hpf_filter(lm2)

      real op_mat(lx1,lx1)
      save op_mat

      integer nel

      integer icalld
      save icalld
      data icalld /0/

!---------------------------------------- 

      nel = nelfld(ifield)
      n = nxyz*nel

      if (hpf_kai.eq.0) then
        if (ifield.eq.1) then
          call opzero(hpfx,hpfy,hpfz)
        else
          call rzero(hpft(1,1,1,1,ifield-1),n)
        endif

        return
      endif


      if (icalld.eq.0) then

        if (hpf_kai.gt.0) then
          if (nio.eq.0) then
            write(6,*) 'Positive filtering is Numerically Unstable.'
            write(6,*) 'Setting to negative value.'
            write(6,*) 'Remove check in hpf.f if this was intentional.'
          endif    
          hpf_kai = -abs(hpf_kai)
        endif

!       Create the filter transfer function
!       Weighted with hpf_kai
        call hpf_trns_fcn(hpf_filter,hpf_kut,hpf_kai)

!       Build the matrix to apply the filter function
!       to an input tensor field
        call build_hpf_mat(op_mat,hpf_filter,hpf_ifboyd)

!       Only initialize once    
        icalld=icalld+1 
      endif

      if (ifield.eq.1) then
!       Apply the filter
!       to velocity fields
        call build_hpf_fld(hpfx,vx,op_mat,nx1,nz1)
        call build_hpf_fld(hpfy,vy,op_mat,nx1,nz1)
        if (if3d) call build_hpf_fld(hpfz,vz,op_mat,nx1,nz1)

!       Multiply by Mass matrix 
!       and add to forcing term 
        call opadd2col (bfx,bfy,bfz,hpfx,hpfy,hpfz,bm1)

      else

!       Apply filter to temp/passive scalar fields      
        call build_hpf_fld(hpft(1,1,1,1,ifield-1),t(1,1,1,1,ifield-1),
     $       op_mat,nx1,nz1)

!       Multiply by Mass matrix    
!       and add to forcing term
        call addcol3(bq(1,1,1,1,ifield-1),hpft(1,1,1,1,ifield-1),bm1,n)

      endif

      return
      end subroutine MAKE_HPF

!----------------------------------------------------------------------
 
      subroutine addhpf

!     Add high pass filtered field to RHS
!     Must be called after userf.
!     Since userf reinitializes the bfx,bfy,bfz field

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'HPF'

!     Multiply by HPF by Mass matrix
!     and add to RHS forcing term
!     BFX = BFX + HPFX*BM1
      CALL OPADD2COL (BFX,BFY,BFZ,HPFX,HPFY,HPFZ,BM1)

      return
      end subroutine addhpf

c-----------------------------------------------------------------------

      subroutine addhpf_q

!     Add high pass filtered fld to RHS
!     for T and passive scalars
!     Must be called after userq.
!     Since userq reinitializes the bq field

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'         ! ifield
      include 'HPF'

      integer ntot,nel
      
      nel=nelfld(ifield)
      ntot=nx1*ny1*nz1*nel

!     Multiply by Mass matrix
!     Add to RHS forcing term
!     BQ = BQ + HPFT*BM1
      CALL ADDCOL3(BQ(1,1,1,1,ifield-1),HPFT(1,1,1,1,ifield-1),BM1,ntot)

      return
      end subroutine addhpf_q

c-----------------------------------------------------------------------

      subroutine build_hpf_mat(op_mat,f_filter,ifboyd)

!     Builds the operator for high pass filtering
!     Transformation matrix from nodal to modal space.
!     Applies f_filter to the the legendre coefficients
!     Transforms back to nodal space
!     Operation: V * f_filter * V^(-1)
!     Where V is the transformation matrix from modal to nodal space

      implicit none

!      include 'SIZE_DEF'
      include 'SIZE'

      logical IFBOYD 
      integer n
      parameter (n=lx1*lx1)
      integer lm, lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)

      real f_filter(lm2)
      real op_mat(lx1,lx1)

      real ref_xmap(lm2)
      real wk_xmap(lm2)

      real wk1(lm2),wk2(lm2)
      real indr(lm),ipiv(lm),indc(lm)

      real rmult(lm)
      integer ierr

      integer i,j

      call spec_coeff_init(ref_xmap,ifboyd)
      
      call copy(wk_xmap,ref_xmap,lm2)
      call copy(wk1,wk_XMAP,lm2)

      call gaujordf  (wk1,lx1,lx1,indr,indc,ipiv,ierr,rmult)  ! xmap inverse

      call mxm  (f_filter,lx1,wk1,lx1,wk2,lx1)        !          -1
      call mxm  (wk_xmap,lx1,wk2,lx1,op_mat,lx1)      !     V D V

      return
      end subroutine build_hpf_mat

c---------------------------------------------------------------------- 

      subroutine build_hpf_fld(v,u,f,nx,nz)

!     Appies the operator f to field u
!     using tensor operations
!     v = f*u

      implicit none

      include 'SIZE'
      include 'TSTEP'         ! ifield

      integer nxyz 
      parameter (nxyz=lx1*ly1*lz1)

      real w1(nxyz*lelt),w2(nxyz*lelt)    ! work arrays
      real v(nxyz,lelt),u(nxyz,lelt)      ! output and input flds

c
      integer nx,nz

      real f(nx,nx),ft(nx,nx)             ! operator f and its transpose
c
      integer e,i,j,k
      integer nel

      nel = nelfld(ifield)

      call copy(v,u,nxyz*nel)
c
      call transpose(ft,nx,f,nx)
c
      if (ndim.eq.3) then
        do e=1,nel
c         Filter
          call copy(w2,v(1,e),nxyz)
          call mxm(f,nx,w2,nx,w1,nx*nx)
          i=1
          j=1
          do k=1,nx
             call mxm(w1(i),nx,ft,nx,w2(j),nx)
             i = i+nx*nx
             j = j+nx*nx
          enddo
          call mxm (w2,nx*nx,ft,nx,w1,nx)

          call sub3(w2,v(1,e),w1,nxyz)
          call copy(v(1,e),w2,nxyz)
        enddo
      else
        do e=1,nel
c         Filter
          call copy(w1,v(1,e),nxyz)
          call mxm(f ,nx,w1,nx,w2,nx)
          call mxm(w2,nx,ft,nx,w1,nx)

          call sub3(w2,v(1,e),w1,nxyz)
          call copy(v(1,e),w2,nxyz)
        enddo
      endif
c
      return
      end subroutine build_hpf_fld

c---------------------------------------------------------------------- 

      subroutine hpf_trns_fcn(diag,kut,wght)

      implicit none

!      include 'SIZE_DEF'
      include 'SIZE'
!      include 'INPUT_DEF'
      include 'INPUT'
!      include 'PARALLEL_DEF'
      include 'PARALLEL'

      real diag(lx1*lx1)
      real intv(lx1*lx1)
      integer nx,k0,kut,kk,k

      real amp,wght

c     Set up transfer function
c
      nx = lx1
      call ident   (diag,nx)
      call rzero   (intv,nx*nx) 
c

      k0 = nx-kut
      do k=k0+1,nx
        kk = k+nx*(k-1)
        amp = (k-k0)*(k-k0)/(kut*kut)     ! Normalized amplitude. quadratic growth
        diag(kk) = 1.-amp
        intv(kk) = wght*amp               ! weighted amplitude          
      enddo

!     Output normalized transfer function
      k0 = lx1+1
      if (nio.eq.0) then
        write(6,6) 'HPF :',(intv(k)/wght,k=1,lx1*lx1,k0)
   6    format(a8,16f9.6,6(/,8x,16f9.6))
      endif

      return
      end subroutine hpf_trns_fcn

!---------------------------------------------------------------------- 

      subroutine spec_coeff_init(ref_xmap,ifboyd)
!     Initialise spectral coefficients
!     For legendre transform

      implicit none

!      include 'SIZE_DEF'
      include 'SIZE'
!      include 'WZ_DEF'
      include 'WZ'

      integer lm, lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)

!     local variables
      integer i, j, k, n, nx, kj
!     Legendre polynomial
      real plegx(lm)
      real z
      real ref_xmap(lm2)
      real pht(lm2)

!     Change of basis
      logical IFBOYD
!---------------------------------------- 

      nx = LX1
      kj = 0
      n  = nx-1
      do j=1,nx
        z = ZGM1(j,1)
        call legendre_poly(plegx,z,n)
        kj = kj+1
        pht(kj) = plegx(1)
        kj = kj+1
        pht(kj) = plegx(2)

        if (IFBOYD) then        ! change basis to preserve element boundary values
          do k=3,nx
             kj = kj+1
             pht(kj) = plegx(k)-plegx(k-2)
          enddo
        else                    ! legendre basis    
          do k=3,nx
             kj = kj+1
             pht(kj) = plegx(k)
          enddo         
        endif
      enddo

      call transpose (ref_xmap,nx,pht,nx)

      return
      end subroutine spec_coeff_init

!=======================================================================

