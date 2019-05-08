!====================================================================== 
!     Author: Prabal Negi
!     Description: Routines to calculate Matrix functions.
!                : Makes use of the wrapper routines for Lapack.
!
!======================================================================       
!----------------------------------------------------------------------
!---------------------------------------------------------------------- 
!     read parameters MATRIX FUNCTIONS 
      subroutine matf_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'MATFCN'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /MATFCN/ ifmatf,matf_ifpr,matf_uzawa,ngs,northo,sstep,
     $                  inistep,matf_omega 
!     default values
      ifmatf         = .FALSE.        ! if perform Matrix functions
      matf_uzawa     = .FALSE.        ! uzawa at the first time step?
      matf_ifpr      = .TRUE.         ! if include pressure in
                                      ! arnoldi vector
      northo         = mfnkryl        ! no of vectors to save
      ngs            = 1              ! no of Gram-Schmidt passes
      sstep          = 100            ! No of steps between
                                      ! successive orthogonalization
      inistep        = sstep          ! Skip first inisteps
      matf_omega     = 1.0            ! Angular frequency of forcing
                                         
!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=MATFCN,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading MATFCN parameters.$')

!     broadcast data
      call bcast(ifmatf       , LSIZE)
      call bcast(matf_uzawa   , LSIZE)
      call bcast(matf_ifpr    , LSIZE)
      call bcast(northo       , ISIZE)
      call bcast(ngs          , ISIZE)
      call bcast(sstep        , ISIZE)
      call bcast(inistep      , ISIZE)
      call bcast(matf_omega   ,WDSIZE)

      return
      end subroutine matf_param_in
!----------------------------------------------------------------------
!     write parameters fluid-structure interaction
      subroutine matf_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'MATFCN'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /MATFCN/ ifmatf,matf_ifpr,matf_uzawa,ngs,northo,sstep, 
     $                  inistep,matf_omega

!     read the file
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=MATFCN,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing MATFCN parameters.$')

      return
      end subroutine matf_param_out
c-----------------------------------------------------------------------

      subroutine MATF_MAIN

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'             ! ifuzawa
      include 'MATFCN'

      include 'MASS_DEF'
      include 'MASS'

      integer icalld
      data icalld /0/
      save icalld

      integer i,j,igs
      complex hj(mfnkryl1)
      complex gj(mfnkryl1)
      complex wkh(mfnkryl1)

!     ZGEMV parameters
      character trans*1
      integer lda       ! leading dimension of Q
      integer m,n       ! Part of the Matrix Q to use m x n
      complex alpha
      complex beta
      integer incx      ! memory skip for x
      integer incy      ! memory skip for y

      complex zdotc     ! BLAS function 

      complex MATF_INPROD

      if (.not.IFMATF) return

      if (northo.gt.mfnkryl) then
        if (nid.eq.0) write(6,*)
     $   'NORTHO > MFNKRYL', northo, mfnkryl

        call exitt
      endif        

!     Skip first inisteps      
      if (istep.lt.inistep) return

      if (mod(istep,sstep).ne.0) then
        return
      endif

      if (icalld.eq.0) then
        call MATF_INIT
        icalld=icalld+1
        return
      endif

!      nsteps=nsteps+1
!      lastep = 0

!     fill up Ax,Axtmp vector
      call MATF_getAx

      call nek_zzero(hj,mfnkryl1)
      call nek_zzero(gj,mfnkryl1)

!     Grahm-Schmidt
      do igs=1,ngs            ! No of GS passes

        call nek_zcopy(MATF_AxW1,MATF_Ax,vlen)
!       Multiply by weights
        call nek_zrcol2(MATF_AxW1,MATF_ArWt,vlen)

!       gj = Q^{H}*AxW1
        alpha = complex(1.0,0)
        beta  = complex(0.,0.)
        trans = 'C'     ! Hermitian transpose
        lda=qlen0
        m=vlen
        n=nkryl
        incx = 1
        incy = 1
        call zgemv(trans,m,n,alpha,MATF_Q,lda,MATF_AxW1,
     $                  incx,beta,gj,incy)

!       Sum over processors
        call gop(gj,wkh,'+  ',2*nkryl)

!       Ax = Ax - Q*gj
        alpha = complex(-1.0,0)
        beta  = complex(1.,0.)
        trans = 'N'
        lda=qlen0
        m=vlen
        n=nkryl
        incx = 1
        incy = 1
        call zgemv(trans,m,n,alpha,MATF_Q,lda,gj,
     $                  incx,beta,MATF_Ax,incy)

        call nek_zadd2(hj,gj,nkryl)
      enddo
     
!     Update Hessenberg
!     beta = sqrt(congj(Ax)*ArWt*Ax)
      call nek_zcopy(MATF_AxW1,MATF_Ax,vlen)
      beta = MATF_INPROD(MATF_AxW1,MATF_Ax,MATF_ArWt,vlen)
      beta = sqrt(beta)
      hj(nkryl+1)=beta
      call nek_zcopy(MATF_HS(1,nkryl),hj,nkryl+1)

!     Normalize vector 
      call nek_zcmult(MATF_Ax,1./beta,vlen)

      nkryl= nkryl + 1

!     Update Orthogonal matrix
      call nek_zcopy(MATF_Q(1,nkryl),MATF_Ax,vlen)

!     Check residuals and restart      
      call MATF_LN_RESTART

      return
      end subroutine MATF_MAIN

!---------------------------------------------------------------------- 

      subroutine MATF_INIT

!     Initialize Arnoldi

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'MATFCN'

!      include 'SOLN_DEF'
!      include 'SOLN'

      integer i

      integer ntot2

      integer n
      complex beta
      integer incx      ! memory skip for x
      integer incy      ! memory skip for y

      complex MATF_INPROD     ! BLAS function


      ntot2=nx2*ny2*nz2*nelv

!     Set weights. 
!     Also sets Arnoldi vector length (vlen)
!      call SET_MATF_WEIGHT
      call SET_MATF_WTONE

!     Create Masks
      call SET_MATF_MSK

!     Zero Hessenberg matrix
      i=mfnkryl1*mfnkryl
      call nek_zzero(MATF_HS,i)

!     Zero Qortho matrix
      i=qlen0*mfnkryl1
      call nek_zzero(MATF_Q,i)

!     get starting vector      
      call MATF_getAx

!     Normalize starting vector
      call nek_zcopy(MATF_AxW1,MATF_Ax,vlen)
      call nek_zrcol2(MATF_AxW1,MATF_ArWt,vlen)
     
!     beta = congj(Ax)*ArWt*Ax
      n=vlen
      incx = 1
      incy = 1
      call nek_zcopy(MATF_AxW1,MATF_Ax,vlen)
      beta = MATF_INPROD(MATF_AxW1,MATF_Ax,MATF_ArWt,vlen)
      beta = sqrt(beta)

      call nek_zcmult(MATF_Ax,1./beta,vlen)     ! normalize

      if (nio.eq.0) write(6,*) 'Initial Norm=',abs(beta)

!     Add starting vector to Krylov space
      call nek_zcopy(MATF_Q(1,1),MATF_Ax,vlen)
      nkryl = 1   ! =0 if skip the first arbitrary initial condition

!     This is our forcing field
      call nek_zcopy(MATF_Forc,MATF_Ax,vlen)

!     Add starting vector to Outer Krylov space
      call nek_zcopy(GMR_Q(1,1),MATF_Ax,vlen)
      gmr_nkryl = nkryl

!     remove inistep condition after initialization      
      inistep=0

!     Restart stepper
      call MATF_LN_RESTART

      return
      end subroutine MATF_INIT

!----------------------------------------------------------------------

      subroutine SET_MATF_MSK

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'
      include 'NLFSI'
      include 'MATFCN'  
        
      integer ifld,nfaces,iel,iface
      character cb*3
      integer i,ntot,ntotp

      real mskfld(lx1,ly1,lz1,lelv)

!     Build a mask to remove only 'v  '
!     and preserving the 'mv ' points.
      call rone(mskfld,lx1*ly1*lz1*lelv)
      ifld = 1
      nfaces=2*ndim
      do iel=1,nelv
        do iface = 1,nfaces
          cb = cbc(iface,iel,ifld)
          if ((cb.eq.'v  ').or.(cb.eq.'W  ').or.(cb.eq.'mv ')) then
!           Put zeros on this face
            call facev(mskfld,iel,iface,0.,nx1,ny1,nz1)
          endif
        enddo ! iel
      enddo   ! iface

      ntot=nx1*ny1*nz1*nelv
      i=1
      call copy(MATF_MSK(i),mskfld,ntot)
      i=i+ntot
      call copy(MATF_MSK(i),mskfld,ntot)
      i=i+ntot
      if (if3d) then
        call copy(MATF_MSK(i),mskfld,ntot)
        i=i+ntot
      endif

      if (matf_ifpr) then
        ntotp=nx2*ny2*nz2*nelv
        call rone(MATF_MSK(i),ntotp)
        i=i+ntotp
      endif  

      return
      end subroutine SET_MATF_MSK
!----------------------------------------------------------------------

      subroutine MATF_getAx

!     Put action of the matrix into the vector Ax

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'
      include 'NLFSI'
      include 'MATFCN'

      integer i,ntot,ntotp
      
      ntot=nx1*ny1*nz1*nelv
      i=1
      call nek_ri2z(MATF_Ax(i),vxp(1,1),vxp(1,2),ntot)           
      i=i+ntot
      call nek_ri2z(MATF_Ax(i),vyp(1,1),vyp(1,2),ntot)           
      i=i+ntot
      if (if3d) then
        call nek_ri2z(MATF_Ax(i),vzp(1,1),vzp(1,2),ntot)           
        i=i+ntot
      endif

      if (matf_ifpr) then
        ntotp=nx2*ny2*nz2*nelv
        call nek_ri2z(MATF_Ax(i),prp(1,1),prp(1,2),ntotp)           
        i=i+ntotp
      endif 

!     Multiply by mask      
      call nek_zrcol2(MATF_Ax,MATF_MSK,vlen)

      return
      end subroutine MATF_getAx
!---------------------------------------------------------------------- 

      subroutine MATF_SETFLDS(cA)

!     assign velocity/fsi fields from MAT_Ax entries

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'ADJOINT_DEF'
      include 'ADJOINT'
      include 'FSI'
      include 'NLFSI'
      include 'MATFCN'

      include 'MASS_DEF'
      include 'MASS'

      integer i,ntot,ntotp
      complex cA(1)
      
      ntot=nx1*ny1*nz1*nelv
      i=1
      call nek_z2ri(vxp(1,1),vxp(1,2),cA(i),ntot)
      call col2(vxp(1,1),MATF_MSK(i),ntot)      ! real
      call col2(vxp(1,2),MATF_MSK(i),ntot)      ! imaginary

      i=i+ntot
      call nek_z2ri(vyp(1,1),vyp(1,2),cA(i),ntot)
      call col2(vyp(1,1),MATF_MSK(i),ntot)      ! real
      call col2(vyp(1,2),MATF_MSK(i),ntot)      ! imaginary

      i=i+ntot
      if (if3d) then
        call nek_z2ri(vzp(1,1),vzp(1,2),cA(i),ntot)           
        call col2(vzp(1,1),MATF_MSK(i),ntot)    ! real
        call col2(vzp(1,2),MATF_MSK(i),ntot)    ! imaginary
        i=i+ntot
      endif

      ntotp=nx2*ny2*nz2*nelv
      if (matf_ifpr) then
        call nek_z2ri(prp(1,1),prp(1,2),cA(i),ntotp)           
        call col2(prp(1,1),MATF_MSK(i),ntotp)    ! real
        call col2(prp(1,2),MATF_MSK(i),ntotp)    ! imaginary
        i=i+ntotp
      else
        call rzero(prp(1,1),ntotp)  ! real
        call rzero(prp(1,2),ntotp)  ! imaginary
      endif


      return
      end subroutine MATF_SETFLDS
!---------------------------------------------------------------------- 

      subroutine SET_MATF_WEIGHT

!     Set weights for the Arnoldi vector

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'MASS_DEF'
      include 'MASS'        ! BM1
      include 'INPUT_DEF'
      include 'INPUT'
      include 'MATFCN'

      integer i,ntot,ntotp
      
      ntot=nx1*ny1*nz1*nelv
      i=1
      call copy(MATF_ArWt(i),BM1,ntot)
      i=i+ntot
      call copy(MATF_ArWt(i),BM1,ntot)
      i=i+ntot
      if (if3d) then
        call copy(MATF_ArWt(i),BM1,ntot)
        i=i+ntot
      endif
      
      if (matf_ifpr) then
        ntotp=nx2*ny2*nz2*nelv
!       If we use the semi-norm then
!       weight of pressure is zero        
        call rzero(MATF_ArWt(i),ntotp)
        i=i+ntotp
      endif  

!     Size of krylov vector
      vlen = i-1

      return
      end subroutine SET_MATF_WEIGHT
!---------------------------------------------------------------------- 

      subroutine SET_MATF_WTONE

!     Set weights for the Arnoldi vector to identity

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'MATFCN'

      integer i,ntot,ntotp
      
      ntot=nx1*ny1*nz1*nelv
      i=1
      call rone(MATF_ArWt(i),ntot)
      i=i+ntot
      call rone(MATF_ArWt(i),ntot)
      i=i+ntot
      if (if3d) then
        call rone(MATF_ArWt(i),ntot)
        i=i+ntot
      endif

      if (matf_ifpr) then
        ntotp=nx2*ny2*nz2*nelv
!       If we use the semi-norm then
!       weight of pressure is zero        
        call rzero(MATF_ArWt(i),ntotp)
       
        i=i+ntotp
      endif  

!     Arnoldi Vector length
      vlen = i-1

!     Check vector length      
      if (vlen.gt.qlen0) then
        if (nid.eq.0) then
          write(6,*) 'Inconsistent Arnoldi Vector Length',
     $      vlen,qlen0
        endif
        call exitt
      endif        


      return
      end subroutine SET_MATF_WTONE
!---------------------------------------------------------------------- 

      subroutine MATF_LN_RESTART

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'         ! ifuzawa
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'MATFCN'


      if (nkryl.gt.1) call MATF_EST_RES

!      call MATF_QORTHO_CHK 

      istep=0
      time=0.

!     set vxp,psi etc from Ax
      call MATF_SETFLDS(MATF_Ax)

      IFUZAWA = MATF_UZAWA

      return
      end subroutine MATF_LN_RESTART
!---------------------------------------------------------------------- 

      complex function MATF_INPROD(x,y,wt,n)

      integer n

      complex x(n)
      complex y(n)
      real wt(n)

      integer incx,incy
      complex beta
      complex wk

      call nek_zrcol2(x,wt,n)

!     beta = congj(x)*ArWt*y
      incx = 1
      incy = 1
      beta = zdotc(n,x,incx,y,incy)
!     Sum over processors
      call gop(beta,wk,'+  ',2)

      MATF_INPROD = beta

      end function MATF_INPROD            
!---------------------------------------------------------------------- 

      subroutine MATF_QORTHO_CHK


      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'MATFCN'

      complex matf_inprod
      complex beta

      integer i,j

      if (nkryl.gt.5) return

      if (nid.eq.0) then
        write(6,*) 'MATF_Q: Orthogonality Check'
      endif

      call nek_zzero(MATF_HWK,mfnkryl1,nkryl) 

!     qq = Q^{H}*Q
      do i=1,nkryl
        do j=1,nkryl
          call nek_zcopy(MATF_AxW1,MATF_Q(1,i),vlen)
          MATF_HWK(i,j) = MATF_INPROD(MATF_AxW1,MATF_Q(1,j),
     $                        MATF_ArWt,vlen)
        enddo          
      enddo

      if (nid.eq.0) then      
        call write_zmat(MATF_HWK,mfnkryl1,nkryl,nkryl,'QHQ')
      endif        


      return
      end subroutine MATF_QORTHO_CHK
!---------------------------------------------------------------------- 

      subroutine MATF_EST_RES

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'MATFCN'

!     ZGEMV parameters
      character trans*1
      integer lda       ! leading dimension of Q
      integer m,n       ! Part of the Matrix Q to use m x n
      complex alpha
      complex beta
      integer incx      ! memory skip for x
      integer incy      ! memory skip for y

      complex zdotc     ! BLAS function

!     Matrix function evaluation      
      integer nc
      character fcn*4
      logical ifinv
      integer pmo       ! Pade approximant order

      complex MATF_INPROD

      complex MATF_HINV(mfnkryl1,mfnkryl) ! Inverse Hessenberg Matrix
      complex HS_WK2(mfnkryl1,mfnkryl)    ! Another work array

      real resid
      real lntol
      parameter (lntol=1.0e-6)


!     Evaluate approximate matrix function
      lda = mfnkryl1
      nc  = nkryl-1
      fcn = 'loge'
      ifinv = .false.
      pmo = 16
      call nek_zcopy(MATF_HINV,MATF_HS,lda*nc)

      call MAT_ZFCN_LN(MATF_HWK,MATF_HINV,HS_WK2,lda,nc,pmo,ifinv)
!      call MAT_ZFCN(MATF_HWK,MATF_HINV,HS_WK2,lda,nc,fcn,ifinv)

!      call write_zmat(MATF_HWK,mfnkryl1,nc,nc,'Hes')      

      resid = abs(MATF_HWK(nc,1))

      if (nid.eq.0) then
!       This is the weight of the new vector
        write(6,*) 'Qk Weight',nc,resid
      endif

      if (nkryl.eq.(northo+1)) then
        if (nid.eq.0) write(6,*) 'f(A)x unconverged',resid
        resid=0.
      endif        


!     If the weight is small we are converged      
      if (resid.lt.lntol) then
!!      Approximate f(A)*x      
        alpha = complex(1.0,0)
        beta  = complex(0.,0.)
        trans = 'N'
        lda=qlen0
        m=vlen
        n=nkryl-1
        incx = 1
        incy = 1
!       Overwrite Ax
!       Ax = Q*f(H)*||b|| 
        call zgemv(trans,m,n,alpha,MATF_Q,lda,MATF_HWK,
     $                  incx,beta,MATF_Ax,incy)

!       Generate new vector for Outer Krylov loop 
        call GMR_EXTEND_KRYLOV

      endif        


      return
      end subroutine MATF_EST_RES            
!---------------------------------------------------------------------- 

      subroutine GMR_EXTEND_KRYLOV

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'MATFCN'

      integer i,j,igs
      complex hj(gmrkryl1)
      complex gj(gmrkryl1)
      complex wkh(gmrkryl1)

!     ZGEMV parameters
      character trans*1
      integer lda       ! leading dimension of Q
      integer m,n       ! Part of the Matrix Q to use m x n
      complex alpha
      complex beta
      integer incx      ! memory skip for x
      integer incy      ! memory skip for y

      complex zdotc     ! BLAS function

!     Matrix function evaluation      
      integer nc
      character fcn*4
      logical ifinv
      integer pmo       ! Pade approximant order

      complex MATF_INPROD


      call nek_zzero(hj,mfnkryl1)
      call nek_zzero(gj,mfnkryl1)

!     Grahm-Schmidt
      do igs=1,ngs            ! No of GS passes

        call nek_zcopy(MATF_AxW1,MATF_Ax,vlen)
!       Multiply by weights
        call nek_zrcol2(MATF_AxW1,MATF_ArWt,vlen)

!       gj = Q^{H}*AxW1
        alpha = complex(1.0,0)
        beta  = complex(0.,0.)
        trans = 'C'     ! Hermitian transpose
        lda=qlen0
        m=vlen
        n=gmr_nkryl
        incx = 1
        incy = 1
        call zgemv(trans,m,n,alpha,GMR_Q,lda,MATF_AxW1,
     $                  incx,beta,gj,incy)

!       Sum over processors
        call gop(gj,wkh,'+  ',2*gmr_nkryl)

!       Ax = Ax - Q*gj
        alpha = complex(-1.0,0)
        beta  = complex(1.,0.)
        trans = 'N'
        lda=qlen0
        m=vlen
        n=gmr_nkryl
        incx = 1
        incy = 1
        call zgemv(trans,m,n,alpha,GMR_Q,lda,gj,
     $                  incx,beta,MATF_Ax,incy)

        call nek_zadd2(hj,gj,gmr_nkryl)
      enddo

!     Update Hessenberg
!     beta = sqrt(congj(Ax)*ArWt*Ax)
      call nek_zcopy(MATF_AxW1,MATF_Ax,vlen)
      beta = MATF_INPROD(MATF_AxW1,MATF_Ax,MATF_ArWt,vlen)
      beta = sqrt(beta)
      hj(gmr_nkryl+1)=beta
      call nek_zcopy(GMR_HS(1,gmr_nkryl),hj,gmr_nkryl+1)

!     Debugging      
!      call write_zmat(GMR_HS,gmrkryl1,gmr_nkryl+1,gmr_nkryl,'GMR')

!     Normalize vector 
      call nek_zcmult(MATF_Ax,1./beta,vlen)

      gmr_nkryl = gmr_nkryl + 1

!     Update Outer Orthogonal matrix
      call nek_zcopy(GMR_Q(1,gmr_nkryl),MATF_Ax,vlen)

!     Reinitialize Orthogonal basis for Matrix log
      nkryl = 1
      call nek_zcopy(MATF_Q(1,1),MATF_Ax,vlen)

!     Just exit for now      
      if (gmr_nkryl.eq.gmrkryl1) call exitt


      

      return
      end subroutine GMR_EXTEND_KRYLOV            
!---------------------------------------------------------------------- 

















