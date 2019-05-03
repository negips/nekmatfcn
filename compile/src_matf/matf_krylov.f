!====================================================================== 
!     Author: Prabal Negi
!     Description: Routines to calculate Matrix functions.
!                : Makes use of the wrapper routines for Lapack.
!
!======================================================================       
!----------------------------------------------------------------------

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

!     ZGEMV parameters
      character trans*1
      integer lda       ! leading dimension of Q
      integer m,n       ! Part of the Matrix Q to use m x n
      complex alpha
      complex beta
      integer incx      ! memory skip for x
      integer incy      ! memory skip for y

      complex zdotc     ! BLAS function 


      if (mod(istep,sstep).ne.0) then
        return
      endif

      if (icalld.eq.0) then
        call MATF_INIT
        icalld=icalld+1
        return
      endif

      nsteps=nsteps+1
      lastep = 0

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

      call nek_zcopy(MATF_AxW1,MATF_Ax,vlen)
      call nek_zrcol2(MATF_AxW1,MATF_ArWt,vlen)
     
!     beta = congj(Ax)*ArWt*Ax
      n=vlen
      incx = 1
      incy = 1
      beta = zdotc(n,MATF_Ax,incx,MATF_AxW1,incy)
      beta = sqrt(beta)
      hj(nkryl+1)=beta
!     Update Hessenberg
      call nek_zcopy(MATF_HS(1,nkryl+1),hj,nkryl+1)

!     Normalize vector 
      call nek_zcmult(MATF_Ax,1./beta,vlen)

!     Update Orthogonal matrix
      call nek_zcopy(MATF_Q(1,nkryl+1),MATF_Ax,vlen)

!     Do something to check residual

      call MATF_RESTART      


      return
      end subroutine MATF_MAIN

!---------------------------------------------------------------------- 

      subroutine MATF_INIT

!     Initialize Arnoldi

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'MATFCN'

      integer i

      integer ntot2

      integer n
      complex beta
      integer incx      ! memory skip for x
      integer incy      ! memory skip for y

      complex zdotc     ! BLAS function 


      ntot2=nx2*ny2*nz2*nelv

!     Set weights. 
!     Also sets Arnoldi vector length (vlen)
      call SET_MATF_WEIGHT
!      call SET_MATF_WTONE

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
      beta = zdotc(n,MATF_Ax,incx,MATF_AxW1,incy)   ! BLAS function
      beta = sqrt(beta)
      call nek_zcmult(MATF_Ax,1./beta,vlen)     ! normalize

      if (nio.eq.0) write(6,*) 'Initial Norm=',beta

!     Add starting vector to Krylov space
      call nek_zcopy(MATF_Q(1,1),MATF_Ax,vlen)
      nkryl = 1       ! we don't want to pick the random initial condition

!     Restart stepper
      call MATF_RESTART

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
          if (cb.eq.'v  ') then
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

      subroutine MATF_SETFLDS

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

      integer i,ntot,ntotp
      
      ntot=nx1*ny1*nz1*nelv
      i=1
      call nek_z2ri(vxp(1,1),vxp(1,2),MATF_Ax(i),ntot)
      call col2(vxp(1,1),MATF_MSK(i),ntot)      ! real
      call col2(vxp(1,2),MATF_MSK(i),ntot)      ! imaginary

      i=i+ntot
      call nek_z2ri(vyp(1,1),vyp(1,2),MATF_Ax(i),ntot)
      call col2(vyp(1,1),MATF_MSK(i),ntot)      ! real
      call col2(vyp(1,2),MATF_MSK(i),ntot)      ! imaginary
      i=i+ntot
      if (if3d) then
        call nek_z2ri(vzp(1,1),vzp(1,2),MATF_Ax(i),ntot)           
        call col2(vzp(1,1),MATF_MSK(i),ntot)    ! real
        call col2(vzp(1,2),MATF_MSK(i),ntot)    ! imaginary
        i=i+ntot
      endif

      ntotp=nx2*ny2*nz2*nelv
      if (matf_ifpr) then
        call nek_z2ri(prp(1,1),prp(1,2),MATF_Ax(i),ntot)           
        call col2(prp(1,1),MATF_MSK(i),ntot)    ! real
        call col2(prp(1,2),MATF_MSK(i),ntot)    ! imaginary
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

      subroutine MATF_RESTART

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'         ! ifuzawa
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'MATFCN'

      istep=0
      time=0.

!     set vxp,psi etc from Ax      
      call MATF_SETFLDS

      IFUZAWA = MATF_UZAWA

      return
      end subroutine MATF_RESTART
!---------------------------------------------------------------------- 









