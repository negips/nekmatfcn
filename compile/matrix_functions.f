!====================================================================== 
!     Author: Prabal Negi
!     Description: Routines to calculate Matrix functions.
!                : Makes use of the wrapper routines for Lapack.
!
!======================================================================       
!---------------------------------------------------------------------- 
      subroutine MAT_ZFCN(fA,A,lda,nc,fcn)

!     A    = U*T*U^{-1}      ! Schur decomposition             
!     f(A) = U*f(T)*U^{-1}

!     Matrix Function
!     Using the simplified Schur-Parlett Method
!     Might suffer instabilities if the eigenvalues are too close            

      implicit none

      integer lda       ! leading dimension of A
      integer ldu       ! leading dimension of U
      integer nc        ! No of columns (and rows considered)

      complex A(lda,nc),U(lda,nc)
      complex fA(lda,nc) ! f(A). Assuming same shape as A.
      complex w(nc)

      character fcn*4

!     Variables for zgemm 
      complex alpha,beta
      character transA*1,transB*1

!     Complex Schur Decomposition (Double-Precision)
      ldu=lda 
      call wrp_zschur(A,lda,nc,U,ldu,w)
!     A now contains the Complex Upper-Triangular Matrix
!     U has the Schur vectors      

!     Debugging
!      call write_zmat(A,lda,nc,nc,'Tri')
     

!     Calculate Matrix function for a Triangular Matrix
!     Using Schur-Parlett      
      call MAT_ZFCN_TRI(fA,A,lda,nc,fcn)

!     A=fA*U^{H}
      alpha = complex(1.0,0)
      beta  = complex(0.,0.)
      transA = 'N'
      transB = 'C'
      call zgemm(transA,transB,nc,nc,nc,alpha,fA,lda,U,lda,beta,A,lda) 

!     fA=U*A
      alpha = complex(1.0,0)
      beta  = complex(0.,0.)
      transA = 'N'
      transB = 'N'
      call zgemm(transA,transB,nc,nc,nc,alpha,U,lda,A,lda,beta,fA,lda) 


      return
      end subroutine MAT_ZFCN
!----------------------------------------------------------------------

      subroutine MAT_ZFCN_TRI(fA,A,lda,nc,fcn)

!     Apply matrix function to an upper-triangular matrix
!     Using the simplified Schur-Parlett method.
!     Might suffer instabilities if the eigenvalues are too close            

      implicit none

      integer lda       ! leading dimension of A
      integer nc        ! No of columns (and rows considered)

      complex A(lda,nc)
      complex fA(lda,nc) ! f(A). Assuming same shape as A.

      character fcn*4
      complex zfcn      ! complex scalar function
     
      integer i,j,k,p
      complex s

!     Zero complex array
      call nek_zzero(fA,lda*nc)
      
      do i=1,nc
        fA(i,i) = zfcn(A(i,i),fcn)
      enddo        
      
      do p=1,nc-1
      do i=1,nc-p
        j=i+p
        s = A(i,j)*(fA(j,j) - fA(i,i))
        do k=i+1,j-1
          s = s + A(i,k)*fA(k,j)-fA(i,k)*A(k,j)
        enddo  
!       Potential division by zero if A(i,i)=-A(j,j) 
        fA(i,j) = s/(A(j,j)-A(i,i))
      enddo
      enddo      


      return
      end subroutine MAT_ZFCN_TRI
!----------------------------------------------------------------------

      complex function ZFCN(z,fcn)
!     Define all the functions we wish to calculate            

      complex z
      character fcn*4

      if (fcn.eq.'sqrt') then
        zfcn=sqrt(z)
      elseif (fcn.eq.'loge') then
        zfcn=log(z)
      else
        write(6,*) 'Unknown function', fcn
        call exitt 
      endif        
            
      end function ZFCN
!----------------------------------------------------------------------
     
      subroutine nek_zzero(x,n)

      implicit none

      integer n
      complex x(n)
      integer i
      complex z0

      z0 = complex(0.,0.)
      do i=1,n
        x(i)=z0
      enddo  

      return
      end subroutine nek_zzero
!---------------------------------------------------------------------- 
      subroutine nek_zadd2(x,y,n)

      implicit none

      integer n
      complex x(n),y(n)
      integer i

      do i=1,n
        x(i)=x(i)+y(i)
      enddo  

      return
      end subroutine nek_zadd2
!----------------------------------------------------------------------
      subroutine nek_zsub2(x,y,n)

      implicit none

      integer n
      complex x(n),y(n)
      integer i

      do i=1,n
        x(i)=x(i)-y(i)
      enddo  

      return
      end subroutine nek_zsub2
!---------------------------------------------------------------------- 
      subroutine nek_zrcmult(x,a,n)

!     Multiply complex array x with real a

      implicit none

      integer n
      complex x(n)
      real a
      integer i

      do i=1,n
        x(i)=a*x(i)
      enddo  

      return
      end subroutine nek_zrcmult
!---------------------------------------------------------------------- 
      subroutine nek_zcopy(x,y,n)

      implicit none

      integer n
      complex x(n),y(n)
      integer i

      do i=1,n
        x(i)=y(i)
      enddo  

      return
      end subroutine nek_zcopy
!---------------------------------------------------------------------- 
!======================================================================
!     Specialized Algorithms for specific functions
!====================================================================== 
      subroutine MAT_ZFCN_LN(fA,A,lda,nc,pmo)

!     Calculating Matrix logarithm 
!     with Schur decomposition and
!     partial fraction Pade Approximant.

!     A = (I+X)
!                         m            
!     Log(I+X) = r_{m} = Sum    w_{i}*X
!                        i=1  --------------
!                             I + x_{i}*X

!     w_{i} - Weights of m point Gauss-Legendre Quadrature
!     x_{i}  - Nodes of m point Gauss-Legendre Quadrature in [0,1]            

!     A    = U*T*U^{-1}      ! Schur decomposition             
!     f(A) = U*f(T)*U^{-1}

!     Matrix Function

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'WRP.inc'

      integer lda       ! leading dimension of A
      integer ldu       ! leading dimension of U
      integer nc        ! No of columns (and rows considered)

      integer pmo       ! order of Pade approximant

      complex A(lda,nc),U(lda,nc)
      complex fA(lda,nc) ! f(A). Assuming same shape as A.
      complex w(nc)

      complex eye(lkryl1,lkryl1)
      complex X(lkryl1,lkryl1)
      complex WK1(lkryl1,lkryl1)
      complex WK2(lkryl1,lkryl1)

      integer nt

      character fcn*4

!     GL weights and nodes
      real wj(lkryl1)           ! just taking a high enough number
      real xj(lkryl1)           ! for memory

!     Variables for zgemm 
      complex a1,b1
      character transA*1,transB*1

      integer info
      integer i,j

!     Calculate Gauss-Legendre Quadrature
      if (pmo.gt.lkryl1) then
        write(6,*) 'Increase Quadrature memory size for Pade Approx.'
        call exitt
      endif

      call rzero(wj,lkryl1)
      call rzero(xj,lkryl1)
      call ZWGL(xj,wj,pmo)

!     Need the weights and points in [0,1]      
      call cadd(xj,1.,pmo)
      call cmult(xj,0.5,pmo)
      call cmult(wj,0.5,pmo)

!     Complex Identity matrix 
      call nek_zzero(eye,lkryl1*lkryl1)
      do i=1,nc
        eye(i,i)=complex(1.0,0)
      enddo 

!     Complex Schur Decomposition (Double-Precision)
      ldu=lda 
      call wrp_zschur(A,lda,nc,U,ldu,w)
!     A now contains the Complex Upper-Triangular Matrix
!     U has the Schur vectors

!     Debugging
!      call write_zmat(A,lda,nc,nc,'PdT')

      nt = lkryl1*lkryl1
      call nek_zzero(X,nt)
      do i=1,nc
      do j=1,nc
        X(i,j)=A(i,j)-eye(i,j)
      enddo        
      enddo 

!     Evaluate Pade Approximant 
      call nek_zzero(fA,lda*nc)
      call nek_zzero(WK1,nt)
      call nek_zzero(WK2,nt)
      do i=1,pmo
        call nek_zcopy(WK2,X,nt)
        call nek_zrcmult(WK2,xj(i),nt)
        call add2(WK2,eye,nt)

!        call write_zmat(WK2,lkryl1,nc,nc,'WK2')

!       Matrix inversion 
        call ztrtri('U','N',nc,WK2,lkryl1,info)
        if (info.ne.0) then
          write(6,*) 'Matrix inversion unsuccesfull in MAT_FCN_LN'
          call exitt
        endif

!        call write_zmat(WK2,lkryl1,nc,nc,'WKI')


        call nek_zcopy(WK1,X,nt)
        call nek_zrcmult(WK1,wj(i),nt)

!       fA= fA + WK2*WK1
        a1 = complex(1.,0.)
        b1 = complex(1.,0.)
        transA = 'N'
        transB = 'N'
        call zgemm(transA,transB,nc,nc,nc,a1,WK2,lkryl1,
     $             WK1,lkryl1,b1,fA,lda) 

      enddo       

!     A=fA*U^{H}
      a1 = complex(1.,0.)
      b1 = complex(0.,0.)
      transA = 'N'
      transB = 'C'
      call zgemm(transA,transB,nc,nc,nc,a1,fA,lda,U,lda,b1,A,lda) 

!     fA=U*A
      a1 = complex(1.,0.)
      b1 = complex(0.,0.)
      transA = 'N'
      transB = 'N'
      call zgemm(transA,transB,nc,nc,nc,a1,U,lda,A,lda,b1,fA,lda) 


      return
      end subroutine MAT_ZFCN_LN
!----------------------------------------------------------------------










