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










