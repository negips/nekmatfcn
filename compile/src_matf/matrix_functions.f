!====================================================================== 
!     Author: Prabal Negi
!     Description: Routines to calculate Matrix functions.
!                : Makes use of the wrapper routines for Lapack.
!     Modified   : 04-05-2019
!
!======================================================================       
!---------------------------------------------------------------------- 
      subroutine MAT_ZFCN(fA,A,wkA,lda,nc,fcn,ifinv)

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
      complex wkA(lda,nc)     ! Work array
      complex w(nc)

      character fcn*4

      logical ifinv

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

!     Multiply by Schur vectors      
      call SCHURVEC_MULT(fA,U,A,wkA,lda,nc,ifinv)


      return
      end subroutine MAT_ZFCN
!----------------------------------------------------------------------

      subroutine MAT_ZFCN_TRI(fT,T,lda,nc,fcn)

!     fT = f(T)

!     Apply matrix function to an upper-triangular matrix
!     Using the simplified Schur-Parlett method.
!     Might suffer instabilities if the eigenvalues are too close            

      implicit none

      integer lda       ! leading dimension of T
      integer nc        ! No of columns (and rows considered)

      complex T(lda,nc)
      complex fT(lda,nc) ! f(T). Assuming same shape as T.

      character fcn*4
      complex zfcn      ! complex scalar function
     
      integer i,j,k,p
      complex s

!     Zero complex array
      call nek_zzero(fT,lda*nc)
      
      do i=1,nc
        fT(i,i) = zfcn(T(i,i),fcn)
      enddo        
      
      do p=1,nc-1
      do i=1,nc-p
        j=i+p
        s = T(i,j)*(fT(j,j) - fT(i,i))
        do k=i+1,j-1
          s = s + T(i,k)*fT(k,j)-fT(i,k)*T(k,j)
        enddo  
!       Potential division by zero if T(i,i)=-T(j,j) 
        fT(i,j) = s/(T(j,j)-T(i,i))
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
      subroutine SCHURVEC_MULT(fA,U,wk1,wk2,lda,nc,ifinv)

      implicit none

      integer lda,nc

      complex fA(lda,nc) ! f(A). Assuming same shape as A.
      complex U(lda,nc)
      complex wk1(lda,nc)     ! Work array
      complex wk2(lda,nc)

!     Variables for zgemm 
      complex a1,b1
      character transA*1,transB*1
      integer m,n,k

      integer info
      integer i,j

      logical ifinv

!     Copy solution      
      call nek_zcopy(wk1,fA,lda*nc)

!     wk2=fA*U^{H}
      a1 = complex(1.,0.)
      b1 = complex(0.,0.)
      transA = 'N'
      transB = 'C'
      m=nc
      n=nc
      k=nc
      call zgemm(transA,transB,m,n,k,a1,fA,lda,U,lda,b1,wk2,lda) 

!     fA=U*wk2
      a1 = complex(1.,0.)
      b1 = complex(0.,0.)
      transA = 'N'
      transB = 'N'
      m=nc
      n=nc
      k=nc
      call zgemm(transA,transB,m,n,k,a1,U,lda,wk2,lda,b1,fA,lda)

!     Calculate inverse if needed      
      if (ifinv) then      
!       Matrix inversion of a Triangular Matrix
        call ztrtri('U','N',nc,wk1,lda,info)
        if (info.ne.0) then
          write(6,*) 'Matrix inversion unsuccesfull in SCHURVEC_MULT'
          call exitt
        endif
!       wk1 now contains f(T)^{-1}

!       wk2=wk1*U^{H}
        a1 = complex(1.0,0)
        b1  = complex(0.,0.)
        transA = 'N'
        transB = 'C'
        m=nc
        n=nc
        k=nc
        call zgemm(transA,transB,m,n,k,a1,wk1,lda,U,lda,b1,wk2,lda) 

!       wk1=U*wk2
        a1 = complex(1.0,0)
        b1 = complex(0.,0.)
        transA = 'N'
        transB = 'N'
        m=nc
        n=nc
        k=nc
        call zgemm(transA,transB,m,n,k,a1,U,lda,wk2,lda,b1,wk1,lda)
!       wk1 now contains f(A)^{-1}      
      endif

      return            
      end subroutine SCHURVEC_MULT            
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
      subroutine nek_zrcol2(x,y,n)

      implicit none

      integer n
      complex x(n)
      real y(n)
      integer i

      do i=1,n
        x(i)=x(i)*y(i)
      enddo  

      return
      end subroutine nek_zrcol2
!----------------------------------------------------------------------

      subroutine nek_zcmult(x,z,n)

!     Multiply complex array x with a complex z

      implicit none

      integer n
      complex x(n)
      complex z
      integer i

      do i=1,n
        x(i)=z*x(i)
      enddo  

      return
      end subroutine nek_zcmult
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
      subroutine nek_ri2z(z,x,y,n)

      implicit none

      integer n
      complex z(n)
      real x(n),y(n)
      integer i

      do i=1,n
        z(i)=complex(x(i),y(i))
      enddo  

      return
      end subroutine nek_ri2z
!---------------------------------------------------------------------- 
      subroutine nek_r2z(z,x,n)

      implicit none

      integer n
      complex z(n)
      real x(n)
      integer i

      do i=1,n
        z(i)=complex(x(i),0.)
      enddo  

      return
      end subroutine nek_r2z
!---------------------------------------------------------------------- 
      subroutine nek_i2z(z,x,n)

      implicit none

      integer n
      complex z(n)
      real x(n)
      integer i

      do i=1,n
        z(i)=complex(0.,x(i))
      enddo  

      return
      end subroutine nek_i2z
!---------------------------------------------------------------------- 

      subroutine nek_z2ri(x,y,z,n)

      implicit none

      integer n
      complex z(n)
      real x(n),y(n)
      integer i

      do i=1,n
        x(i)=real(z(i))
        y(i)=imag(z(i))
      enddo  

      return
      end subroutine nek_z2ri
!----------------------------------------------------------------------
!======================================================================
!     Specialized Algorithms for some functions
!     Slowly growing this list      
!====================================================================== 
      subroutine MAT_ZFCN_LN(fA,A,wkA,lda,nc,pmo,ifinv)

!     A    = U*T*U^{-1}      ! Schur decomposition             
!     f(A) = U*f(T)*U^{-1}

!     Calculating Matrix logarithm 
!     with Schur decomposition and
!     partial fraction Pade Approximant.

!     Algorithm 11.9, pp 277-278 in Higham (2008) Functions of Matrices
!     - Theory and Computation
           
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'MFN.inc'

      integer lda       ! leading dimension of A
      integer ldu       ! leading dimension of U
      integer nc        ! No of columns (and rows considered)

      integer pmo       ! order of Pade approximant
      integer pm
      logical ifinv     ! If get inverse (fA)^{-1}

      complex A(lda,nc),U(lda,nc)
      complex fA(lda,nc) ! f(A). Assuming same shape as A.
      complex w(nc)
      complex wkA(lda,nc)     ! Work array
      complex wk2(mfnkryl1,mfnkryl1)     ! Work array

      complex eye(mfnkryl1,mfnkryl1)      ! identity

      integer i

      integer maxroots
      integer nroots          ! No of square roots taken
      real scal               ! 2^{nroots}

      real zlange       ! LAPACK. 1 norm
      real tau
      character norm*1
      real wkNorm(mfnkryl1)   ! not referenced for 1 norm

      character fcn*4

      real theta(16)

!     These are the reference values from the book Higham (2008).
!     They are based on the forward error propagation analysis and
!     matrix norm ||X||
      parameter (theta = [1.10e-5, 1.82e-3, 1.62e-2, 5.39e-2,1.14e-1,
     $                      1.87e-1, 2.64e-1, 3.40e-1, 4.11e-1, 4.75e-1,
     $                      5.31e-1, 5.81e-1, 6.24e-1, 6.62e-1, 6.95e-1,
     $                      7.24e-1])

!!     These are the reference values from the paper 
!      Al-Mohy & Higham (2012) Improved Inverse Scaling and Squaring 
!      Algorithms for the Matrix Logarithm. They are based on the 
!      backward error propagation analysis and the norm ||X^{p}||^{1/p}
!      parameter (theta= [1.59e-5, 2.31e-3, 1.94e-2, 6.21e-2,
!     $                   1.28e-1, 2.06e-1, 2.88e-1, 3.67e-1, 
!     $                   4.39e-1, 5.03e-1, 5.60e-1, 6.09e-1,
!     $                   6.52e-1, 6.89e-1, 7.21e-1, 7.49e-1])
       

      if (pmo.gt.16) then
        if (nid.eq.0) write(6,*) 'Setting Pade Order to 16'
        pmo=16
      endif        

!     Complex Schur Decomposition (Double-Precision)
      ldu=lda 
      call wrp_zschur(A,lda,nc,U,ldu,w)
!     A now contains the Complex Upper-Triangular Matrix
!     U has the Schur vectors

!!     Debugging      
!      call write_zmat(A,lda,nc,nc,'Tri')

!     make a copy      
      call nek_zcopy(wkA,A,lda*nc)

!     Complex Identity matrix 
      call nek_zzero(eye,mfnkryl1*mfnkryl1)
      do i=1,nc
        eye(i,i) = complex(1.0,0)
      enddo 

      do i=1,nc
        wkA(i,i)   = A(i,i)-eye(i,i)
      enddo

!!     Debugging      
!      call write_zmat(wkA,lda,nc,nc,'wkA')


!     Find the ||A-I||_{1}
      norm = 'F'
      tau  = zlange(norm,nc,nc,wkA,lda,wkNorm)

      if (nid.eq.0) then
        write(6,'(A16,1x,E25.16E2)') 'Matrix Log norm:',tau
!        write(6,*) theta
      endif        

!     If the norm is already small
      pm = pmo 
      do i=1,pmo
        if (tau.lt.theta(i)) then
          pm=i
          exit
        endif
      enddo 

!     Matrix Square roots
      nroots = 0
      norm   = '1'
      maxroots = 30
      do while (tau.ge.theta(pm).and.nroots.le.maxroots)

        call TRIMAT_SQRT(fA,A,lda,nc)
        nroots = nroots+1
        call nek_zcopy(A,fA,lda*nc)
        call nek_zcopy(wkA,A,lda*nc)
        do i=1,nc
          wkA(i,i) = A(i,i)-eye(i,i)
        enddo

!       Find the ||A-I||_{1}
        tau  = zlange(norm,nc,nc,wkA,lda,wkNorm)

      enddo
      scal = 2.**nroots
      if (nid.eq.0) then
        write(6,'(A22,1x,I3,1x,2(E25.16E2,1x))') 'k,||T^{1/k}||,theta_k'
     $    ,nroots,tau,theta(pm)
      endif        

!     Evaluate Pade Approximant
      call trimat_ln_pade(fA,A,lda,nc,pmo)

!     Multiply by Schur vectors
!     fA = U*fA*U^{H}
!     ifinv=.true. : A = U*A^{-1}*U^{H}
      call SCHURVEC_MULT(fA,U,A,wkA,lda,nc,ifinv)

!     Multiply by scalar factor      
      call nek_zrcmult(fA,scal,lda*nc)
     

      return
      end subroutine MAT_ZFCN_LN
!----------------------------------------------------------------------
      subroutine TRIMAT_LN_PADE(fT,T,lda,nc,pmo)

!     T - Triangular Matrix after Schur decomposition
!     fT = f(T)            

!     T  = (I+X)
!                         m            
!     Log(I+X) = r_{m} = Sum    w_{i}*X
!                        i=1  --------------
!                             I + x_{i}*X

!     w_{i} - Weights of m point Gauss-Legendre Quadrature
!     x_{i}  - Nodes of m point Gauss-Legendre Quadrature in [0,1]            

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'MFN.inc'

      integer lda       ! leading dimension of T
      integer nc        ! No of columns (and rows considered)

      integer pmo       ! order of Pade approximant

      complex T(lda,nc)
      complex fT(lda,nc) ! f(T). Assuming same shape as T

      complex eye(mfnkryl1,mfnkryl1)
      complex X(mfnkryl1,mfnkryl1)
      complex WK1(mfnkryl1,mfnkryl1)
      complex WK2(mfnkryl1,mfnkryl1)

      integer nt

!     GL weights and nodes
      real wj(mfnkryl1)           ! just taking a high enough number
      real xj(mfnkryl1)           ! for memory

!     Variables for zgemm 
      complex a1,b1
      character transA*1,transB*1

      integer info
      integer i,j

!     Calculate Gauss-Legendre Quadrature
      if (pmo.gt.mfnkryl1) then
        write(6,*) 'Increase Quadrature memory size for Pade Approx.'
        call exitt
      endif

      call rzero(wj,mfnkryl1)
      call rzero(xj,mfnkryl1)
      call ZWGL(xj,wj,pmo)

!     Need the weights and points in [0,1]      
      call cadd(xj,1.,pmo)
      call cmult(xj,0.5,pmo)
      call cmult(wj,0.5,pmo)

!     Complex Identity matrix 
      call nek_zzero(eye,mfnkryl1*mfnkryl1)
      do i=1,nc
        eye(i,i)=complex(1.0,0)
      enddo 

      nt = mfnkryl1*mfnkryl1
      call nek_zzero(X,nt)
      do i=1,nc
      do j=1,nc
        X(i,j)=T(i,j)-eye(i,j)
      enddo        
      enddo 

!     Evaluate Pade Approximant 
      call nek_zzero(fT,lda*nc)
      call nek_zzero(WK1,nt)
      call nek_zzero(WK2,nt)
      do i=1,pmo
        call nek_zcopy(WK2,X,nt)
        call nek_zrcmult(WK2,xj(i),nt)
        call nek_zadd2(WK2,eye,nt)

!       Matrix inversion of a Triangular Matrix 
        call ztrtri('U','N',nc,WK2,mfnkryl1,info)
        if (info.ne.0) then
          write(6,*) 'Matrix inversion unsuccesfull in MAT_LN_PADE'
          call exitt
        endif

        call nek_zcopy(WK1,X,nt)
        call nek_zrcmult(WK1,wj(i),nt)

!       fT= fT + WK2*WK1
        a1 = complex(1.,0.)
        b1 = complex(1.,0.)
        transA = 'N'
        transB = 'N'
        call zgemm(transA,transB,nc,nc,nc,a1,WK2,mfnkryl1,
     $             WK1,mfnkryl1,b1,fT,lda) 

      enddo       

      return            
      end subroutine TRIMAT_LN_PADE            

!======================================================================
!======================================================================

      subroutine MAT_ZFCN_SQRT(fA,A,wkA,lda,nc,ifinv)

!     A    = U*T*U^{-1}      ! Schur decomposition             
!     f(A) = U*f(T)*U^{-1}

!     Calculating Matrix Square root
!     with Schur decomposition and
!     Algorithm 6.3 pg 136 in 
!     Higham (2008) Functions of Matrices - Theory and Computation

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'MFN.inc'

      integer lda       ! leading dimension of A
      integer ldu       ! leading dimension of U
      integer nc        ! No of columns (and rows considered)

      integer pmo       ! order of Pade approximant
      logical ifinv     ! If get inverse (fA)^{-1}

      complex A(lda,nc),U(lda,nc)
      complex fA(lda,nc) ! f(A). Assuming same shape as A.
      complex w(nc)
      complex wkA(lda,nc)     ! Work array
       

!     Complex Schur Decomposition (Double-Precision)
      ldu=lda 
      call wrp_zschur(A,lda,nc,U,ldu,w)
!     A now contains the Complex Upper-Triangular Matrix
!     U has the Schur vectors

!     Debugging
!      call write_zmat(A,lda,nc,nc,'SqT')

!     Square root of upper-triangular matrix      
      call trimat_sqrt(fA,A,lda,nc) 

!     Multiply by Schur vectors
      call SCHURVEC_MULT(fA,U,A,wkA,lda,nc,ifinv)



      return
      end subroutine MAT_ZFCN_SQRT
!----------------------------------------------------------------------
      subroutine TRIMAT_SQRT(fT,T,lda,nc)

!     f(T) = sqrt(T)
           
!     Calculating Matrix Square root for a triangular matrix
!     Algorithm 6.3 pg 136 in 
!     Higham (2008) Functions of Matrices - Theory and Computation

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'MFN.inc'

      integer lda       ! leading dimension of A
      integer nc        ! No of columns (and rows considered)

      complex T(lda,nc)
      complex fT(lda,nc) ! f(A). Assuming same shape as A.

      integer i,j,k

      complex uijk

!     Algorithm 6.3. Higham (2008) Functions of Matrices -
!     Theory and Computation. pg 136      
      call nek_zzero(fT,lda*nc)
      do i=1,nc
        fT(i,i)=sqrt(T(i,i))
      enddo        
      
      do j=2,nc
        do i=j-1,1,-1
          uijk = complex(0.,0.)
          do k=i+1,j-1
            uijk=uijk+fT(i,k)*fT(k,j)
          enddo
          fT(i,j)=(T(i,j) - uijk)/(fT(i,i) + fT(j,j))
        enddo
      enddo        


      return
      end subroutine TRIMAT_SQRT            
!======================================================================       
!---------------------------------------------------------------------- 

      subroutine MATF_ALPHAP(A,wk1,wk2,lda,nc,p)

      implicit none

      integer lda,nc
      integer p

      complex A(lda,nc)
      complex wk1(lda,nc)
      complex wk2(lda,nc)


      return
      end subroutine MATF_ALPHAP
!----------------------------------------------------------------------
      subroutine write_zmat(A,lda,r,c,nam)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'

      integer lda,r,c
      complex A(lda,c)
      character outfmt*38
      character nam*3

      integer i,j

      call blank(outfmt,38)
      write(outfmt,'(A7,I3.3,A27)') '(A3,1x,',c,
     $                '(E12.4E2,1x,E12.4E2,A1,1x))'

      if (nid.eq.0) then
        do i=1,r
          write(6,outfmt) nam, (real(A(i,j)),imag(A(i,j)),'i',j=1,c)
        enddo
      endif

      return
      end subroutine write_zmat

c-----------------------------------------------------------------------

      subroutine fwrite_zmat(A,lda,r,c,nam)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'

      integer lda,r,c
      complex A(lda,c)
      character outfmt*38
      character nam*3
      character fname*6

      logical ifexist

      integer i,j

      call blank(outfmt,38)
      write(outfmt,'(A1,I3.3,A27)') '(',c,
     $                '(E25.16E2,1x,E25.16E2,1x))'

      i=0
      ifexist=.true.
      do while (ifexist)
        i=i+1
        write(fname,'(A3,I3.3)') nam,i
        inquire(file=fname,exist=ifexist)
      enddo 
      
      if (nid.eq.0) then
        open(unit=121,file=fname,form='formatted',status='new')
        do i=1,r
          write(121,outfmt) (real(A(i,j)),imag(A(i,j)),j=1,c)
        enddo
        close(121)
      endif


      return
      end subroutine fwrite_zmat

c-----------------------------------------------------------------------
     








