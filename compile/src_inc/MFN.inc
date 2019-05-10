!	Author      : Prabal Negi
!	Description : Matrix Functions 
!                 : Memory controlling parameters are set here
!
!---------------------------------------------------------------------- 

!     Krylov sub-space for Matrix Function evaluation
      integer mfnkryl     ! No. of Krylov vectors
      integer mfnkryl1    ! MFNKRYL+1

      parameter (mfnkryl=100)
      parameter (mfnkryl1=mfnkryl+1)

!     Krylov sub-space for Outer Krylov space
      integer gmrkryl,gmrkryl1
      parameter (gmrkryl=20)
      parameter (gmrkryl1=gmrkryl+1)


!---------------------------------------------------------------------- 

