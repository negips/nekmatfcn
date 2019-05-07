!	Author      : Prabal Negi
!	Description : Matrix Functions 
!                 : Memory controlling parameters are set here
!
!---------------------------------------------------------------------- 

!     Krylov sub-space for Matrix Function evaluation
      integer mfnkryl     ! No. of Krylov vectors
      integer mfnkryl1    ! MFNKRYL+1

      parameter (mfnkryl=20)
      parameter (mfnkryl1=mfnkryl+1)


!---------------------------------------------------------------------- 

