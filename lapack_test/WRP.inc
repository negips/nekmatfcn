!	Author      : Prabal Negi
!	Description : LAPACK INTERFACE 
!                 : Memory controlling parameters are set here
!
!---------------------------------------------------------------------- 

!     Krylov sub-space size
      integer lkryl     ! No. of Krylov vectors
      integer lkryl1    ! lKRYL+1

      parameter (lkryl=10)
      parameter (lkryl1=lkryl+1)

!     Length of Krylov vector     
      integer vlen,prlen,tlen,fsilen
      parameter (vlen   = 10)
      parameter (prlen  = 1)
      parameter (tlen   = 1)   ! Modify for passive scalars
      parameter (fsilen = 1)

      integer qlen
      parameter (qlen = vlen + prlen + tlen + fsilen)

!---------------------------------------------------------------------- 

