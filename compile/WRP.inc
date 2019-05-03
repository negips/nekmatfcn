!	Author      : Prabal Negi
!	Description : LAPACK INTERFACE 
!                 : Memory controlling parameters are set here
!
!---------------------------------------------------------------------- 

!     Krylov sub-space size
      integer lkryl     ! No. of Krylov vectors
      integer lkryl1    ! lKRYL+1

      parameter (lkryl=20)
      parameter (lkryl1=lkryl+1)

!     Length of Krylov vector     
      integer vlen,prlen,tlen,fsilen
      parameter (vlen   = lx1*ly1*lz1*lelv*ldim)
      parameter (prlen  = lx2*ly2*lz2*lelv)
      parameter (tlen   = lx1*ly1*lz1*lelv*0)   ! Modify for passive scalars
      parameter (fsilen = 0)

      integer qlen
      parameter (qlen = vlen + prlen + tlen + fsilen)

!---------------------------------------------------------------------- 

