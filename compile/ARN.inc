!	Author:	Prabal
!	Description: Testing Arnoldi method
!
!---------------------------------------------------------------------- 

	integer arnkryl
	parameter (arnkryl=10)

      integer vlen2,prlen2,tlen2,fsilen2
      parameter (vlen2  = lx1*ly1*lz1*lelv*ldim)
      parameter (prlen2 = lx2*ly2*lz2*lelv)
      parameter (tlen2  = lx1*ly1*lz1*lelv*0)
      parameter (fsilen2 = 0)

	integer qlen2
	parameter (qlen2 = vlen2 + prlen2 + tlen2 + fsilen2)


