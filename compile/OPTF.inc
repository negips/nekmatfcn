!	Author      : Prabal
!	Description : Optimal forcing
!
!---------------------------------------------------------------------- 

	integer opfkryl
	parameter (opfkryl=10)

      integer vlen3,prlen3,tlen3,fsilen3
      parameter (vlen3  = lx1*ly1*lz1*lelv*ldim)
      parameter (prlen3 = lx2*ly2*lz2*lelv)
      parameter (tlen3  = lx1*ly1*lz1*lelv*0)
      parameter (fsilen3 = 0)

	integer qlen3
	parameter (qlen3 = vlen3 + prlen3 + tlen3 + fsilen3)

