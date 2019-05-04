c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      udiff  = 0
      utrans = 0

      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)


      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'NEKUSE_DEF'
      include 'NEKUSE'
      include 'RTFILTER'

      integer ix,iy,iz,ieg,iel
      
      iel=gllel(ieg)

      ffx = 0.
      ffy = 0.
      ffz = 0.

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      qvol   = 0.0
      source = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'NEKUSE'
      include 'PARALLEL'
      include 'WING_MVMSH'
      include 'FSI'

      integer iel

      iel=gllel(ieg)
      ux = 0.
      uy = 0.
      uz = 0.

      if (cbu.eq.'v  ') then
        if (.not.ifpert) then
          ux = 1.
          uy = 0.
          uz = 0.
        endif  
      elseif (cbu.eq.'mv ') then
        ux = umeshx(ix,iy,iz,iel)
        uy = umeshy(ix,iy,iz,iel)
        if (ndim.eq.3) uz = umeshz(ix,iy,iz,iel)

      endif  

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'NEKUSE_DEF'
      include 'NEKUSE'
      include 'INPUT_DEF'
      include 'INPUT'

      integer ix,iy,iz,ieg
      real amp, ran

c     velocity
c     random distribution

      amp = param(110)

      ran = 3.e4*(ieg+x*sin(y)) - 1.5e3*ix*iy + .5e5*ix
      ran = 1.e3*sin(ran)
      ran = 1.e3*sin(ran)
      ran = cos(ran)
      ux  = ran*amp

      ran = 2.3e4*(ieg+x*sin(y)) + 2.3e3*ix*iy - 2.e5*ix
      ran = 1.e3*sin(ran)
      ran = 1.e3*sin(ran)
      ran = cos(ran)
      uy  = ran*amp

      uz = 0.0


      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat

      implicit none  

      include 'SIZE_DEF'
      include 'SIZE'
      include 'ADJOINT_DEF'
      include 'ADJOINT'

!      if (param(33) .eq. 0.0) param(33) = 1.0
!      if (param(34) .eq. 0.0) param(34) = 1.0e-2

      call uprm_read

!     Adjoints
!      ifadj = .true.

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'

!      call my_mvmsh_init

!      param(59) = 1  ! all elements deformed

      ifxyo = .true.
!      ifusermv = .true.  ! define our own mesh velocity
!      ifstrs = .true.

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk

      implicit none

      INCLUDE 'SIZE_DEF'      
      INCLUDE 'SIZE'
      INCLUDE 'SOLN_DEF'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'
      INCLUDE 'INPUT_DEF'
      INCLUDE 'INPUT'

      integer lt,lt2
      parameter (lt=lx1*ly1*lz1*lelv)
      parameter (lt2=lx2*ly2*lz2*lelv)

      real vort,w1,w2
      common /scrns/ vort(lt,3), w1(lt), w2(lt)

      integer n,m
      parameter (n=10)
      real matA(n,n),matV(n,n)

      complex cA(n,n),cV(n,n),cB(n,n)
      integer pord                          ! Pade Order

      integer seed

      integer i,j
      logical ifinv
      real rr,ri


      m=5

      if (istep.eq.0) then

        call rzero(matA,m)
        seed = 86458
        call srand(seed)

        call nek_zzero(cA,n*n)
        call nek_zzero(cB,n*n)

        do i=1,m
        do j=1,m
          rr = rand()
          ri = rand()
          cA(i,j)   = complex(rr,ri)
          cB(i,j)   = complex(rr,ri)
        enddo
        enddo

        write(6,*) ' '
        call write_zmat(cA,n,m,m,'Ain')

        ifinv = .false.
        call MAT_ZFCN(cV,cA,n,m,'loge',ifinv)

        call write_zmat(cV,n,m,m,'fAo')

        pord = 16
        ifinv = .false.
        call MAT_ZFCN_LN(cV,cB,n,m,pord,ifinv) 
        call write_zmat(cV,n,m,m,'Pde')
      

      endif

      call exitt


      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c
c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)
      return
      end
c
c automatically added by makenek
      subroutine cmt_switch ! to set IFCMT logical flag
      include 'SIZE'
      include 'INPUT'
      IFCMT=.false.
      return
      end
c
c automatically added by makenek
      subroutine usrflt(rmult) ! user defined filter
      include 'SIZE'
      real rmult(lx1)
      call rone(rmult,lx1)
      return
      end
c
c automatically added by makenek
      subroutine userflux ! user defined flux
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      real fluxout(lx1*lz1)
      return
      end
c
c automatically added by makenek
      subroutine userEOS ! user defined EOS 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      return
      end
