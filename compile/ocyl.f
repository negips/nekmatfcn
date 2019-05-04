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

!     Just need the pointer from SOLN      
      integer jp
      common /ppointr/ jp

      integer ix,iy,iz,ieg,iel

!     Reshaped array for convenience      
      real optfx(lx1,ly1,lz1,lelt,2)
      real optfy(lx1,ly1,lz1,lelt,2)
      real optfz(lx1,ly1,lz1,lelt,2)
      common /optfxyz/ optfx,optfy,optfz
     
      iel=gllel(ieg)

      ffx = 0.
      ffy = 0.
      ffz = 0.

      if ((jp.gt.0).and.(jp.le.2)) then
        ffx = optfx(ix,iy,iz,iel,jp)
        ffy = optfy(ix,iy,iz,iel,jp)
        if (ndim.eq.3) ffz = optfz(ix,iy,iz,iel,jp)
      endif        
            
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
        ux = 0. ! umeshx(ix,iy,iz,iel)
        uy = 0. ! umeshy(ix,iy,iz,iel)
        if (ndim.eq.3) uz = 0. ! umeshz(ix,iy,iz,iel)

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
      include 'USERPAR'

      integer ix,iy,iz,ieg
      real amp, ran

c     velocity
c     random distribution

      amp = UPRM_DAMPL

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
      INCLUDE 'MVGEOM_DEF'
      INCLUDE 'MVGEOM'
      INCLUDE 'MATFCN'

      integer lt,lt2
      parameter (lt=lx1*ly1*lz1*lelv)
      parameter (lt2=lx2*ly2*lz2*lelv)

      real vort,w1,w2
      common /scrns/ vort(lt,3), w1(lt), w2(lt)

      real optfx(lt,2),optfy(lt,2),optfz(lt,2)
      common /optfxyz/ optfx,optfy,optfz

      real Omega

!!     Hard-coding parameters for now
!      ifmatf = .true.
!      matf_ifpr = .true.
!      matf_uzawa = .false.
!       
!      ngs=1        ! no of Gram-Schmidt Orthogonalizations
!      northo=90    ! no of Krylov vectors to save
!      sstep=100    ! no of iterations between Re-Ortho
!      Omega = 0.5  ! Angular frequency of Forcing

      call opzero(wx,wy,wz)

      if (ifmatf) then
        call MATF_MAIN

      else
        if (mod(istep,sstep).eq.0) then
          time=0.
          istep=0
        endif
      endif

      call opcopy(optfx(1,1),optfy(1,1),optfz(1,1),
     $                vxp(1,2),vyp(1,2),vzp(1,2))
      call opcmult(optfx(1,1),optfy(1,1),optfz(1,1),matf_omega)

      call opcopy(optfx(1,2),optfy(1,2),optfz(1,2),
     $                vxp(1,1),vyp(1,1),vzp(1,1))
      call opcmult(optfx(1,2),optfy(1,2),optfz(1,2),-matf_omega)


      if (istep.eq.0) then
        call outpost(vxp(1,1),vyp(1,1),vzp(1,1),prp(1,1),
     $      tp(1,1,1),'pr1') 
        call outpost(vxp(1,2),vyp(1,2),vzp(1,2),prp(1,2),
     $      tp(1,1,2),'pr2')

        call outpost(optfx(1,1),optfy(1,1),optfz(1,1),prp(1,1),
     $      tp(1,1,1),'pf1') 

      endif        

      return
      end
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
