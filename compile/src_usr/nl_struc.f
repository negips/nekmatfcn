!======================================================================
!     Routines for fluid structure interaction implicit solutions
!     Author: Prabal S. Negi
!     Based on fixed-point method with dynamic relaxation
!
!======================================================================       
!====================================================================== 
      subroutine nlfsi_main

      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'ADJOINT_DEF'
      include 'ADJOINT'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'MVGEOM_DEF'
      include 'MVGEOM'
      include 'NLFSI'
      include 'CHKPOINT'

      real scale
      real x0(3)
      logical ifdout,iftout
      integer ierr

      integer icalld
      save icalld
      data icalld /0/

!     Takes from Calc_torque_axis/pert 
      real dragx, dragpx, dragvx,
     $     dragy, dragpy, dragvy,
     $     dragz, dragpz, dragvz,

     $     torqx, torqpx, torqvx,
     $     torqy, torqpy, torqvy,
     $     torqz, torqpz, torqvz,

     $     dpdx_mean,dpdy_mean,dpdz_mean,
     $     dgtq

      common /ctorq/ dragx(0:maxobj),dragpx(0:maxobj),dragvx(0:maxobj)
     $             , dragy(0:maxobj),dragpy(0:maxobj),dragvy(0:maxobj)
     $             , dragz(0:maxobj),dragpz(0:maxobj),dragvz(0:maxobj)
c
     $             , torqx(0:maxobj),torqpx(0:maxobj),torqvx(0:maxobj)
     $             , torqy(0:maxobj),torqpy(0:maxobj),torqvy(0:maxobj)
     $             , torqz(0:maxobj),torqpz(0:maxobj),torqvz(0:maxobj)
c
     $             , dpdx_mean,dpdy_mean,dpdz_mean
     $             , dgtq(3,4)

      integer lvp,lvp2
      parameter(lvp=lpx1*lpy1*lpz1*lpelv)
      parameter(lvp2=lpx2*lpy2*lpz2*lpelv)
      real wk_vxp(lvp),wk_vyp(lvp),wk_vzp(lvp),wk_prp(lvp2)

      real wk_vx0(lvp),wk_vy0(lvp),wk_vz0(lvp),wk_pr0(lvp2)
      real tmp_pr(lvp)

      real rot_sn(3,3)
      real rot_cs(3,3)
      real dxyz(3),r0(3)
      real psi_sx,psi_sy,psi_sz

      real*8 dnekclock          ! function. timing
      real*8 steptime           ! timing for step

      real rxyz 
      common /scrsf/ rxyz(lx1*ly1*lz1*lelt,3)

      integer i

      if (.not.IFNLFSI) then
        return
      endif  

      if (ifpert.and.ifadj) then
        call nlfsi_main_adj
        return
      endif

      if (icalld.eq.0.and..not.nlfsi_ifinit) then
!       Initialize timer 
        nlfsi_timea = 0.
        nlfsi_timei = 0.
           
        call nlfsi_init
        
        icalld=icalld+1
      endif  


!     Start timing for this step   
      steptime = dnekclock()

      if (ifchkptrst.and.istep.lt.chkptnrsf) then
!       Read saved values for restart
!       If the simulation has been restarted
        call nlfsi_rstread

        call nlfsi_output

        if (ifusermv) then
          if (ifpert) then
            call nlfsi_meshv_pert
          else  
            call nlfsi_meshv
          endif
        endif  

        nlfsi_timei=dnekclock()-steptime
        nlfsi_timea=nlfsi_timea+nlfsi_timei
        return

      elseif (istep.eq.0) then
!       If there are no restarts              
!       At first time step the projected velocity is
!       just the initialized velocity
        psiv_s = psiv_ini
        psiv_m = psiv_ini
        if (ifadj) psiv_m = 0. 
        psiv   = psiv_ini
        psia   = 0.
        psia_s = 0.    
        psi    = psi_ini
        psi_s  = psi + psiv_m*DT

        call nlfsi_output

        if (ifusermv) then
          if (ifpert) then
            call nlfsi_meshv_pert
          else  
            call nlfsi_meshv
          endif
        endif  

        nlfsi_timei=dnekclock()-steptime
        nlfsi_timea=nlfsi_timea+nlfsi_timei
        return
      endif

      x0(1) = nlfsi_x0
      x0(2) = nlfsi_y0
      if (ndim.eq.3) x0(3) = 1000.         ! set outside the domain

      ifdout=.false.
      iftout=.false.

      scale=nlfsi_rescale
      
      if (ifpert) then
        call torque_calc_axis(scale,x0,ifdout,iftout,
     $      vxp(1,jp),vyp(1,jp),vzp(1,jp),prp(1,jp))
      else            
        call torque_calc_axis(scale,x0,ifdout,iftout,vx,vy,vz,pr)
      endif  

!     Who came up with this code? 
      call copy(nlfs_dragx (0),dragx (0),maxobj+1)
      call copy(nlfs_dragpx(0),dragpx(0),maxobj+1)
      call copy(nlfs_dragvx(0),dragpx(0),maxobj+1)
      call copy(nlfs_dragy (0),dragy (0),maxobj+1)
      call copy(nlfs_dragpy(0),dragpy(0),maxobj+1)
      call copy(nlfs_dragvy(0),dragvy(0),maxobj+1)
      call copy(nlfs_dragz (0),dragz (0),maxobj+1)
      call copy(nlfs_dragpz(0),dragpz(0),maxobj+1)
      call copy(nlfs_dragvz(0),dragvz(0),maxobj+1)
      call copy(nlfs_torqx (0),torqx (0),maxobj+1)
      call copy(nlfs_torqpx(0),torqpx(0),maxobj+1)
      call copy(nlfs_torqvx(0),torqpx(0),maxobj+1)
      call copy(nlfs_torqy (0),torqy (0),maxobj+1)
      call copy(nlfs_torqpy(0),torqpy(0),maxobj+1)
      call copy(nlfs_torqvy(0),torqvy(0),maxobj+1)
      call copy(nlfs_torqz (0),torqz (0),maxobj+1)
      call copy(nlfs_torqpz(0),torqpz(0),maxobj+1)
      call copy(nlfs_torqvz(0),torqvz(0),maxobj+1)

      call nlfsi_struc_solve_acc

      nlfsi_timei=dnekclock()-steptime
      nlfsi_timea=nlfsi_timea+nlfsi_timei

      return
      end subroutine nlfsi_main 
!---------------------------------------------------------------------- 
      subroutine nlfsi_struc_solve_acc

!     Calculate superposition of Solutions to get implicit update
!     In this routine we start with an assumed or extrapolated value of
!     psi (psi_s). Then solve the structural equations to obtain the new
!     value to psi and psiv.
!     Convergence is achieved if the new and extrapolated values are
!     consistent within a given tolerance.            

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'GEOM_DEF'      
      include 'GEOM'
      include 'NLFSI'

      integer icalld
      save icalld
      data icalld /0/

      real ext_k,bd_psi,bd_psiv,bd_psia
      real ab0,ab1,ab2
      integer iobj
      integer ilag
      real const
      real solid_area
      real fb

      real lhs,rhs
!      real alpha
      real Fs,Fd,Fk,Ft

      integer irst,NBDMSH,NABMSH
      integer istep_old

      integer it
      real rnorm,dr(2),dr_1(2)
      real omg,tmp
      real omegamin
      parameter (omegamin=1.0e-10)
      logical ifbisection

      real inert,stiff,damp

      inert=nlfsi_inertia
      stiff=nlfsi_stiff
      damp=nlfsi_damp

      it = nlfsi_iter

      irst=0
!     ABBD is not set if there is no moving boundary condition
      IF (.not.IFMVBD) THEN
         NBDMSH = 1
         NABMSH = PARAM(28)
         IF (NABMSH.GT.ISTEP .and. irst.le.0) NABMSH = ISTEP
         IF (IFSURT)          NABMSH = NBD
         CALL RZERO   (ABMSH,10)
         CALL SETABBD (ABMSH,DTLAG,NABMSH,NBDMSH)
      ENDIF

!     Only doing for one object right now 
      iobj=0
      if (nlfsi_ifrot) then
!       Rotational coordinate system is opposite to the one used in Nek 
        Fs  = -nlfs_torqz(iobj) 
      else  
        Fs  = nlfs_dragy(iobj) 
      endif

!     Backward differentiation for psi
      bd_psi = psi*bd(2)
      do ilag=2,NBD
        bd_psi = bd_psi + bd(ilag+1)*psilag(ilag-1)
      enddo

!     Backward differentiation for psiv
      bd_psiv = psiv*bd(2)
      do ilag=2,NBD
        bd_psiv = bd_psiv + bd(ilag+1)*psivlag(ilag-1)
      enddo

!     Backward differentiation for psia
      bd_psia = psia*bd(2)
      do ilag=2,NBD
        bd_psia = bd_psia + bd(ilag+1)*psialag(ilag-1)
      enddo
     
      lhs = bd(1)/DT*inert + stiff*DT/bd(1) - damp
      rhs = -stiff*bd_psi/bd(1) + inert*bd_psiv/DT +
     $       Fs

!     Velocity
      psiv_iter(it) = rhs/lhs
!     Position
      psi_iter(it) = (psiv_iter(it)*DT + bd_psi)/bd(1)
!     Acceleration
      psia_iter(it) = (bd(1)*psiv_iter(it)-bd_psiv)/DT

!     Calculate residual
      psi_res(it) = psi_iter(it) - psi_s 
      psiv_res(it) = psiv_iter(it) - psiv_s
      psia_res(it) = psia_iter(it) - psia_s

!     Check convergence
      rnorm = abs(psia_res(it))

      nlfsi_rnorm(it) = rnorm      

      if (rnorm.le.nlfsi_tol) then
        nlfsi_ifconv = .true.
      endif  

      if (it.eq.nlfsi_miters.and..not.nlfsi_ifconv) then
        if (nio.eq.0) write(6,*) 'FSI unconverged', rnorm
        nlfsi_ifconv = .true.
      endif  

!     update omega if not converged
      ifbisection=.false.            ! use simple bisection
      if (nlfsi_ifconv) then
        
!       Copy lag array for psi.
        do ilag=lorder-1,2,-1
          psilag(ilag)=psilag(ilag-1)
        enddo
        psilag(1)=psi
        psi = psi_s

!       Copy lag array for psiv 
        do ilag=lorder-1,2,-1
          psivlag(ilag)=psivlag(ilag-1)
        enddo
        psivlag(1)=psiv
        psiv = psiv_s

!       Copy lag array for psia
        do ilag=lorder-1,2,-1
          psialag(ilag)=psialag(ilag-1)
        enddo
        psialag(1)=psia
        psia = psia_s

!       Copy lag array for psivm
        do ilag=lorder-1,2,-1
          psivmlag(ilag)=psivmlag(ilag-1)
        enddo
        psivmlag(1)=psiv_m

!       Extrapolate accleration for next time step 
        psia_s = ab(1)*psia + ab(2)*psialag(1) + ab(3)*psialag(2)

        nlfsi_omegai(it) = 0.
        omega_prev = nlfsi_omegai(it-1)

!       Recalculate Backward difference terms.
!       Needed for estimating vel and pos for next step            
!       Backward differentiation for psi
        bd_psi = psi*bd(2)
        do ilag=2,NBD
          bd_psi = bd_psi + bd(ilag+1)*psilag(ilag-1)
        enddo

!       Backward differentiation for psiv
        bd_psiv = psiv*bd(2)
        do ilag=2,NBD
          bd_psiv = bd_psiv + bd(ilag+1)*psivlag(ilag-1)
        enddo

      else        ! if not converged
        if (it.le.1) then
          omg = omega_prev          ! From last time step
          nlfsi_omegai(it) = omg
        elseif (ifbisection) then
!         using simple bisection for now  
          omg = 0.5
          nlfsi_omegai(it) = omg
        else

!           This one works somehow              
!!          tmp = psia_iter(it) - psia_iter(it-1)
!!          omg = tmp/(psia_res(it-1) - psia_res(it))
!!
!!!         Can't be negative
!!          omg = max(omg,1.e-6) 
!!          nlfsi_omegai(it) = omg


          tmp = -nlfsi_omegai(it-1)*psia_res(it)
          omg = tmp/(psia_res(it) - psia_res(it-1))
          nlfsi_omegai(it) = omg

        endif

!       Acceleration for next iteration            
        psia_s = psia_iter(it) + omg*(-psia_res(it))

      endif 

!     Calculate velocity using backward difference
      psiv_s = (psia_s*DT + bd_psiv)/bd(1)

!     Calculate position using backward difference
      psi_s = (psiv_s*DT + bd_psi)/bd(1)

!     Estimate mesh velocity for this position
!     This is different from surface velocity because mesh movement is
!     done using explicit integration. So we apply velocities on the
!     mesh such that the extrapolated surface position are correct
      psiv_m =  1/ABMSH(1)*((psi_s-psi)/DT-(ABMSH(2)*psivmlag(1)
     $     + ABMSH(3)*psivmlag(2)))

      Fk = -stiff*psi
      Fd = damp*psiv
      Ft = Fs

!     Save forces for output 
      nlfsi_Ff  = Fs
      nlfsi_Fk  = Fk
      nlfsi_Fd  = Fd
      nlfsi_Ft  = Ft

!     prabal
!      if (nio.eq.0) write(6,12) 'iters',it,psia_iter(it),
!     $         psia_res(it),Fs,omg
!   12 format(A5,1x,I3,1x,4(1pE14.7,1x))

      return
      end subroutine nlfsi_struc_solve_acc

!---------------------------------------------------------------------- 

      subroutine nlfsi_struc_solve_vel

!     Calculate superposition of Solutions to get implicit update
!     In this routine we start with an assumed or extrapolated value of
!     psi (psi_s). Then solve the structural equations to obtain the new
!     value to psi and psiv.
!     Convergence is achieved if the new and extrapolated values are
!     consistent within a given tolerance.
!     In this routine velocity is extrapolated 

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'GEOM_DEF'      
      include 'GEOM'
      include 'NLFSI'

      integer icalld
      save icalld
      data icalld /0/

      real ext_k,bd_psi,bd_psiv,bd_psia
      real ab0,ab1,ab2
      integer iobj
      integer ilag
      real const
      real solid_area
      real fb

      real lhs,rhs
      real Fs,Fk,Fd,Ft

      integer irst,NBDMSH,NABMSH
      integer istep_old

      integer it
      real rnorm,dr(2),dr_1(2)
      real omg,tmp
      logical ifbisection

      real inert,stiff,damp

      inert=nlfsi_inertia
      stiff=nlfsi_stiff
      damp=nlfsi_damp

      it = nlfsi_iter

      irst=0
!     ABBD is not set if there is no moving boundary condition
      IF (.not.IFMVBD) THEN
         NBDMSH = 1
         NABMSH = PARAM(28)
         IF (NABMSH.GT.ISTEP .and. irst.le.0) NABMSH = ISTEP
         IF (IFSURT)          NABMSH = NBD
         CALL RZERO   (ABMSH,10)
         CALL SETABBD (ABMSH,DTLAG,NABMSH,NBDMSH)
      ENDIF

!     Only doing for one object right now 
      iobj=0
      if (nlfsi_ifrot) then
!       Rotational coordinate system is opposite to the one used in Nek 
        Fs  = -nlfs_torqz(iobj) 
      else  
        Fs  = nlfs_dragy(iobj) 
      endif

!     Backward differentiation for psi
      bd_psi = psi*bd(2)
      do ilag=2,NBD
        bd_psi = bd_psi + bd(ilag+1)*psilag(ilag-1)
      enddo

!     Backward differentiation for psiv
      bd_psiv = psiv*bd(2)
      do ilag=2,NBD
        bd_psiv = bd_psiv + bd(ilag+1)*psivlag(ilag-1)
      enddo

!     Backward differentiation for psia
      bd_psia = psia*bd(2)
      do ilag=2,NBD
        bd_psia = bd_psia + bd(ilag+1)*psialag(ilag-1)
      enddo

      lhs = bd(1)/DT*inert + stiff*DT/bd(1) - damp
      rhs = -stiff*bd_psi/bd(1) + inert*bd_psiv/DT +
     $       Fs

!     Velocity
      psiv_iter(it) = rhs/lhs
!     Position
      psi_iter(it) = (psiv_iter(it)*DT + bd_psi)/bd(1)
!     Acceleration
      psia_iter(it) = (bd(1)*psiv_iter(it)-bd_psiv)/DT

!     Calculate residual
      psi_res(it) = psi_iter(it) - psi_s 
      psiv_res(it) = psiv_iter(it) - psiv_s
      psia_res(it) = psia_iter(it) - psia_s

!     Check convergence
      rnorm = abs(psia_res(it))

      nlfsi_rnorm(it) = rnorm      

      if (rnorm.le.nlfsi_tol) then
        nlfsi_ifconv = .true.
      endif  

!     prabal
!      nlfsi_ifconv=.true.
!      psi_s = 0.
!      psiv_s = 0.
!      psia_s = 0.
!      psiv_m = 0.

      if (it.eq.nlfsi_miters.and..not.nlfsi_ifconv) then
        if (nio.eq.0) write(6,*) 'FSI unconverged', rnorm
        nlfsi_ifconv = .true.
      endif  

!     update omega if not converged
      ifbisection=.false.
      if (nlfsi_ifconv) then
        
!       Copy lag array for psi.
        do ilag=lorder-1,2,-1
          psilag(ilag)=psilag(ilag-1)
        enddo
        psilag(1)=psi
        psi = psi_s

!       Copy lag array for psiv 
        do ilag=lorder-1,2,-1
          psivlag(ilag)=psivlag(ilag-1)
        enddo
        psivlag(1)=psiv
        psiv = psiv_s

!       Copy lag array for psia
        do ilag=lorder-1,2,-1
          psialag(ilag)=psialag(ilag-1)
        enddo
        psialag(1)=psia
        psia = psia_s

!       Copy lag array for psivm
        do ilag=lorder-1,2,-1
          psivmlag(ilag)=psivmlag(ilag-1)
        enddo
        psivmlag(1)=psiv_m

!       Extrapolate velocity for next time step 
        psiv_s = ab(1)*psiv + ab(2)*psivlag(1) + ab(3)*psivlag(2)

        nlfsi_omegai(it) = 0.
        omega_prev = nlfsi_omegai(it-1)

!       Recalculate Backward difference terms.
!       Needed for estimating vel and pos for next step            
!       Backward differentiation for psi
        bd_psi = psi*bd(2)
        do ilag=2,NBD
          bd_psi = bd_psi + bd(ilag+1)*psilag(ilag-1)
        enddo

!       Backward differentiation for psiv
        bd_psiv = psiv*bd(2)
        do ilag=2,NBD
          bd_psiv = bd_psiv + bd(ilag+1)*psivlag(ilag-1)
        enddo

      else        ! if not converged
        if (it.le.1) then
          omg = omega_prev
          nlfsi_omegai(it) = omg
        elseif (ifbisection) then
!         using simple bisection for now  
          omg = 0.5
          nlfsi_omegai(it) = omg
        else
!          dr(1) = psi_res(it) - psi_res(it-1)
!          dr(2) = psiv_res(it) - psiv_res(it-1)
!          tmp = psi_res(it-1)*dr(1) + psiv_res(it-1)*dr(2)
!          tmp = tmp/(dr(1)**2 + dr(2)**2)
!          omg = -nlfsi_omegai(it-1)*tmp

          tmp = psiv_iter(it) - psiv_iter(it-1)
          omg = tmp/(psiv_res(it-1) - psiv_res(it))

!          omg = min(omg,1.)
          omg = max(omg,1.0-6) 
          nlfsi_omegai(it) = omg  
        endif

!       Velocity for next iteration 
        psiv_s = psiv_iter(it) + omg*(-psiv_res(it))

      endif 

!     Calculate acceleration
      psia_s = (bd(1)*psiv_s-bd_psiv)/DT

!     Calculate position using backward difference
      psi_s = (psiv_s*DT + bd_psi)/bd(1)

!     Estimate mesh velocity for this position
!     This is different from surface velocity because mesh movement is
!     done using explicit integration. So we apply velocities on the
!     mesh such that the extrapolated surface position are correct
      psiv_m =   1/ABMSH(1)*((psi_s-psi)/DT-(ABMSH(2)*psivmlag(1)
     $     + ABMSH(3)*psivmlag(2)))

      Fk = -stiff*psi
      Fd = damp*psiv
      Ft = Fs + Fk + Fd

!     Save forces for output 
      nlfsi_Ff  = Fs
      nlfsi_Fk  = Fk
      nlfsi_Fd  = Fd
      nlfsi_Ft  = Ft

!     prabal
      if (nio.eq.0) write(6,12) 'iters',it,psia_iter(it),
     $         psia_res(it),tmp,omg
   12 format(A5,1x,I3,1x,4(1pE14.7,1x))

      return
      end subroutine nlfsi_struc_solve_vel

!---------------------------------------------------------------------- 

      subroutine nlfsi_rstsave

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'NLFSI'
      include 'CHKPOINT'

      integer i,is
      character*132 fname
      character*1   fname1(132)
      equivalence (fname1,fname)      

      character*1 exten0,exten1
      save exten0,exten1
      data exten0 /"0"/
      data exten1 /"1"/

      integer exten
      save exten
      data exten /0/

      character*1 dot
      save dot
      data dot /"."/
      character*16 outfmt

      integer iunit,ierr

      if (istep.lt.chkptnrsf) return

      is = mod(istep,chkptstep)
      if (is.ge.chkptstep) return

      ierr=0
      if (is.lt.chkptnrsf) then
        i=is+1
        nlfsi_rst_psi(i)  = psi
        nlfsi_rst_psiv(i) = psiv
      endif

      if (nid.eq.0) then
        if (is.eq.chkptnrsf-1) then
!         Create filename
          call blank(fname1,132)
          call chcopy(fname1,nlfsi_rst_fname,13)
          call chcopy(fname1(14),dot,1)
          if (exten.eq.0) then
            call chcopy(fname1(15),exten0,1)
            exten=1  
          else
            call chcopy(fname1(15),exten1,1)
            exten=0  
          endif     ! exten.eq.0

!         Write output
          call IO_freeid(iunit,ierr)     
          open(iunit,file=fname,status='unknown',action='write',
     $      iostat=ierr)
          if (ierr.eq.0) then
            write(6,'(A26,1x,A13)') 'NLFSI Restart: Saving file',
     $               fname
            write(iunit,*) (nlfsi_rst_psi(i),i=1,nlfsi_nrst)
            write(iunit,*) (nlfsi_rst_psiv(i),i=1,nlfsi_nrst)
            close(iunit)
            write(6,'(A4)')  'Done'
          endif   ! ierr.eq.0 
        endif     ! is.eq.chkptnrst-1
      endif       ! nid.eq.0

      call err_chk(ierr,'Error opening nlfsi_restart file.')
      
      return 
      end subroutine nlfsi_rstsave

!---------------------------------------------------------------------- 
      subroutine nlfsi_meshv

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'          ! xm1,ym1,zm1
      include 'MVGEOM_DEF'
      include 'MVGEOM'        ! wx,wy,wz
      include 'WING_MVMSH'
      include 'NLFSI'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'

      real ucx,ucy,ucz        ! mesh velocities
      real dx,dy,dz

      integer i,n

      n = nx1*ny1*nz1*nelv

      do i=1,n                          ! Translational velocity
!       Current time step mesh velocity at the wall is the same
!       as the velocity at the wall
        if (.not.nlfsi_ifrot) then

          ucx =  0.
          ucy =  psiv_m
          ucz =  0.0
  
          wx(i,1,1,1) = basev(i)*ucx
          wy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) wz(i,1,1,1) = basev(i)*ucz

          ucx = 0.
          ucy = psiv_s
          ucz = 0.

          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz

        else
          dx = xm1(i,1,1,1) - nlfsi_x0
          dy = ym1(i,1,1,1) - nlfsi_y0
          dz = 0.
          
          ucx =  psiv_m*dy
          ucy = -psiv_m*dx
          if (ndim.eq.3) ucz =  0.0
  
          wx(i,1,1,1) = basev(i)*ucx
          wy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) wz(i,1,1,1) = basev(i)*ucz

          ucx =  psiv_s*dy
          ucy = -psiv_s*dx
          if (ndim.eq.3) ucz =  0.0
         
          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz
         
        endif
!       This is the prescribed wall velocity for the next step
!       (userbc)
       
      enddo

      return
      end subroutine nlfsi_meshv
!---------------------------------------------------------------------- 
      subroutine nlfsi_meshv_pert

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'          ! xm1,ym1,zm1
      include 'MVGEOM_DEF'
      include 'MVGEOM'        ! wx,wy,wz
      include 'WING_MVMSH'
      include 'NLFSI'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'

      real ucx,ucy,ucz        ! mesh velocities
      real dx,dy,dz
      real psi_sx,psi_sy,psi_sz

      integer i,n
      real rot_sn(3,3)
      real rot_cs(3,3)
      real dxyz(3),r0(3)


      call rzero(rot_sn,9)            ! sine terms
      call rzero(rot_cs,9)            ! cosine terms
!     Matrix of sine terms
!     For a clockwise rotation for a positive psi
      if (nlfsi_ifrot) then
        rot_sn(1,2) =  psi_s
        rot_sn(2,1) = -psi_s
      else
        continue
      endif  

!     Matrix of cosine terms
!     In a linear case this is just Identity      
      rot_cs(1,1) = 1
      rot_cs(2,2) = 1
      rot_cs(3,3) = 1

      n = nx1*ny1*nz1*nelv

      do i=1,n                          ! Translational velocity
!       Current time step mesh velocity at the wall is the same
!       as the velocity at the wall
        if (.not.nlfsi_ifrot) then

!         Mesh motion                
          wx(i,1,1,1) = 0.
          wy(i,1,1,1) = 0.
          if (ndim.eq.3) wz(i,1,1,1) = 0.

!         Boundary condition            
          psi_sx = 0.
          psi_sy = psi_s
          psi_sz = 0.  

          ucx =  0. - (nlfsi_grvx0(i,1,1,1,1)*psi_sx + 
     $           nlfsi_grvx0(i,1,1,1,2)*psi_sy +
     $           nlfsi_grvx0(i,1,1,1,3)*psi_sz) 
          ucy =  psiv_s - (nlfsi_grvy0(i,1,1,1,1)*psi_sx + 
     $           nlfsi_grvy0(i,1,1,1,2)*psi_sy +
     $           nlfsi_grvy0(i,1,1,1,3)*psi_sz)
          if (if3d) then
            ucz =  0. - (nlfsi_grvz0(i,1,1,1,1)*psi_sx + 
     $             nlfsi_grvz0(i,1,1,1,2)*psi_sy +
     $             nlfsi_grvz0(i,1,1,1,3)*psi_sz)
          else
            ucz = 0.
          endif 

          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz

        else
!         Rotational FSI

!         Distance from Rotational axis                
          dx = xm1(i,1,1,1) - nlfsi_x0
          dy = ym1(i,1,1,1) - nlfsi_y0
          dz = 0.

!         Mesh motion 
          wx(i,1,1,1) = 0.
          wy(i,1,1,1) = 0.
          if (ndim.eq.3) wz(i,1,1,1) = 0.

!         Calculate distance from axis   
          r0(1) = dx 
          r0(2) = dy
          r0(3) = dz

!         perturbation in position
          call mxm(rot_sn,3,r0,3,dxyz,1)
          psi_sx = dxyz(1)
          psi_sy = dxyz(2)
          psi_sz = dxyz(3) 

          ucx =  psiv_s*dy
          ucy =  -psiv_s*dx
          ucz =  0.

          ucx =  ucx - (nlfsi_grvx0(i,1,1,1,1)*psi_sx + 
     $           nlfsi_grvx0(i,1,1,1,2)*psi_sy +
     $           nlfsi_grvx0(i,1,1,1,3)*psi_sz) 
          ucy =  ucy - (nlfsi_grvy0(i,1,1,1,1)*psi_sx + 
     $           nlfsi_grvy0(i,1,1,1,2)*psi_sy +
     $           nlfsi_grvy0(i,1,1,1,3)*psi_sz)
          if (if3d) then
            ucz =  ucz - (nlfsi_grvz0(i,1,1,1,1)*psi_sx + 
     $             nlfsi_grvz0(i,1,1,1,2)*psi_sy +
     $             nlfsi_grvz0(i,1,1,1,3)*psi_sz)
          else
            ucz = 0.
          endif 
         
          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz
         
        endif
       
      enddo

      return
      end subroutine nlfsi_meshv_pert
!---------------------------------------------------------------------- 

      subroutine nlfsi_rstread

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'NLFSI'
      include 'CHKPOINT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      
      integer i,is
      character*1 exten0,exten1
      save exten0,exten1
      data exten0 /'0'/
      data exten1 /'1'/

      integer exten
      save exten
      data exten /0/

      character*132 fname
      character*1   fname1(132)
      equivalence (fname1,fname)      


      character*1 dot
      save dot
      data dot /'.'/
      character*16 outfmt

      integer ltrunc
      integer len

      real bd_psiv
      integer ilag
      integer iunit,ierr

      character*132 tmp

!     Do nothing if we are past the checkpoint files.        
      if (istep.ge.chkptnrsf) return

      ierr=0
!     Which set to read?        
!      exten=fsi_rst_fli
!     Using the value from the casename.restart
      if (istep.eq.0) then
        call blank(fname1,132)
        len=ltrunc(session,132)
        call chcopy(fname1,session,len)
        call chcopy(fname1(len+1),'.restart',8)
        if (nid.eq.0) then
          call IO_freeid(iunit,ierr)
          open(iunit,file=fname,status='old',action='read',
     $      iostat=ierr)
          if (ierr.eq.0) read(iunit,fmt=*,iostat=ierr) nlfsi_rst_fli
          close(iunit,iostat=ierr)
        endif
        call err_chk(ierr,'Error reading .restart file in fsi_rstread')
        call bcast(nlfsi_rst_fli, ISIZE)
      endif
      exten=nlfsi_rst_fli 

      if (istep.eq.0) then
!       Create filename
        call blank(fname1,132)
        call chcopy(fname1,nlfsi_rst_fname,13)
        call chcopy(fname1(14),dot,1)
        if (exten.eq.0) then
          call chcopy(fname1(15),exten0,1)
        else
          call chcopy(fname1(15),exten1,1)
        endif     ! exten.eq.0

!       Read restart file 
        if (nid.eq.0) then
          call IO_freeid(iunit,ierr)
          open(iunit, file=fname, status='old',action='read',
     $      iostat=ierr)
          if (ierr.eq.0) then
            write(6,'(A27,1x,A13)') 'NLFSI Restart: Opening file',
     $               fname
!            call blank(outfmt,16)  
!            write(outfmt,'(A1,I1,A14)') '(',fsi_nrst,'(E18.12E2,1x))'
!            read(1022,outfmt) (fsi_rst_psi(i),i=1,fsi_nrst)
!            read(1022,outfmt) (fsi_rst_psiv(i),i=1,fsi_nrst)
            read(iunit,*) (nlfsi_rst_psi(i),i=1,nlfsi_nrst)
            read(iunit,*) (nlfsi_rst_psiv(i),i=1,nlfsi_nrst)
            close(iunit)
            write(6,'(A4)')  'Done'
          endif   ! ierr.eq.0 
        endif     ! nid.eq.0
        call err_chk(ierr,'Error reading fsi_restart file.')

        call bcast(nlfsi_rst_psi, nlfsi_nrst*WDSIZE)
        call bcast(nlfsi_rst_psiv,nlfsi_nrst*WDSIZE)
      endif       ! istep.eq.0

      is     = istep+1
      if (is.gt.1) then
!       lag array for psi/psiv 
        do ilag=lorder-1,2,-1
          psilag(ilag)=psilag(ilag-1)
          psivlag(ilag)=psivlag(ilag-1)
        enddo
        psilag(1)=psi
        psivlag(1)=psiv
      endif  
      psi    = nlfsi_rst_psi(is) 
      psiv   = nlfsi_rst_psiv(is)

      if (istep.gt.0) then
        bd_psiv = psiv*bd(2)
        do ilag=2,NBD
          bd_psiv=bd_psiv+bd(ilag+1)*psivlag(ilag-1)
        enddo
        psiv_s=-nlfsi_inertia/DT*bd_psiv/
     $       (nlfsi_damp-nlfsi_inertia*bd(1)/DT)
      else
        psiv_s=psiv
      endif

!      if (nid.eq.0)
!     $   write(6,'(A6,1x,I10,1x,5(E15.8E2,1x))') 'Alpha:',istep,time,
!     $            0.000,psi,psiv,psiv_s

      if (nid.eq.0) then
!       Output psi,psiv,alpha            
        write(nlfsi_iunit1,'(I10,1x,4(E19.11E3,1x))') istep, 
     $       time,psi,psiv
        flush(nlfsi_iunit1) 
      endif             ! nid.eq.0


      return 
      end subroutine nlfsi_rstread

!---------------------------------------------------------------------- 
      subroutine nlfsi_map12(pm2,pm1)

      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'

      real pm1(lx1*ly1*lz1,lelv)
      real pm2(lx2*ly2*lz2,lelv)
      integer e

      do e=1,nelv
        call map12 (pm2(1,e),pm1(1,e),e)
      enddo
   
      return
      end
c-----------------------------------------------------------------------

