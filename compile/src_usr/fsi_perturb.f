!======================================================================
!     Routines for fluid structure interaction implicit solutions
!     Author: Prabal S. Negi
!     Description: Solver routines for the AMP Stokes' solve.
!
!     Based on Fischer P., Schmitt M. and Tomboulides A. (2017) Recent
!     Developments in Spectral Element Simulations of Moving-Domain
!     Problems.
!======================================================================       
c-----------------------------------------------------------------------

      subroutine fsi_main_pert

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
      include 'GEOM_DEF'
      include 'GEOM'
      include 'MVGEOM_DEF'
      include 'MVGEOM'
      include 'FSI'
      include 'CHKPOINT'
      include 'TIME_STEPPERD'
      include 'ADJOINT_DEF'
      include 'ADJOINT'
      

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

      integer i,lvp,lvp2
      parameter(lvp=lpx1*lpy1*lpz1*lpelv)
      parameter(lvp2=lpx2*lpy2*lpz2*lpelv)
      real wk_vxp(lvp),wk_vyp(lvp),wk_vzp(lvp),wk_prp(lvp2)

      real wk_vx0(lvp),wk_vy0(lvp),wk_vz0(lvp),wk_pr0(lvp2)
      real tmp_pr(lvp)

      integer jpold

      real*8 dnekclock          ! function. timing
      real*8 steptime           ! timing for step

      real rot_sn(3,3)
      real rot_cs(3,3)
      real dxyz(3),r0(3)
      real eta_sx,eta_sy,eta_sz


!     Start timing for this step   
      steptime = dnekclock()

      if (icalld.eq.0) then
!       Initialize timer            
        fsi_timea = 0.
        fsi_timei = 0.

        call fsi_init

        icalld=icalld+1
      endif  

      if (.not.IFFSI) return  

      if (TSTMODE.gt.0.and.ifchkptrst) then

!       do nothing. fsi_arnoldi_rstread is called within arn_rst_read
!       Need to verify that restarts are correct for Arnoldi 
        if (istep.eq.0) then
          if (ifusermv) then
            if (ifadj) then 
              call fsi_meshv_adj
            else
              call fsi_meshv_pert
            endif
          endif
          return
        endif

      else             

        if (ifchkptrst.and.istep.lt.chkptnrsf) then
!         Read saved values for restart
!         If the simulation has been restarted
          call fsi_rstread
          if (ifusermv) then
            if (ifadj) then 
              call fsi_meshv_adj
            else
              call fsi_meshv_pert
            endif
          endif

!         placeholder             
!          etaa = 0.
!          eta_s= eta

          call fsi_output

          fsi_timei=dnekclock()-steptime
          fsi_timea=fsi_timea+fsi_timei
          return

        elseif (istep.eq.0) then
!         If there are no restarts 
!         At first time step the projected velocity is
!         just the initialized velocity
          etav_s = etav_ini
          etav   = etav_ini
          eta    = eta_ini

          etaa   = 0.
          eta_s  = eta + etav_s*DT

          if (ifusermv) then
            if (ifadj) then 
              call fsi_meshv_adj
            else
              call fsi_meshv_pert
            endif
          endif

          call fsi_output

          fsi_timei=dnekclock()-steptime
          fsi_timea=fsi_timea+fsi_timei
          return
        endif

      endif       ! tstmode.gt.0.and.ifchkptrst


      x0(1) = fsi_x0
      x0(2) = fsi_y0
      if (ndim.eq.3) x0(3) = 1000.         ! set outside the domain

      ifdout=.false.
      iftout=.false.

!     Assuming perturbations are on jp=1
      scale=fsi_rescale
      call torque_calc_axis(scale,x0,ifdout,iftout,
     $      vxp(1,1),vyp(1,1),vzp(1,1),prp(1,1))

      call copy(fs_dragx (0),dragx (0),maxobj+1)
      call copy(fs_dragpx(0),dragpx(0),maxobj+1)
      call copy(fs_dragvx(0),dragpx(0),maxobj+1)
      call copy(fs_dragy (0),dragy (0),maxobj+1)
      call copy(fs_dragpy(0),dragpy(0),maxobj+1)
      call copy(fs_dragvy(0),dragvy(0),maxobj+1)
      call copy(fs_dragz (0),dragz (0),maxobj+1)
      call copy(fs_dragpz(0),dragpz(0),maxobj+1)
      call copy(fs_dragvz(0),dragvz(0),maxobj+1)
      call copy(fs_torqx (0),torqx (0),maxobj+1)
      call copy(fs_torqpx(0),torqpx(0),maxobj+1)
      call copy(fs_torqvx(0),torqpx(0),maxobj+1)
      call copy(fs_torqy (0),torqy (0),maxobj+1)
      call copy(fs_torqpy(0),torqpy(0),maxobj+1)
      call copy(fs_torqvy(0),torqvy(0),maxobj+1)
      call copy(fs_torqz (0),torqz (0),maxobj+1)
      call copy(fs_torqpz(0),torqpz(0),maxobj+1)
      call copy(fs_torqvz(0),torqvz(0),maxobj+1)

      if (ifadj) then
        call torque_calc_adj(scale,x0,ifdout,iftout,
     $           vxp(1,jp),vyp(1,jp),vzp(1,jp),prp(1,jp),
     $           fsi_grvx0,fsi_grvy0,fsi_grvz0)

!       Adjoint forces for eta^+ equation
        call copy(fsadj_dragx (0),dragx (0),maxobj+1)
        call copy(fsadj_dragpx(0),dragpx(0),maxobj+1)
        call copy(fsadj_dragvx(0),dragpx(0),maxobj+1)
        call copy(fsadj_dragy (0),dragy (0),maxobj+1)
        call copy(fsadj_dragpy(0),dragpy(0),maxobj+1)
        call copy(fsadj_dragvy(0),dragvy(0),maxobj+1)
        call copy(fsadj_dragz (0),dragz (0),maxobj+1)
        call copy(fsadj_dragpz(0),dragpz(0),maxobj+1)
        call copy(fsadj_dragvz(0),dragvz(0),maxobj+1)
        call copy(fsadj_torqx (0),torqx (0),maxobj+1)
        call copy(fsadj_torqpx(0),torqpx(0),maxobj+1)
        call copy(fsadj_torqvx(0),torqpx(0),maxobj+1)
        call copy(fsadj_torqy (0),torqy (0),maxobj+1)
        call copy(fsadj_torqpy(0),torqpy(0),maxobj+1)
        call copy(fsadj_torqvy(0),torqvy(0),maxobj+1)
        call copy(fsadj_torqz (0),torqz (0),maxobj+1)
        call copy(fsadj_torqpz(0),torqpz(0),maxobj+1)
        call copy(fsadj_torqvz(0),torqvz(0),maxobj+1)
      endif

!     Calculate Impulsive response (Green's function)      
      call fsi_plan3

!     AMP forces on the object
      call torque_calc_axis(scale,x0,ifdout,iftout,
     $      amp_vx,amp_vy,amp_vz,amp_pr)

      call copy(fg_dragx (0),dragx (0),maxobj+1)
      call copy(fg_dragpx(0),dragpx(0),maxobj+1)
      call copy(fg_dragvx(0),dragpx(0),maxobj+1)
      call copy(fg_dragy (0),dragy (0),maxobj+1)
      call copy(fg_dragpy(0),dragpy(0),maxobj+1)
      call copy(fg_dragvy(0),dragvy(0),maxobj+1)
      call copy(fg_dragz (0),dragz (0),maxobj+1)
      call copy(fg_dragpz(0),dragpz(0),maxobj+1)
      call copy(fg_dragvz(0),dragvz(0),maxobj+1)
      call copy(fg_torqx (0),torqx (0),maxobj+1)
      call copy(fg_torqpx(0),torqpx(0),maxobj+1)
      call copy(fg_torqvx(0),torqpx(0),maxobj+1)
      call copy(fg_torqy (0),torqy (0),maxobj+1)
      call copy(fg_torqpy(0),torqpy(0),maxobj+1)
      call copy(fg_torqvy(0),torqvy(0),maxobj+1)
      call copy(fg_torqz (0),torqz (0),maxobj+1)
      call copy(fg_torqpz(0),torqpz(0),maxobj+1)
      call copy(fg_torqvz(0),torqvz(0),maxobj+1)

      if (ifadj) then
        call torque_calc_adj(scale,x0,ifdout,iftout,
     $           amp_vx,amp_vy,amp_vz,amp_pr,
     $           fsi_grvx0,fsi_grvy0,fsi_grvz0)

!       Adjoint forces for eta^+ equation
        call copy(fgadj_dragx (0),dragx (0),maxobj+1)
        call copy(fgadj_dragpx(0),dragpx(0),maxobj+1)
        call copy(fgadj_dragvx(0),dragpx(0),maxobj+1)
        call copy(fgadj_dragy (0),dragy (0),maxobj+1)
        call copy(fgadj_dragpy(0),dragpy(0),maxobj+1)
        call copy(fgadj_dragvy(0),dragvy(0),maxobj+1)
        call copy(fgadj_dragz (0),dragz (0),maxobj+1)
        call copy(fgadj_dragpz(0),dragpz(0),maxobj+1)
        call copy(fgadj_dragvz(0),dragvz(0),maxobj+1)
        call copy(fgadj_torqx (0),torqx (0),maxobj+1)
        call copy(fgadj_torqpx(0),torqpx(0),maxobj+1)
        call copy(fgadj_torqvx(0),torqpx(0),maxobj+1)
        call copy(fgadj_torqy (0),torqy (0),maxobj+1)
        call copy(fgadj_torqpy(0),torqpy(0),maxobj+1)
        call copy(fgadj_torqvy(0),torqvy(0),maxobj+1)
        call copy(fgadj_torqz (0),torqz (0),maxobj+1)
        call copy(fgadj_torqpz(0),torqpz(0),maxobj+1)
        call copy(fgadj_torqvz(0),torqvz(0),maxobj+1)
      endif


      if (ifadj) then
        call fsi_calc_implicit_adj
        call fsi_meshv_adj
      else         
        call fsi_calc_implicit_pert
        call fsi_meshv_pert
      endif  

      fsi_timei=dnekclock()-steptime
      fsi_timea=fsi_timea+fsi_timei

      return
      end subroutine fsi_main_pert 
!---------------------------------------------------------------------- 
      subroutine fsi_meshv_pert

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'          ! xm1,ym1,zm1
      include 'MVGEOM_DEF'
      include 'MVGEOM'        ! wx,wy,wz
      include 'WING_MVMSH'
      include 'FSI'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'

      real ucx,ucy,ucz        ! mesh velocities
      real dx,dy,dz
      real eta_sx,eta_sy,eta_sz

      integer i,n
      real rot_sn(3,3)
      real rot_cs(3,3)
      real dxyz(3),r0(3)

      real bcx,bcy,bcz


      call rzero(rot_sn,9)            ! sine terms
      call rzero(rot_cs,9)            ! cosine terms
!     Matrix of sine terms
!     For a clockwise rotation for a positive eta
      if (fsi_ifrot) then
        rot_sn(1,2) =  eta_s
        rot_sn(2,1) = -eta_s
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
        if (.not.fsi_ifrot) then

!         No mesh motion for linear case
          wx(i,1,1,1) = 0. 
          wy(i,1,1,1) = 0.
          if (ndim.eq.3) wz(i,1,1,1) = 0. 

!         Boundary condition            
          eta_sx = 0.
          eta_sy = eta_s
          eta_sz = 0.  
            
          bcx = (fsi_grvx0(i,1,1,1,1)*eta_sx + 
     $           fsi_grvx0(i,1,1,1,2)*eta_sy +
     $           fsi_grvx0(i,1,1,1,3)*eta_sz)
          bcy = (fsi_grvy0(i,1,1,1,1)*eta_sx + 
     $           fsi_grvy0(i,1,1,1,2)*eta_sy +
     $           fsi_grvy0(i,1,1,1,3)*eta_sz)
          
          ucx = 0. - bcx
          ucy = etav_s - bcy
          if (if3d) then
            bcz = (fsi_grvz0(i,1,1,1,1)*eta_sx + 
     $             fsi_grvz0(i,1,1,1,2)*eta_sy +
     $             fsi_grvz0(i,1,1,1,3)*eta_sz)
            ucz = 0. - bcz
          else
            ucz = 0.
          endif 

!         Must ensure basev at the boundary == 1            
          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz

        else
!         Rotational FSI

!         Distance from Rotational axis                
          dx = xm1(i,1,1,1) - fsi_x0
          dy = ym1(i,1,1,1) - fsi_y0
          dz = 0.

!         Mesh motion 
          ucx =  0.
          ucy =  0.
          ucz =  0.
  
          wx(i,1,1,1) = 0. 
          wy(i,1,1,1) = 0.
          if (ndim.eq.3) wz(i,1,1,1) = 0. 

!!        Boundary condition            

!         Distance from axis   
          r0(1) = dx 
          r0(2) = dy
          r0(3) = dz

!         perturbation in position
          call mxm(rot_sn,3,r0,3,dxyz,1)
          eta_sx = dxyz(1)
          eta_sy = dxyz(2)
          eta_sz = dxyz(3) 

          bcx = (fsi_grvx0(i,1,1,1,1)*eta_sx + 
     $           fsi_grvx0(i,1,1,1,2)*eta_sy +
     $           fsi_grvx0(i,1,1,1,3)*eta_sz) 
          bcy = (fsi_grvy0(i,1,1,1,1)*eta_sx + 
     $           fsi_grvy0(i,1,1,1,2)*eta_sy +
     $           fsi_grvy0(i,1,1,1,3)*eta_sz)

          ucx =  etav_s*dy - bcx
          ucy = -etav_s*dx - bcy
          if (if3d) then
            bcz = (fsi_grvz0(i,1,1,1,1)*eta_sx + 
     $             fsi_grvz0(i,1,1,1,2)*eta_sy +
     $             fsi_grvz0(i,1,1,1,3)*eta_sz)
            ucz = 0. - bcz
          else
            ucz = 0.
          endif 

!         Must ensure basev at the boundary == 1            
          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz
        endif
       
      enddo

      return
      end subroutine fsi_meshv_pert
!---------------------------------------------------------------------- 

      subroutine fsi_calc_implicit_pert

!     Calculate superposition of Solutions to get implicit update

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
      include 'FSI'

      real uy                 ! output vertical velocity

      integer icalld
      save icalld
      data icalld /0/

      real ext_k,bd_eta,bd_etav
      real ab0,ab1,ab2
      integer iobj
      integer ilag
      real const
      real solid_area

      real lhs,rhs,rhs2
!      real alpha
      real Fs,Fb,Fk_ext,Fg,Fdx,Ft

      integer irst,NBDMSH,NABMSH
      integer istep_old

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
      if (.not.fsi_ifrot) then
        Fs = fs_dragy(iobj) 
        Fg = fg_dragy(iobj)     ! Green's function
      else  
!       Rotational coordinate system is opposite to the one used in Nek            
        Fs = -fs_torqz(iobj) 
        Fg = -fg_torqz(iobj)    ! Green's function
      endif                            

!     Position
!     Copy lag array for eta 
      do ilag=lorder-1,2,-1
        etalag(ilag)=etalag(ilag-1)
      enddo
      etalag(1)=eta
      eta  = eta + (ABMSH(1)*etav + ABMSH(2)*etavlag(1)
     $       + ABMSH(3)*etavlag(2))*DT

      Fk_ext = -fsi_stiff*eta

!     Backward differentiation
      bd_etav = etav*bd(2)
      do ilag=2,NBD
        bd_etav = bd_etav + bd(ilag+1)*etavlag(ilag-1)
      enddo

!     Copy lag array for etav 
      do ilag=lorder-1,2,-1
        etavlag(ilag)=etavlag(ilag-1)
      enddo
      etavlag(1)=etav

!     Backward differentiation coefficients already assume the terms
!     have been moved to the rhs. So a -ve sign is already included.      
      rhs = Fs + Fk_ext - fsi_inertia/DT*bd(1)*etav_s
     $         + fsi_inertia/DT*bd_etav + fsi_damp*etav_s
      lhs = fsi_inertia/DT*bd(1)*etav_g - Fg -fsi_damp*etav_g
      fsi_alpha = rhs/lhs               ! Eqn (49) in paper.
      
      etav = etav_s + fsi_alpha*etav_g

      etaa = (bd(1)*etav - bd_etav)/DT
     
!     This assumes jp=1      
      call opadd2cm(vxp,vyp,vzp,amp_vx,amp_vy,amp_vz,fsi_alpha)
      call add2s2(prp,amp_pr,fsi_alpha,nx2*ny2*nz2*nelt)

!     Get approximate velocity for the next time-step.
!     Calculation is independent since we correct this in the Stokes' step anyway.
!     Good projection allows the correction to be small.
!     According to paper:
      bd_etav=etav*bd(2)
      do ilag=2,NBD
        bd_etav=bd_etav+bd(ilag+1)*etavlag(ilag-1)
      enddo
      etav_s=-fsi_inertia/DT*bd_etav/(fsi_damp-fsi_inertia*bd(1)/DT)

!     Find position for next step.
!     Here I'm assuming that DT is constant.
!     Otherwise we need to update dtlag array as well.
!     And then revert it back for the normal nek solve.      
      if (ifpert) then
         istep_old = istep
         istep = istep+1
         NBDMSH = 1
         NABMSH = PARAM(28)
         IF (NABMSH.GT.ISTEP .and. irst.le.0) NABMSH = ISTEP
         IF (IFSURT)          NABMSH = NBD
         CALL RZERO   (ABMSH,10)
         CALL SETABBD (ABMSH,DTLAG,NABMSH,NBDMSH)

!        This is the position for the next step.         
         eta_s  = eta + (ABMSH(1)*etav + ABMSH(2)*etavlag(1)
     $       + ABMSH(3)*etavlag(2))*DT
         istep = istep_old
      endif 

!     All forces      
      fsi_Fs = Fs

      fsi_Fg = Fg

      fsi_Ff = fsi_Fs + fsi_alpha*fsi_Fg

      fsi_Fd = fsi_damp*etav

      fsi_Fk = Fk_ext

      fsi_Ft = fsi_Ff+fsi_Fd+fsi_Fk

      call fsi_output

      call fsi_rstsave 

      return
      end subroutine fsi_calc_implicit_pert

!---------------------------------------------------------------------- 

      subroutine fsi_calc_implicit_adj

!     Calculate superposition of Solutions to get implicit update

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
      include 'FSI'

      real uy                 ! output vertical velocity

      integer icalld
      save icalld
      data icalld /0/

      real ext_k,bd_eta,bd_etav
      real ab0,ab1,ab2
      integer iobj
      integer ilag
      real const
      real solid_area

      real lhs,rhs,rhs2
      real Fs,Fb,Fg,Fdx,Ft
      real FsAdj,FgAdj


!     Only doing for one object right now      
      iobj=0
      if (.not.fsi_ifrot) then
!       Forces for etav^+            
        Fs    = fs_dragy(iobj) 
        Fg    = fg_dragy(iobj)        ! Green's function
!       Forces for eta^+        
        FsAdj = fsadj_dragy(iobj) 
        FgAdj = fgadj_dragy(iobj)     ! Green's function
      else  
!       Rotational coordinate system is opposite to the one used in Nek
!       Forces for etav^+            
        Fs    = -fs_torqz(iobj) 
        Fg    = -fg_torqz(iobj)       ! Green's function
!       Forces for eta^+        
        FsAdj = -fsadj_dragy(iobj) 
        FgAdj = -fgadj_dragy(iobj)    ! Green's function
      endif                            

!     Backward differentiation etav
      bd_etav = etav*bd(2)
      do ilag=2,NBD
        bd_etav = bd_etav + bd(ilag+1)*etavlag(ilag-1)
      enddo

!     Backward differentiation eta
      bd_eta = eta*bd(2)
      do ilag=2,NBD
        bd_eta = bd_eta + bd(ilag+1)*etalag(ilag-1)
      enddo

!     Copy lag array for etav 
      do ilag=lorder-1,2,-1
        etavlag(ilag)=etavlag(ilag-1)
      enddo
      etavlag(1)=etav

!     Copy lag array for eta 
      do ilag=lorder-1,2,-1
        etalag(ilag)=etalag(ilag-1)
      enddo
      etalag(1)=eta

!     Backward differentiation coefficients already assume the terms
!     have been moved to the rhs. So a -ve sign is already included.      
      rhs =  (Fs + fsi_inertia/DT*bd_etav)
     $     + (FsAdj*DT/bd(1) + bd_eta/bd(1))
     $     - (fsi_inertia*bd(1)/DT*etav_s - fsi_damp*etav_s 
     $        +fsi_stiff*DT/bd(1)*etav_s)

      lhs = (fsi_inertia/DT*bd(1)*etav_g - fsi_damp*etav_g
     $       + fsi_stiff*DT/bd(1)*etav_g - Fg - FgAdj*DT/bd(1))
      fsi_alpha = rhs/lhs
      
      etav = etav_s + fsi_alpha*etav_g
      eta  = DT/bd(1)*(FsAdj + fsi_alpha*FgAdj
     $                 + bd_eta/DT - fsi_stiff*etav)

      etaa = (bd(1)*etav - bd_etav)/DT
     
!     This assumes jp=1      
      call opadd2cm(vxp,vyp,vzp,amp_vx,amp_vy,amp_vz,fsi_alpha)
      call add2s2(prp,amp_pr,fsi_alpha,nx2*ny2*nz2*nelt)

!     Get approximate velocity for the next time-step.
!     Test case when holding acceleration constant.
!     Corrected later anyway      
      bd_etav=etav*bd(2)
      do ilag=2,NBD
        bd_etav=bd_etav+bd(ilag+1)*etavlag(ilag-1)
      enddo
      etav_s=DT/bd(1)*(etaa+bd_etav/DT)

!     All forces      
      fsi_Fs = Fs

      fsi_Fg = Fg

      fsi_Ff = fsi_Fs + fsi_alpha*fsi_Fg

      fsi_Fd = FgAdj

      fsi_Fk = FsAdj

      fsi_Ft = fsi_Ff

      call fsi_output

      call fsi_rstsave 

      return
      end subroutine fsi_calc_implicit_adj

!---------------------------------------------------------------------- 
      subroutine fsi_meshv_adj

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'          ! xm1,ym1,zm1
      include 'MVGEOM_DEF'
      include 'MVGEOM'        ! wx,wy,wz
      include 'WING_MVMSH'
      include 'FSI'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'

      real ucx,ucy,ucz        ! mesh velocities
      real dx,dy,dz

      integer i,n


      n = nx1*ny1*nz1*nelv

      do i=1,n                          ! Translational velocity
        if (.not.fsi_ifrot) then

!         No mesh motion for linear case
          wx(i,1,1,1) = 0. 
          wy(i,1,1,1) = 0.
          if (ndim.eq.3) wz(i,1,1,1) = 0. 

!         Boundary condition
          ucx = 0.
          ucy = etav_s
          ucz = 0.

!         Must ensure basev at the boundary == 1            
          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz

        else
!         Rotational FSI Adjoint

!         Needs to be derived              
          if (nid.eq.0) write(6,*) 'Rotational Adjoint not implemented'
          call exitt
  
        endif
       
      enddo

      return
      end subroutine fsi_meshv_adj
!---------------------------------------------------------------------- 

