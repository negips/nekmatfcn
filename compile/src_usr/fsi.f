!======================================================================
!     Routines for fluid structure interaction implicit solutions
!     Author: Prabal S. Negi
!     Based on Fischer P., Schmitt M. and Tomboulides A. (2017) Recent
!     Developments in Spectral Element Simulations of Moving-Domain
!     Problems.
!     Big thanks to Paul Fischer for his help.
!
!======================================================================       
!     read parameters fluid-structure interaction 
      subroutine fsi_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /FSI/ iffsi,eta_ini,etav_ini,fsi_stiff,fsi_damp,
     $               fsi_inertia,fsi_x0,fsi_y0,fsi_rescale,fsi_rst_fli,
     $               fsi_iftermso 

!     default values
      iffsi             = .FALSE.   ! if FSI
      eta_ini           =  0.       ! initial position
      etav_ini          =  0.       ! initial velocity
      fsi_stiff         =  1.       ! structural stiffness
      fsi_damp          = -1.       ! structural damping
      fsi_inertia       =  1.0      ! inertia
      fsi_x0            =  0.       ! elastic axis x0
      fsi_y0            =  0.       ! elastic axis y0
      fsi_rst_fli       =  0        ! restart file no.
      fsi_rescale       =  1.0      ! scale the fluid forces
      fsi_iftermso      = .FALSE.   ! if output terms for debugging

!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=FSI,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading FSI parameters.$')

!     broadcast data
      call bcast(iffsi        , LSIZE)
      call bcast(eta_ini      ,WDSIZE)
      call bcast(etav_ini     ,WDSIZE)
      call bcast(fsi_stiff    ,WDSIZE)
      call bcast(fsi_damp     ,WDSIZE)
      call bcast(fsi_inertia  ,WDSIZE)
      call bcast(fsi_x0       ,WDSIZE)
      call bcast(fsi_y0       ,WDSIZE)
      call bcast(fsi_rst_fli  , ISIZE)
      call bcast(fsi_rescale  ,WDSIZE)
      call bcast(fsi_iftermso , LSIZE)

      return
      end
!----------------------------------------------------------------------
!     write parameters fluid-structure interaction
      subroutine fsi_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'FSI'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

      logical ifrot

!     namelists
      namelist /FSI/ iffsi,eta_ini,etav_ini,fsi_stiff,fsi_damp,
     $               fsi_inertia,fsi_x0,fsi_y0,fsi_rescale,fsi_rst_fli,
     $               fsi_iftermso,ifrot

      ifrot = fsi_ifrot

!     write to log
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=FSI,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing FSI parameters.$')

      return
      end
!----------------------------------------------------------------------
!====================================================================== 
      subroutine fsi_main

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
      include 'MVGEOM_DEF'
      include 'MVGEOM'
      include 'FSI'
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

      logical add_displ

      real*8 dnekclock          ! function. timing
      real*8 steptime           ! timing for step

      if (ifpert) then
        call fsi_main_pert
        return
      endif  

      if (icalld.eq.0) then
!       Initialize timer 
        fsi_timea = 0.
        fsi_timei = 0.
           
        call fsi_init
        
        icalld=icalld+1
      endif  

      if (.not.IFFSI) then
        etav_s = 0.
        etav   = 0.
        eta_s  = 0.
        if (ifusermv) call fsi_meshv
        return
      endif  

!     Start timing for this step   
      steptime = dnekclock()

      if (ifchkptrst.and.istep.lt.chkptnrsf) then
!       Read saved values for restart
!       If the simulation has been restarted
        call fsi_rstread

        if (ifusermv) call fsi_meshv

        call fsi_output

        fsi_timei=dnekclock()-steptime
        fsi_timea=fsi_timea+fsi_timei
        return

      elseif (istep.eq.0) then
!       If there are no restarts              
!       At first time step the projected velocity is
!       just the initialized velocity
        etav_s = etav_ini
        etav   = etav_ini
        eta    = eta_ini
        etaa   = 0.
        eta_s  = etav_s*DT            ! Should probably change this

        if (ifusermv) call fsi_meshv

        call fsi_output

        fsi_timei=dnekclock()-steptime
        fsi_timea=fsi_timea+fsi_timei
        return
      endif

      x0(1) = fsi_x0
      x0(2) = fsi_y0
      if (ndim.eq.3) x0(3) = 1000.         ! set outside the domain

      ifdout=.false.
      iftout=.false.

      scale=fsi_rescale
      call torque_calc_axis(scale,x0,ifdout,iftout,vx,vy,vz,pr)

!     Who came up with this code? 
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

!     Calculate Impulsive response (Green's function)      
      call fsi_plan3

!     AMP forces on the object
      call torque_calc_axis(scale,x0,ifdout,iftout,amp_vx,amp_vy,
     $                       amp_vz,amp_pr)

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

      call fsi_calc_implicit_soln

      if (ifusermv) call fsi_meshv

      fsi_timei=dnekclock()-steptime
      fsi_timea=fsi_timea+fsi_timei

      return
      end subroutine fsi_main      
!---------------------------------------------------------------------- 

      subroutine fsi_calc_implicit_soln

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
      real Fs,Fb,Fk_ext,Fg,Fdx,Fd,Ft,F0
      real Fk_paul,y_del

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
        F0 = fs0_dragy(iobj)  
      else  
!       Rotational coordinate system is opposite to the one used in Nek            
        Fs = -fs_torqz(iobj) 
        Fg = -fg_torqz(iobj)    ! Green's function
        F0 = -fs0_torqz(iobj)  
      endif

!     These are zero in the non-linear case

      if (fsi_remove_base) then
        continue
      else
        F0 = 0.
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
!      rhs = Fs + 0*Fk_ext - fsi_inertia/DT*bd(1)*etav_s
!     $         + Fk_paul  + fsi_inertia/DT*bd_etav + fsi_damp*etav_s
      rhs = Fs + Fk_ext - fsi_inertia/DT*bd(1)*etav_s
     $         + fsi_inertia/DT*bd_etav + fsi_damp*etav_s - F0
      lhs = fsi_inertia/DT*bd(1)*etav_g - Fg -fsi_damp*etav_g
      fsi_alpha = rhs/lhs               ! Eqn (49) in paper.

      etav = etav_s + fsi_alpha*etav_g

      etaa = (bd(1)*etav - bd_etav)/DT
      
      call opadd2cm(vx,vy,vz,amp_vx,amp_vy,amp_vz,fsi_alpha)
      call add2s2(pr,amp_pr,fsi_alpha,nx2*ny2*nz2*nelt)

!     Get approximate velocity for the next time-step.
!     Calculation is independent since we correct this in the Stokes' step anyway.
!     Good projection allows the correction to be small.
!     According to paper:
      bd_etav = etav*bd(2)
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

      fsi_Fb = Fb 

      fsi_Fg = Fg

      fsi_Ff = fsi_Fs + fsi_Fb + fsi_alpha*fsi_Fg

      fsi_Fd = fsi_damp*etav

!      fsi_Fdx = Fdx

      fsi_Fk = Fk_ext

      fsi_Ft = fsi_Ff+fsi_Fd+fsi_Fk-F0

      call fsi_output

      call fsi_rstsave 

      return
      end subroutine fsi_calc_implicit_soln

!---------------------------------------------------------------------- 

      subroutine fsi_rstsave

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'FSI'
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
        fsi_rst_eta(i)  = eta
        fsi_rst_etav(i) = etav
      endif

      if (nid.eq.0) then
        if (is.eq.chkptnrsf-1) then
!         Create filename
          call blank(fname1,132)
          call chcopy(fname1,fsi_rst_fname,11)
          call chcopy(fname1(12),dot,1)
          if (exten.eq.0) then
            call chcopy(fname1(13),exten0,1)
            exten=1  
          else
            call chcopy(fname1(13),exten1,1)
            exten=0  
          endif     ! exten.eq.0

!         Write output
          call IO_freeid(iunit,ierr)     
          open(iunit,file=fname,status='unknown',action='write',
     $      iostat=ierr)
          if (ierr.eq.0) then
            write(6,'(A24,1x,A13)') 'FSI Restart: Saving file',
     $               fname
!            call blank(outfmt,16)  
!            write(outfmt,'(A1,I1,A14)') '(',fsi_nrst,'(E18.12E2,1x))'
!            write(1021,outfmt) (fsi_rst_eta(i),i=1,fsi_nrst)
!            write(1021,outfmt) (fsi_rst_etav(i),i=1,fsi_nrst)
            write(iunit,*) (fsi_rst_eta(i),i=1,fsi_nrst)
            write(iunit,*) (fsi_rst_etav(i),i=1,fsi_nrst)
            close(iunit)
            write(6,'(A4)')  'Done'
          endif   ! ierr.eq.0 
        endif     ! is.eq.chkptnrst-1
      endif       ! nid.eq.0

      call err_chk(ierr,'Error opening fsi_restart file.')
      
      return 
      end subroutine fsi_rstsave

!---------------------------------------------------------------------- 
      subroutine fsi_meshv

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
      real dx,dy

      integer i,n


      if (ifpert) then
        call fsi_meshv_pert
        return
      endif  

      n = nx1*ny1*nz1*nelv

      do i=1,n                          ! Translational velocity
!       Current time step mesh velocity at the wall is the same
!       as the velocity at the wall
        if (.not.fsi_ifrot) then

          ucx =  0.
          ucy =  etav
          ucz =  0.0
  
          wx(i,1,1,1) = basev(i)*ucx
          wy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) wz(i,1,1,1) = basev(i)*ucz

          ucx =  0.
          ucy =  etav_s
          ucz =  0.

          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz

        else
          dx = xm1(i,1,1,1) - fsi_x0
          dy = ym1(i,1,1,1) - fsi_y0
          
          ucx =  etav*dy
          ucy =  -etav*dx
          if (ndim.eq.3) ucz =  0.0
  
          wx(i,1,1,1) = basev(i)*ucx
          wy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) wz(i,1,1,1) = basev(i)*ucz

          ucx =   etav_s*dy
          ucy =  -etav_s*dx
          if (ndim.eq.3) ucz =  0.0
         
          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz
         
        endif
!       This is the prescribed wall velocity for the next step
!       (userbc)
       
      enddo

      return
      end subroutine fsi_meshv
!---------------------------------------------------------------------- 

      subroutine fsi_rstread

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'FSI'
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

      real bd_etav
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
          if (ierr.eq.0) read(iunit,fmt=*,iostat=ierr) fsi_rst_fli
          close(iunit,iostat=ierr)
        endif
        call err_chk(ierr,'Error reading .restart file in fsi_rstread')
        call bcast(fsi_rst_fli, ISIZE)
      endif
      exten=fsi_rst_fli       

      if (istep.eq.0) then
!       Create filename
        call blank(fname1,132)
        call chcopy(fname1,fsi_rst_fname,11)
        call chcopy(fname1(12),dot,1)
        if (exten.eq.0) then
          call chcopy(fname1(13),exten0,1)
        else
          call chcopy(fname1(13),exten1,1)
        endif     ! exten.eq.0

!       Read restart file 
        if (nid.eq.0) then
          call IO_freeid(iunit,ierr)
          open(iunit, file=fname, status='old',action='read',
     $      iostat=ierr)
          if (ierr.eq.0) then
            write(6,'(A25,1x,A13)') 'FSI Restart: Opening file',
     $               fname
!            call blank(outfmt,16)  
!            write(outfmt,'(A1,I1,A14)') '(',fsi_nrst,'(E18.12E2,1x))'
!            read(1022,outfmt) (fsi_rst_eta(i),i=1,fsi_nrst)
!            read(1022,outfmt) (fsi_rst_etav(i),i=1,fsi_nrst)
            read(iunit,*) (fsi_rst_eta(i),i=1,fsi_nrst)
            read(iunit,*) (fsi_rst_etav(i),i=1,fsi_nrst)
            close(iunit)
            write(6,'(A4)')  'Done'
          endif   ! ierr.eq.0 
        endif     ! nid.eq.0
        call err_chk(ierr,'Error reading fsi_restart file.')

        call bcast(fsi_rst_eta, fsi_nrst*WDSIZE)
        call bcast(fsi_rst_etav,fsi_nrst*WDSIZE)
      endif       ! istep.eq.0

      is     = istep+1
      if (is.gt.1) then
!       lag array for eta/etav 
        do ilag=lorder-1,2,-1
          etalag(ilag)=etalag(ilag-1)
          etavlag(ilag)=etavlag(ilag-1)
        enddo
        etalag(1)=eta
        etavlag(1)=etav
      endif  
      eta    = fsi_rst_eta(is) 
      etav   = fsi_rst_etav(is)

      if (istep.gt.0) then
        bd_etav = etav*bd(2)
        do ilag=2,NBD
          bd_etav=bd_etav+bd(ilag+1)*etavlag(ilag-1)
        enddo
        etav_s =-fsi_inertia/DT*bd_etav/(fsi_damp-fsi_inertia*bd(1)/DT)
      else
        etav_s = etav
        etaa   = 0.
      endif


      return 
      end subroutine fsi_rstread

!---------------------------------------------------------------------- 

      SUBROUTINE CHKDIV (UX,UY,UZ)

      implicit none            

      INCLUDE 'SIZE_DEF'
      INCLUDE 'SIZE'
      INCLUDE 'MASS_DEF'
      INCLUDE 'MASS'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'

      REAL DIVFLD,WORK
      COMMON /SCRMG/ DIVFLD (LX2,LY2,LZ2,LELV)
     $ ,             WORK   (LX2,LY2,LZ2,LELV)
      INTEGER NTOT2
      REAL GLSUM
      REAL DIVV

      REAL UX(LX1,LY1,LZ1,LELT),UY(LX1,LY1,LZ1,LELT),
     $     UZ(LX1,LZ1,LZ1,LELT)  

      NTOT2 = NX2*NY2*NZ2*NELV
      CALL OPDIV   (DIVFLD,UX,UY,UZ)
      CALL COL3    (WORK,DIVFLD,BM2INV,NTOT2)
      CALL COL2    (WORK,DIVFLD,NTOT2)
      DIVV  = SQRT(GLSUM(WORK,NTOT2)/VOLVM2)
C
      if (nid.eq.0) 
     $  write(6,'(A17,2x,E15.8E2)') 'Soln. Divergence:', DIVV

      RETURN
      END
!----------------------------------------------------------------------            

      subroutine fsi_init

      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'

      integer lt2
      parameter (lt2=lx2*ly2*lz2*lelv)

      real scale
      real x0(3)
      logical ifdout,iftout
      integer ierr

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

      integer iel,ifld,iface,nfaces
      character cb*3

!     Initialize timer            
      fsi_timea = 0.
      fsi_timei = 0.

!      Done in userchk
!      call set_obj
!      Put a check for number of objects. Prabal

!     This is ridiculous    
      call rzero(fs_dragx (0),maxobj+1)
      call rzero(fs_dragpx(0),maxobj+1)
      call rzero(fs_dragvx(0),maxobj+1)
      call rzero(fs_dragy (0),maxobj+1)
      call rzero(fs_dragpy(0),maxobj+1)
      call rzero(fs_dragvy(0),maxobj+1)
      call rzero(fs_dragz (0),maxobj+1)
      call rzero(fs_dragpz(0),maxobj+1)
      call rzero(fs_dragvz(0),maxobj+1)
      call rzero(fs_torqx (0),maxobj+1)
      call rzero(fs_torqpx(0),maxobj+1)
      call rzero(fs_torqvx(0),maxobj+1)
      call rzero(fs_torqy (0),maxobj+1)
      call rzero(fs_torqpy(0),maxobj+1)
      call rzero(fs_torqvy(0),maxobj+1)
      call rzero(fs_torqz (0),maxobj+1)
      call rzero(fs_torqpz(0),maxobj+1)
      call rzero(fs_torqvz(0),maxobj+1)
      
!     More ridiculousness        
      call rzero(fg_dragx (0),maxobj+1)
      call rzero(fg_dragpx(0),maxobj+1)
      call rzero(fg_dragvx(0),maxobj+1)
      call rzero(fg_dragy (0),maxobj+1)
      call rzero(fg_dragpy(0),maxobj+1)
      call rzero(fg_dragvy(0),maxobj+1)
      call rzero(fg_dragz (0),maxobj+1)
      call rzero(fg_dragpz(0),maxobj+1)
      call rzero(fg_dragvz(0),maxobj+1)
      call rzero(fg_torqx (0),maxobj+1)
      call rzero(fg_torqpx(0),maxobj+1)
      call rzero(fg_torqvx(0),maxobj+1)
      call rzero(fg_torqy (0),maxobj+1)
      call rzero(fg_torqpy(0),maxobj+1)
      call rzero(fg_torqvy(0),maxobj+1)
      call rzero(fg_torqz (0),maxobj+1)
      call rzero(fg_torqpz(0),maxobj+1)
      call rzero(fg_torqvz(0),maxobj+1)

!     ...        
      call rzero(fs0_dragx (0),maxobj+1)
      call rzero(fs0_dragpx(0),maxobj+1)
      call rzero(fs0_dragvx(0),maxobj+1)
      call rzero(fs0_dragy (0),maxobj+1)
      call rzero(fs0_dragpy(0),maxobj+1)
      call rzero(fs0_dragvy(0),maxobj+1)
      call rzero(fs0_dragz (0),maxobj+1)
      call rzero(fs0_dragpz(0),maxobj+1)
      call rzero(fs0_dragvz(0),maxobj+1)
      call rzero(fs0_torqx (0),maxobj+1)
      call rzero(fs0_torqpx(0),maxobj+1)
      call rzero(fs0_torqvx(0),maxobj+1)
      call rzero(fs0_torqy (0),maxobj+1)
      call rzero(fs0_torqpy(0),maxobj+1)
      call rzero(fs0_torqvy(0),maxobj+1)
      call rzero(fs0_torqz (0),maxobj+1)
      call rzero(fs0_torqpz(0),maxobj+1)
      call rzero(fs0_torqvz(0),maxobj+1)

      call rzero(etalag,lorder-1)
      call rzero(etavlag,lorder-1)

      if (ifpert) then
        call gradm1(fsi_grvx0(1,1,1,1,1),fsi_grvx0(1,1,1,1,2),
     $              fsi_grvx0(1,1,1,1,3),vx)
        call gradm1(fsi_grvy0(1,1,1,1,1),fsi_grvy0(1,1,1,1,2),
     $              fsi_grvy0(1,1,1,1,3),vy)
        if (if3d) then
          call gradm1(fsi_grvz0(1,1,1,1,1),fsi_grvz0(1,1,1,1,2),
     $                fsi_grvz0(1,1,1,1,3),vz)
        endif

      endif  


      if (nid.eq.0) then

!       Open file to output eta,etav,alpha               
        call IO_freeid(fsi_iunit1,ierr)
        open(fsi_iunit1,file=fsi_fname1,status='unknown',
     $       action='write',iostat=ierr)
        
        write(fsi_iunit1,'(A10,1x,7(A19,1x))') 'istep',
     $           'time','alpha','eta','etav','etaa','eta_s','etav_s'
        flush(fsi_iunit1)

!        Open file to output terms of the sturctural solve
!        Mostly for debugging.
         if (fsi_iftermso) then          
           call IO_freeid(fsi_iunit2,ierr)
           open(fsi_iunit2,file=fsi_fname2,status='unknown',
     $          action='write',iostat=ierr)
           write(fsi_iunit2,'(A10,1x,8(A16,1x))') 'ist','Time',
     $         'Fs','Fg','alpha','Ff','Fk','Fd','Ft'
           flush(fsi_iunit2)
         endif  
      endif

!     Forces due to Baseflow state
      x0(1) = fsi_x0
      x0(2) = fsi_y0
      if (ndim.eq.3) x0(3) = 1000.         ! set outside the domain

      ifdout=.true.
      iftout=.true.

      scale=fsi_rescale
      call torque_calc_axis(scale,x0,ifdout,iftout,vx,vy,vz,pr)

      call copy(fs0_dragx (0),dragx (0),maxobj+1)
      call copy(fs0_dragpx(0),dragpx(0),maxobj+1)
      call copy(fs0_dragvx(0),dragpx(0),maxobj+1)
      call copy(fs0_dragy (0),dragy (0),maxobj+1)
      call copy(fs0_dragpy(0),dragpy(0),maxobj+1)
      call copy(fs0_dragvy(0),dragvy(0),maxobj+1)
      call copy(fs0_dragz (0),dragz (0),maxobj+1)
      call copy(fs0_dragpz(0),dragpz(0),maxobj+1)
      call copy(fs0_dragvz(0),dragvz(0),maxobj+1)
      call copy(fs0_torqx (0),torqx (0),maxobj+1)
      call copy(fs0_torqpx(0),torqpx(0),maxobj+1)
      call copy(fs0_torqvx(0),torqpx(0),maxobj+1)
      call copy(fs0_torqy (0),torqy (0),maxobj+1)
      call copy(fs0_torqpy(0),torqpy(0),maxobj+1)
      call copy(fs0_torqvy(0),torqvy(0),maxobj+1)
      call copy(fs0_torqz (0),torqz (0),maxobj+1)
      call copy(fs0_torqpz(0),torqpz(0),maxobj+1)
      call copy(fs0_torqvz(0),torqvz(0),maxobj+1)

      etav_g = 1.0           ! currently hard-coded

!     Initialize Added-Mass-Partition arrays
!     Set masks for the AMP solve
      call set_amp_mask

      call opzero(amp_vx,amp_vy,amp_vz)
      call rzero(amp_pr,nx2*ny2*nz2*nelt)

!     Build a mask to remove only 'v  '
!     and preserving the 'mv ' points.
      if (ifpert) then
        call rone(fsi_mask,lx1*ly1*lz1*lelv)
        ifld = 1
        nfaces=2*ndim
        do iel=1,nelv
          do iface = 1,nfaces
            cb = cbc(iface,iel,ifld)
            if (cb.eq.'v  ') then
!             Put zeros on this face
              call facev(fsi_mask,iel,iface,0.,nx1,ny1,nz1)
            endif
          enddo ! iel
        enddo   ! iface

!       Weights for inner-product in Arnoldi        
        fsi_eta_wt = 1.0
        fsi_etav_wt = fsi_inertia
      endif     ! ifpert


      return 
      end subroutine fsi_init
!----------------------------------------------------------------------       

      subroutine fsi_output

      implicit none

      INCLUDE 'SIZE_DEF'
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'
      INCLUDE 'FSI'

      if (nid.eq.0) then
!       Output eta,etav,alpha            
        write(fsi_iunit1,'(i10,1x,7(1pe19.11e3,1x))') istep, 
     $       time,fsi_alpha,eta,etav,etaa,eta_s,etav_s
        flush(fsi_iunit1) 

!       If forces output is on
        if (fsi_iftermso) then
          write(fsi_iunit2,'(I10,1x,8(1pE16.7E3,1x))') istep, 
     $       time,fsi_Fs,fsi_Fg,fsi_alpha,fsi_Ff,
     $       fsi_Fk,fsi_Fd,fsi_Ft
          flush(fsi_iunit2)
        endif           ! iftermso 
      endif             ! nid.eq.0

      return
      end subroutine fsi_output            
!----------------------------------------------------------------------       



