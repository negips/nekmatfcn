!=======================================================================
!     Adam Peplinski; 2015.12.02
!     Set of subroutines to use time steppers for power iterations or
!     solution of eigenvalue problem with Arnoldi algorithm
!
!     Parameters used by this set of subroutines:
!     TIME_STEPPERD:
!     TSTMODE - 0 - no time stepper
!               1 - direct mode
!               2 - adjoint mode
!               3 - initial optimal condition
!     TSTSTEP - frequency of calling stepper_vsolve (number of time steps)
!     TSTCMAX - maximal number of stepper cycles (meaning of this parameter
!               depends on the applied package; arnoldi vs power iterations)
!     TSTTOL - convergence tolerance (e.g. PARPACK tolerance for Ritz values)
!     TSTIFUZAWA - first step run with UZAWA
!
!     For conjugated heat transfer
!     CONHT:
!     CHCST_SC, CHCFF_V, CHCFF_T - velocity and temperature scaling
!                                  factors
!     CHGRX, CHGRY, CHGRZ - gravitational acceleration
!     CHRA, CHPR - Rayleight and Prandtl numbers
!=======================================================================
!***********************************************************************
!     read parameters for time stepper
      subroutine tstpr_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'TIME_STEPPERD'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr, len
      integer itmp(3)

!     namelists
      namelist /TSTEPPER/ TSTMODE, TSTSTEP, TSTCMAX,
     $     TSTTOL, TSTIFUZAWA, TSTIFPR
!-----------------------------------------------------------------------
!     default values
      TSTMODE = 1
      TSTSTEP = 40
      TSTCMAX = 40
      TSTTOL = 1.0d-6
      TSTIFUZAWA = .TRUE.
      TSTIFPR    = .TRUE.
!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=TSTEPPER,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading TSTEPPER parameters.$')

!     broadcast data
      if (NID.eq.0) then
         itmp(1) = TSTMODE
         itmp(2) = TSTSTEP
         itmp(3) = TSTCMAX
      endif
      len = 3*ISIZE
      call bcast(itmp,len)
      if (NID.ne.0) then
         TSTMODE = itmp(1)
         TSTSTEP = itmp(2)
         TSTCMAX = itmp(3)
      endif

      call bcast(TSTTOL,WDSIZE)
      call bcast(TSTIFUZAWA,LSIZE)
      call bcast(TSTIFPR,LSIZE)

      return
      end
!***********************************************************************
!     write parameters for time stepper
      subroutine tstpr_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'TIME_STEPPERD'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /TSTEPPER/ TSTMODE, TSTSTEP, TSTCMAX,
     $     TSTTOL, TSTIFUZAWA,TSTIFPR
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=TSTEPPER,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing TSTEPPER parameters.$')

      return
      end
!***********************************************************************
!     main time stepper interface
      subroutine tstepper
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO, NDIM, NPERT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP
!-----------------------------------------------------------------------
      if(ISTEP.eq.0) then
         call tstpr_init
      else
         call tstpr_solve
      endif

      return
      end
!***********************************************************************
!     initialise time stepper
      subroutine tstpr_init
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO, NDIM, NPERT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! NSTEPS
      include 'INPUT_DEF'
      include 'INPUT'           ! IFPERT, IFBASE, IFTRAN, PARAM, 
                                ! IFHEAT, IFUZAWA
      include 'MASS_DEF'
      include 'MASS'            ! BM1
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, TP, PRP
      include 'ADJOINT_DEF'
      include 'ADJOINT'         ! IFADJ
      include 'TIME_STEPPERD'
      include 'FSI'             !
      include 'NLFSI'           !

!     local variables
      integer i, ifield_tmp

!     functions
      real dnekclock, cht_glsc2_wt, op_glsc2_wt
      real fsi_glsc2_wt, nlfsi_glsc2_wt
!-----------------------------------------------------------------------
!     set parameters
      IFTST = .TRUE.
      if (TSTMODE.eq.0) IFTST = .FALSE.

      if (IFTST) then
!        timing
         TSTIME1=dnekclock()

!        initialise related packages
!        conjuagate heat transfer init
         if (IFHEAT) call cht_init

!        check nek5000 parameters
         if (.not.IFTRAN) then
            if (NIO.eq.0) write(6,*)
     $           'ERROR: time stepper requres transient equations'
            call exitt
         endif

         if (NSTEPS.eq.0) then
            if (NIO.eq.0) write(6,*)
     $           'ERROR: time stepper requires NSTEPS>0'
            call exitt
         endif

         if (PARAM(12).ge.0) then
            if (NIO.eq.0) write(6,*)
     $           'ERROR: time stepper assumes const. dt'
            call exitt
         endif

         if (.not.IFPERT) then
            if (NIO.eq.0) write(6,*)
     $         'ERROR: time stepper has to be run in perturbation mode'
            call exitt
         endif

         if (IFBASE) then
            if (NIO.eq.0) write(6,*)
     $           'ERROR: time stepper assumes constatnt base flow'
            call exitt
         endif

         if (NPERT.ne.1) then
            if (NIO.eq.0) write(6,*)
     $           'ERROR: time stepper requires NPERT=1'
            call exitt
         endif

!        make sure NSTEPS is bigger than the possible number of iterations
!        in time stepper phase
         NSTEPS = max(NSTEPS,TSTSTEP*TSTCMAX+10)

C        initialise cycle counters
         TSTISTEP = 0
         TSTVSTEP = 0

c        for timing
         TSTIMET=0.0
         TSTIMES=0.0
         TSTIMESA=0.0

c        vector length
         NVECAV  = NX1*NY1*NZ1*NELV ! velocity single component
         if(IFHEAT) then        !temperature
            NVECAT  = NX1*NY1*NZ1*NELT
         endif
         
!        We only keep the FSI variables on nid=0         
!        FSI/NLFSI         
         NVECAFSI=0
         NVECA_NLFSI=0
         if (IFFSI) then
            if (nid.eq.0) NVECAFSI = 1
         elseif (IFNLFSI) then
            if (nid.eq.0) NVECA_NLFSI = 1 
         endif

         NVECAP  = NX2*NY2*NZ2*NELV ! presure

!        print info
         if (NIO.eq.0) then
            write(6,*)
            write(6,*) 'TIME STEPPER initialised'
            if (TSTMODE.eq.1) then
               write(6,*) 'DIRECT mode'
            elseif (TSTMODE.eq.2) then
               write(6,*) 'ADJOINT mode'
            elseif (TSTMODE.eq.3) then
               write(6,*) 'Opt. init. cond.'
            endif
            write(6,*)
            write(6,*) 'Parameters:'

            write(6,'(A15,G13.5)') 'TOL = ',TSTTOL
            write(6,'(A15,I13)') 'Nstep = ',TSTSTEP
            write(6,'(A15,I13)') 'Ncmax = ',TSTCMAX
            write(6,*)
         endif

!        Only direct solver implemented for AMP based FSI
         if (IFFSI) then
           if (TSTMODE.eq.2) then
             if (nio.eq.0)
     $         write(6,*) 
     $           'Adjoint mode not implemented for Linear FSI for (AMP)'
             call exitt
           elseif (TSTMODE.eq.3) then
             if (nio.eq.0)
     $         write(6,*) 
     $           'Optimal init. cond. not implemented for Linear FSI'
             call exitt
           endif
         elseif (IFNLFSI) then
!          Optimal initial condition not yet implemented for NLFSI
           if (TSTMODE.eq.3) then
             if (nio.eq.0)
     $         write(6,*) 
     $     'Optimal init. cond. not implemented for Linear FSI (NLFSI)'
             call exitt
           endif

         endif 

!        place to initialise vector solver (arpack, power iteration)
         call stepper_vinit

!        setup simulation parameters
!        set time and iteration number
         TIME=0.0

!        should be the first step of every cycle performed with Uzawa 
!        turned on?
         IFUZAWA = TSTIFUZAWA

!        zero presure
         call rzero(PRP,NVECAP)

         IFADJ = .FALSE.
         if (TSTMODE.eq.2) then
!           Is it adjoint mode
            IFADJ = .TRUE.
         elseif  (TSTMODE.eq.3) then
!           If it is optimal initial condition save initial L2 norm
            if (IFFSI) then
              TSTL2INI = fsi_glsc2_wt()
            elseif (IFNLFSI) then
              TSTL2INI = nlfsi_glsc2_wt()
            elseif (ifheat) then 
!              TSTL2INI = cht_glsc2_wt(VXP,VYP,VZP,TP,
!     $                       VXP,VYP,VZP,TP,BM1)
            else
              TSTL2INI = op_glsc2_wt(VXP,VYP,VZP,
     $                       VXP,VYP,VZP,BM1)
            endif  

            if (TSTL2INI.eq.0.0) then
               if (NIO.eq.0)
     $              write(*,*) 'ERROR; tstpr_init, TSTL2INI = 0'
               call exitt
            endif
         endif

!        set cpfld for conjugated heat transfer
         if (IFHEAT) call cht_cpfld_set

!        timing
         TSTIME2=dnekclock()
         TSTIMET=TSTIME2-TSTIME1

      endif                     ! IFTST

!100  FORMAT('  SELECT(',I2,') = ',L)

      if (TSTMODE.eq.3.and.NIO.eq.0) write(6,200)
 200  FORMAT('TSTEPPER: opt. init. cond. direct phase start')

      return
      end
c***********************************************************************
c     subroutine to control time stepper after every nek5000 step
      subroutine tstpr_solve
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP, TIME
      include 'INPUT_DEF'
      include 'INPUT'           ! IFHEAT, IF3D, IFUZAWA
      include 'MASS_DEF'
      include 'MASS'            ! BM1
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, PRP, TP, VMULT, V?MASK
      include 'ADJOINT_DEF'
      include 'ADJOINT'         ! IFADJ
      include 'TIME_STEPPERD'
      include 'FSI'             ! eta,fsi_stiff
      include 'NLFSI'           ! psi,psi_stiff

!     global comunication in nekton
      integer NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
      common /nekmpi/ NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL

!     local variables
      real grw                  ! growth rate

!     functions
      real dnekclock, op_glsc2_wt, cht_glsc2_wt
      real fsi_glsc2_wt, nlfsi_glsc2_wt
!-----------------------------------------------------------------------
      if (IFTST) then

!        turn off Uzawa after first step
         IFUZAWA = .FALSE.

!        step counting
         TSTISTEP = TSTISTEP + 1

         if (mod(TSTISTEP,TSTSTEP).eq.0) then
!           stepper phase end
!           Generate new Krylov vector

!           check for the calculation mode
            if (TSTMODE.eq.3.and.(.not.IFADJ)) then
!              optimal initial condition
!              stamp output
               if (NIO.eq.0) write(6,200)
 200           FORMAT('TSTEPPER: opt. init. cond. adjoint phase start')

               IFADJ = .TRUE.

!              itaration count
               TSTISTEP = 0

!              should be the first step of every cycle performed with Uzawa 
!              turned on?
               IFUZAWA = TSTIFUZAWA

!              set time and iteration number
               TIME=0.0
               ISTEP=0

!              get L2 norm after direct phase
               if (IFFSI) then
                 TSTL2DIR = fsi_glsc2_wt()
               elseif (IFNLFSI) then
                 TSTL2DIR = nlfsi_glsc2_wt()
               elseif (ifheat) then 
!                 TSTL2DIR = cht_glsc2_wt(VXP,VYP,VZP,TP,
!     $                          VXP,VYP,VZP,TP,BM1)
               else
                 TSTL2DIR = op_glsc2_wt(VXP,VYP,VZP,
     $                          VXP,VYP,VZP,BM1)
               endif  

!              normalise vector
               grw = sqrt(TSTL2INI/TSTL2DIR)

               if (IFFSI) then
                 call fsi_opcmult(vxp,vyp,vzp,tp,eta,etav,grw)
               elseif (IFNLFSI) then
                 call nlfsi_opcmult(vxp,vyp,vzp,tp,psi,psiv,grw)
               elseif (ifheat) then
!                 call cht_opcmult(VXP,VYP,VZP,TP,grw)
               else 
                 call opcmult(VXP,VYP,VZP,grw)
               endif  

!              zero presure
               call rzero(PRP,NVECAP)

!              set cpfld for conjugated heat transfer
               if (IFHEAT) call cht_cpfld_set
            else
c              stepper phase counting
               TSTISTEP = 0
               TSTVSTEP = TSTVSTEP +1

!              timing
               TSTIME1=dnekclock()
!              average stepper phase length
               TSTIMES=TSTIME1-TSTIME2
               TSTIME2=TSTIME1
               TSTIMESA=(TSTIMESA*(TSTVSTEP-1)+TSTIMES)/real(TSTVSTEP)

c              stamp output
               if (NIO.eq.0) write(6,210) TSTVSTEP
 210           FORMAT('TSTEPPER: after',I5,' stepper phase')

               if (TSTMODE.eq.3) then
!                optimal initial condition
!                stamp output
                 if (NIO.eq.0) write(6,220)
 220             FORMAT('TSTEPPER: opt. init. cond. soln. rescaling')

!                get L2 norm after direct phase
                 if (iffsi) then
                   TSTL2DIR = fsi_glsc2_wt()
                 elseif (ifnlfsi) then
                   TSTL2DIR = nlfsi_glsc2_wt()
                 elseif (ifheat) then 
!                   TSTL2DIR = cht_glsc2_wt(VXP,VYP,VZP,TP,
!     $                            VXP,VYP,VZP,TP,BM1)
                 else
                   TSTL2DIR = op_glsc2_wt(VXP,VYP,VZP,
     $                            VXP,VYP,VZP,BM1)
                 endif  

!                normalise vector after whole cycle
                 grw = sqrt(TSTL2DIR/TSTL2INI)! add direct growth

                 if (iffsi) then
                   call fsi_opcmult(vxp,vyp,vzp,tp,eta,etav,grw)
                 elseif (ifnlfsi) then
                   call nlfsi_opcmult(vxp,vyp,vzp,tp,psi,psiv,grw)
                 elseif (ifheat) then
!                   call cht_opcmult (VXP,VYP,VZP,TP,grw)
                 else 
                   call opcmult (VXP,VYP,VZP,grw)
                 endif  

                 if (NIO.eq.0) then
                    write(*,*) 'Scaling factors:'
                    write(6,'(A15,G13.5)') 'TSTL2INI = ',TSTL2INI
                    write(6,'(A15,2G13.5)') 'TSTL2DIR = ',TSTL2DIR,
     $                   TSTL2DIR/TSTL2INI
                    write(6,'(A15,2G13.5)') 'TSTL2ADJ = ',TSTL2ADJ,
     $                   TSTL2ADJ/TSTL2INI
                 endif
               endif    ! tstmode.eq.3

!              place to run vector solver (arpack, power iteration)
               call stepper_vsolve

               if (LASTEP.eq.1) then
!                 finalise stepper
!                 timing
                  TSTIME3=dnekclock()
                  TSTIMET=TSTIMET+TSTIME3-TSTIME1

                  call tstpr_end
               else
!                 stepper restart;
!                 set time and iteration number
                  TIME=0.0
                  ISTEP=0

!                 should be the first step of every cycle performed with Uzawa 
!                 turned on?
                  IFUZAWA = TSTIFUZAWA

                  if (IFFSI) then
                    call fsi_arnoldi_restart
                  elseif (IFNLFSI) then
                    call nlfsi_arnoldi_restart
                  else
!                   zero presure
                    call rzero(PRP,NVECAP)
                  endif  

                  if (TSTMODE.eq.3) then
!                    optimal initial condition
!                    stamp output
                     if (NIO.eq.0) write(6,230)
 230            FORMAT('TSTEPPER: opt. init. cond. direct phase start')

                     IFADJ = .FALSE.

!                    get initial L2 norm
                     if (iffsi) then
                       TSTL2INI = fsi_glsc2_wt()
                     elseif (IFNLFSI) then
                       TSTL2INI = nlfsi_glsc2_wt()
                     elseif (ifheat) then 
!                       TSTL2INI = cht_glsc2_wt(VXP,VYP,VZP,TP,
!     $                                VXP,VYP,VZP,TP,BM1)
                     else
                       TSTL2INI = op_glsc2_wt(VXP,VYP,VZP,
     $                                VXP,VYP,VZP,BM1)
                     endif  

!                    set cpfld for conjugated heat transfer
                     if (IFHEAT) call cht_cpfld_set

                  endif
               endif

!              timing
               TSTIME3=dnekclock()
               TSTIMET=TSTIMET+TSTIME3-TSTIME1

            endif               ! TSTMODE.eq.3.and.(.not.IFADJ)
         endif                  ! mod(TSTISTEP,TSTSTEP).eq.0
      endif                     ! IFTST

      return
      end
c***********************************************************************
c     subroutine to finalise time stepper
      subroutine tstpr_end
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NID, NDIM, NPERT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP, NSTEPS, LASTEP
      include 'TIME_STEPPERD'
!-----------------------------------------------------------------------
      if (IFTST.and.(LASTEP.eq.1)) then
         if (NIO.eq.0) then
            write(6,*) ''
            write(6,*) 'TSTEPPER: finalize'
            write(6,*) '   Number of stepper cycles   ',TSTVSTEP
            write(6,*) '   Time spent in tstpr_solve  ',TSTIMET
            write(6,*) '   Average stepper phase time ',TSTIMESA
         endif
      endif

      return
      end
c***********************************************************************
