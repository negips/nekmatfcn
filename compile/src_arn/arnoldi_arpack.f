!=======================================================================
!     Adam Peplinski; 2015.10.09
!     Set of subroutines to solve eigenvalue problem with Arnoldi 
!     algorithm using PARPACK/ARPACK
!
!
!     Parameters used by this set of subroutines:
!     CHKPOINT:  
!     IFCHKPTRST - if restart
!
!     TIME_STEPPERD:
!     TSTMODE - 0 - no arnoldi
!               1 - direct mode
!               2 - adjoint mode
!               3 - initial optimal condition
!     TSTSTEP - frequency of calling arn_vsolve (number of time steps)
!     TSTCMAX - maximal number of arnoldi cycles
!     TSTTOL - PARPACK tolerance for Ritz values
!
!     ARNOLDI_ARPACKD:   
!     ARNKRYLOV - size of Krylov space
!     ARNEGV - number of calculated eigenvalues
!     ARNISTART - checkpoint file number for restart (reading)
!     ARNISTOP - checkpoint file number for saving data
!
!     For conjugated heat transfer
!     CONHT:
!     CHCST_SC, CHCFF_V, CHCFF_T - velocity and temperature scaling
!                                  factors
!     CHGRX, CHGRY, CHGRZ - gravitational acceleration
!     CHRA, CHPR - Rayleight and Prandtl numbers
!
!     Two possible modes are supported:
!     direct: A*x = lambda*x
!     inverse: A*x = lambda*M*x, M-mass matrix defining inner product
!     Inverse mode is preferred as it takes into account inner product
!     Simulation with temperature or passive scalars has to be performed
!     in inverse mode due to speciffic inner product.
!     To define ARPACK mode: direct or inverse
!#define ARPACK_DIRECT
#undef ARPACK_DIRECT
!
!     IMPORTANT!! No restart option for serial ARPACK version 
!     (only for PARPACK)
!=======================================================================
!***********************************************************************
!     read parameters ARPACK
      subroutine arp_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'ARNOLDI_ARPACKD'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr, len
      integer itmp(4)

!     namelists
      namelist /ARPACK/ ARNKRYLOV, ARNEGV, ARNISTART, ARNISTOP
!-----------------------------------------------------------------------
!     default values
      ARNKRYLOV = 100
      ARNEGV = 10
      ARNISTART = 1
      ARNISTOP = 2
!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=ARPACK,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading ARPACK parameters.$')

!     broadcast data
      if (NID.eq.0) then
         itmp(1) = ARNKRYLOV
         itmp(2) = ARNEGV
         itmp(3) = ARNISTART
         itmp(4) = ARNISTOP
      endif
      len = 4*ISIZE
      call bcast(itmp,len)
      if (NID.ne.0) then
         ARNKRYLOV = itmp(1)
         ARNEGV = itmp(2)
         ARNISTART = itmp(3)
         ARNISTOP = itmp(4)
      endif

      return
      end
!***********************************************************************
!     write parameters ARPACK
      subroutine arp_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'      
      include 'ARNOLDI_ARPACKD'         

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /ARPACK/ ARNKRYLOV, ARNEGV, ARNISTART, ARNISTOP
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=ARPACK,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing ARPACK parameters.$')

      return
      end
!***********************************************************************
!     initialise ARPACK (routine called by stepper)
      subroutine stepper_vinit
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO, NDIM, NPERT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! NSTEPS
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFHEAT
      include 'SOLN_DEF'
      include 'SOLN'            ! V?MASK, TMASK, V[XYZ]P, TP
      include 'CHKPOINT'        ! IFCHKPTRST
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'
      include 'FSI'
      include 'NLFSI'

!     ARPACK include file
      INCLUDE 'debug.h'

!     local variables
      integer i, ifield_tmp

      integer vecptr

!     functions
      real dnekclock
!-----------------------------------------------------------------------
!     set parameters
!     timing
      ARPTIME1=dnekclock()

!     check nek5000 parameters
!     Simulation with temperature or passive scalars has to be performed
!     in inverse mode due to speciffic inner product.
#ifdef ARPACK_DIRECT
!     standard eigenvalue problem A*x = lambda*x
      if(IFHEAT) then
         if (NIO.eq.0) write(6,*)
     $        'ERROR: IFHEAT requires #undef ARPACK_DIRECT'
         call exitt
      endif
#endif

      if (ARNKRYLOV.gt.LDIMA) then
         if (NIO.eq.0) write(6,*) 'ERROR: ARNKRYLOV too big'
         call exitt
      endif

      if (ARNEGV.ge.(ARNKRYLOV/2)) then
         if (NIO.eq.0) write(6,*)'ERROR: ARNEGV is >ARNKRYLOV/2'
         call exitt
      endif

c     make sure NSTEPS is bigger than the possible number of iteration in arnoldi
      NSTEPS = max(NSTEPS,TSTSTEP*ARNKRYLOV*TSTCMAX+10)

!     for timing
      ARPTIMET=0.0

!     related to restart
      NPARP = 0
      NCARP = 0
      RNMARP= 0.0

c     initialize ARPACK parameters
#ifdef ARPACK_DIRECT
c     standard eigenvalue problem A*x = lambda*x
      BMATARP='I'
#else
c     generalized eigenvalue problem A*x = lambda*B*x
      BMATARP='G'
#endif
c     eigenvalues of largest magnitude
      WHICHARP='LM'

      call izero(IPARP,11)
      call izero(IPNTARP,14)
c     exact shifts with respect to the current Hessenberg matrix
      IPARP(1)=1
c     maximum number of Arnoldi update iterations allowed
      IPARP(3)=TSTCMAX
#ifdef ARPACK_DIRECT
c     A*x = lambda*x
      IPARP(7)=1
#else
c     A*x = lambda*M*x, M symmetric positive definite; BMATARP='G'
      IPARP(7)=2
#endif
c     used size of WORKLA
      NWLARP = (3*ARNKRYLOV+6)*ARNKRYLOV

c     user supplied initial conditions
      INFARP=1

c     get eigenvectors
      RVARP=.true.
c     compute Ritz vectors
      HOWARP='A'
c     select should be specified for HOWARP='S'

c     no shift
      SIGARP(1) = 0.0
      SIGARP(2) = 0.0

c     vector length

!     single vector length in Krylov space
!     velocity
      NVECAS = NVECAV*NDIM

!     if we are adding pressure to the arnoldi vector      
      if (TSTIFPR) NVECAS = NVECAS + NVECAP

!     temperature
      if (IFHEAT) then
         NVECAS = NVECAS + NVECAT
      endif

!     FSI degreees of freedom only on nid=0 
!     FSI/NLFSI 
      if (IFFSI) then
        NVECAS = NVECAS + 2*NVECAFSI
      elseif (IFNLFSI) then
        NVECAS = NVECAS + 2*NVECA_NLFSI
      endif

      if (NVECAS.gt.LVAS) then
        if (NIO.eq.0) then
          write(6,*) 'ERROR: NVECAS too big'
          write(6,*) 'NVECAS = ', NVECAS
          write(6,*) 'LVAS   = ', LVAS
        endif
        call exitt
      else
        if (NIO.eq.0) then
          write(6,*) 'KRYLOV VECTOR LENGTH'
          write(6,*) 'NVECAS = ', NVECAS
          write(6,*) 'LVAS   = ', LVAS
        endif
      endif

c     initialize arrays
      call rzero(WORKDA,WDDIMA)
      call rzero(WORKLA,WLDIMA)
      call rzero(WORKEA,WEDIMA)
      call rzero(VBASEA,LVAS*LDIMA)
      call rzero(RESIDA,LVAS)
      call rzero(DRIARP,LDIMA*4)

c     info level from ARPACK
      ndigit = -3
      logfil = 6
      mngets = 0
      mnaitr = 2
      mnapps = 0
      mnaupd = 2
      mnaup2 = 2
      mneupd = 0

c     PLACE FOR RESTART
      if (IFCHKPTRST) then
c        read checkpoint
         call arn_rst_read
      else
c        if no restart fill RESIDA with initial conditions
c        V?MASK removes points at the wall and inflow
#ifdef ARPACK_DIRECT
c        A*x = lambda*x
!        velocity
         vecptr=1   
         call col3(RESIDA(vecptr),VXP,V1MASK,NVECAV)
         vecptr=vecptr+NVECAV
         call col3(RESIDA(vecptr),VYP,V2MASK,NVECAV)
         vecptr=vecptr+NVECAV
         if (IF3D) then
            call col3(RESIDA(vecptr),VZP,V3MASK,NVECAV)
            vecptr=vecptr+NVECAV
         endif 
!        no temperature here

#else
!        To the best of what I can understand, this is not really used.
c        x0
         if (IFFSI) then
           call fsi_flds2resida
         elseif (IFNLFSI) then
           call nlfsi_flds2resida
         else  
           call arn_flds2resida
         endif      

#endif

c        initialize rest of variables
c        first call
         IDOARP=0
      endif

!     ARPACK interface
      call arn_naupd

c     we should start stepper here
      if (IDOARP.ne.-1.and.IDOARP.ne.1) then
        if (NIO.eq.0) then
          write(6,*) 'ARNOLDI_ARPACK: '
          write(6,*) ' Error with arn_naupd, IDOARP = ', IDOARP
        endif
        call exitt
      endif
         
!     print info
      if (NIO.eq.0) then
         write(6,*)
         write(6,*) 'ARPACK initialised'
         write(6,*) 'Parameters:'
         write(6,'(A15,A1)') 'BMAT = ',BMATARP
         write(6,'(A15,A2)') 'WHICH = ',WHICHARP
         write(6,'(A15,G13.5)') 'TOL = ',TSTTOL
         write(6,'(A15,I13)') 'NEV = ',ARNEGV
         write(6,'(A15,I13)') 'NCV = ',ARNKRYLOV
         write(6,'(A15,I13)') 'IPARAM(1) = ',IPARP(1)
         write(6,'(A15,I13)') 'IPARAM(3) = ',IPARP(3)
         write(6,'(A15,I13)') 'IPARAM(7) = ',IPARP(7)
         write(6,'(A15,L)') 'RVEC = ',RVARP
         write(6,'(A15,A1)') 'HOWMNY = ',HOWARP
         if (HOWARP.eq.'S') then
            do i=1,LDIMA
               write(6,100) i,SELARP(i)
            enddo
         endif
      endif

c     timing
      ARPTIME2=dnekclock()
      ARPTIMET=ARPTIME2-ARPTIME1


 100  FORMAT('  SELECT(',I2,') = ',L)

      return
      end
c***********************************************************************
c     subroutine to create Krylov space, get Ritz values and restart
c     stepper phase (routine called by stepper)
      subroutine stepper_vsolve

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'INPUT_DEF'
      include 'INPUT'           ! IFHEAT, IF3D
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, TP, V?MASK
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !
      include 'FSI'                     ! eta
      include 'NLFSI'                   ! psi

!     functions
      real dnekclock

      integer vecptr
!-----------------------------------------------------------------------
!     timing
      ARPTIME1=dnekclock()

c     fill work array with velocity 
c     V?MASK removes points at the boundary
#ifdef ARPACK_DIRECT
c     A*x = lambda*x
!     velocity
      call col3(WORKDA(IPNTARP(2)),VXP,V1MASK,NVECAV)
      call col3(WORKDA(IPNTARP(2)+NVECAV),VYP,V2MASK,NVECAV)
      if (IF3D)
     $     call col3(WORKDA(IPNTARP(2)+2*NVECAV),VZP,V3MASK,NVECAV)
!     no temperature here

#else

c     y = A*x
      if (IFFSI) then
        call fsi_flds2workda
      elseif (IFNLFSI) then
        call nlfsi_flds2workda
      else 
        call arn_flds2workda
      endif  

#endif

!     ARPACK interface
      call arn_naupd

      if (IDOARP.eq.-2) then
!        checkpoint
         call arn_rst_save
         call arn_end
      elseif (IDOARP.eq.99) then
!        finalise
         call arn_esolve
         call arn_end
      elseif (IDOARP.eq.-1.or.IDOARP.eq.1) then
!        stepper restart, nothing to do
      else
         if (NIO.eq.0) then
            write(6,*) 'ARNOLDI_ARPACK: '
            write(6,*) ' Error with arn_naupd, IDOARP = ',
     $           IDOARP
         endif
         call exitt
      endif

c     timing
      ARPTIME2=dnekclock()
      ARPTIMET=ARPTIMET+ARPTIME2-ARPTIME1

      return
      end
c***********************************************************************
c     ARPACK postprocessing
      subroutine arn_esolve
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO, NID, LDIMT1
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP, DT, LASTEP
      include 'SOLN_DEF'
      include 'SOLN'            ! VX, VY, VZ, VMULT, V?MASK
      include 'INPUT_DEF'
      include 'INPUT'           ! IFXYO,IFPO,IFVO,IFTO,IFPSO,IF3D,IFHEAT
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !
      include 'FSI'
      include 'NLFSI'

!     local variables
      integer i,n, iunit, ierror
      real dumm
      logical lifxyo, lifpo, lifvo, lifto, lifpso(LDIMT1)
!     functions
      real dnekclock

c     global comunication in nekton
      integer NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
      common /nekmpi/ NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL

      integer vecptr
!-----------------------------------------------------------------------
      if (IDOARP.eq.99) then

c        timing
         ARPTIME1=dnekclock()

         if (NIO.eq.0) write(6,405) IPARP(5)
 405     FORMAT('ARNOLDI_ARPACK: ',I4,' eigenvectors converged;
     $postprocessing')

#ifdef MPI
         call pdneupd(NEKCOMM,RVARP,HOWARP,SELARP,DRIARP,DRIARP(1,2),
     $        VBASEA,LVAS,SIGARP(1),SIGARP(2),WORKEA,BMATARP,NVECAS,
     $        WHICHARP,ARNEGV,TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,
     $        IPARP,IPNTARP,WORKDA,WORKLA,NWLARP,IERRARP)
#else
         call dneupd(RVARP,HOWARP,SELARP,DRIARP,DRIARP(1,2),
     $        VBASEA,LVAS,SIGARP(1),SIGARP(2),WORKEA,BMATARP,NVECAS,
     $        WHICHARP,ARNEGV,TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,
     $        IPARP,IPNTARP,WORKDA,WORKLA,NWLARP,IERRARP)
#endif

c        timing
         ARPTIME2=dnekclock()
         ARPTIMET=ARPTIMET+ARPTIME2-ARPTIME1

         if (IERRARP.eq.0) then

            if (NIO.eq.0) then
               write(6,*) 'ARNOLDI_ARPACK:'
               write(6,*) 'writing eigenvalues and eigenvectors.'
            endif

            ierror=0
!           open file 
            if (NID.eq.0) then
!              find free unit
               call IO_freeid(iunit, ierror)
               if (ierror.eq.0) then
                  open (unit=iunit,file='eigenvalues.txt',
     $                 action='write', iostat=ierror)
                  write(unit=iunit,FMT=410,iostat=ierror)
     $               'I','re(RITZ)','im(RITZ)','(ln|RITZ|)/DT_KRYLOV',
     $               '(arg(RITZ))/DT_KRYLOV'
! 410              FORMAT(9x,'I',17x,'re(RITZ)',17x,'im(RITZ)',17x,
!     $                 'ln|RITZ|',16x,'arg(RITZ)')
 410              FORMAT(5x,A5,4(1x,A24))

               endif
            endif
c           error check
            call err_chk(ierror,'Error opening eigenv. file.$')

!           integration time
            dumm = DT*TSTSTEP
            dumm = 1.0/dumm

!           copy and set output parameters
            lifxyo = IFXYO
            IFXYO = .TRUE.
            lifvo= IFVO
            IFVO = .true.
            lifpo= IFPO
            if (TSTIFPR) then
              IFPO = .true.
            else            
              IFPO = .false.
            endif
            lifto= IFTO
            if (IFHEAT) then
               IFTO = .TRUE.
            else
               IFTO = .FALSE.
            endif
            do i=1,LDIMT1
               lifpso(i)= IFPSO(i)
               IFPSO(i) = .false.
            enddo

c     we have to take into account storrage of imaginary and real
c     parts of eigenvectors in arpack.
c     The complex Ritz vector associated with the Ritz value 
c     with positive imaginary part is stored in two consecutive 
c     columns.  The first column holds the real part of the Ritz 
c     vector and the second column holds the imaginary part.  The 
c     Ritz vector associated with the Ritz value with negative 
c     imaginary part is simply the complex conjugate of the Ritz
c     vector associated with the positive imaginary part.

            ierror=0
            do i=1,IPARP(5)
!              copy eigenvectors to perturbation variables

!              Velocities            
               vecptr=1   
               call copy(VXP,VBASEA(vecptr,i),NVECAV)
               vecptr=vecptr+NVECAV
               call copy(VYP,VBASEA(vecptr,i),NVECAV)
               vecptr=vecptr+NVECAV
               if (IF3D) then
                  call copy(VZP,VBASEA(vecptr,i),NVECAV)
                  vecptr=vecptr+NVECAV
               endif
 
!              Pressure. (if we are using a Semi-norm)
               if (TSTIFPR) then
                 call copy(prp,VBASEA(vecptr,i),NVECAP)
                 vecptr=vecptr+NVECAP
               endif  

!              Temperature   
               if(IFHEAT) then
                 call copy(TP,VBASEA(vecptr,i),NVECAT)
                 vecptr=vecptr+NVECAT
               endif

               if (IFFSI) then
                  if (nid.eq.0) then                     
!                   Position                     
                    call copy(eta,VBASEA(vecptr,i),NVECAFSI)
                    vecptr=vecptr+NVECAFSI
!                   Velocity                  
                    call copy(etav,VBASEA(vecptr,i),NVECAFSI)
                    vecptr=vecptr+NVECAFSI
                  endif  
               elseif (IFNLFSI) then
                  if (nid.eq.0) then                     
!                   Position                     
                    call copy(psi,VBASEA(vecptr,i),NVECA_NLFSI)
                    vecptr=vecptr+NVECA_NLFSI
!                   Velocity                  
                    call copy(psiv,VBASEA(vecptr,i),NVECA_NLFSI)
                    vecptr=vecptr+NVECA_NLFSI
                  endif  
               endif   


!              get growth rate; get eigenvalues of continuous operator
               DRIARP(i,3) = log(sqrt(DRIARP(i,1)**2+
     $              DRIARP(i,2)**2))*dumm
               DRIARP(i,4) = atan2(DRIARP(i,2),DRIARP(i,1))*dumm

               if (NID.eq.0)  write(unit=iunit,FMT=411,iostat=ierror)
     $              i,DRIARP(i,1),DRIARP(i,2),DRIARP(i,3),DRIARP(i,4)

 411           FORMAT(5x,I5,4(1x,E24.15E3))

!              Output eigenvectors
               istep=i            ! putting cyle no as identifier for eigen vector no
               time =DRIARP(i,4)  ! time as the identifier for angular frequency
               if(IFHEAT) then
                  call outpost2(VXP,VYP,VZP,PRP,TP,1,'egv')
               else
                  call outpost2(VXP,VYP,VZP,PRP,TP,0,'egv')
               endif

!              Outpost structural part of the eigenvector               
               if (iffsi) then
                 call fsi_outpostegv(i)
               elseif (ifnlfsi) then 
                 call nlfsi_outpostegv(i)
               endif   

               istep=0
               time=0.

!              possible place to test error
            enddo
!           error check
            call err_chk(ierror,'Error witing to eigenv. file.$')

!           put output variables back
            IFXYO = lifxyo
            IFPO = lifpo
            IFVO = lifvo
            IFTO = lifto
            do i=1,LDIMT1
               IFPSO(i) = lifpso(i)
            enddo

!           close eigenvalue file
            if (NID.eq.0)  close(unit=iunit)

         else                   ! IERRARP
            if (NIO.eq.0) then
               write(6,*) 'ARNOLDI_ARPACK:'
               write(6,*) ' Error with _neupd, info = ', IERRARP
               write(6,*) ' Check the documentation of _naupd.'
            endif
            call exitt
         endif                  ! IERRARP

c        finish run
         LASTEP=1

      endif

      return
      end
c***********************************************************************
c     interface to pdnaupd
      subroutine arn_naupd
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO, NDIM, N[XYZ]1
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFHEAT
      include 'SOLN_DEF'
      include 'SOLN'            ! V?MASK, TMASK, V[XYZ]P, TP, PRP
      include 'MASS_DEF'
      include 'MASS'            ! BM1
      include 'TSTEP_DEF'
      include 'TSTEP'           ! IFIELD
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'FSI'                     ! eta, fsi_stiff
      include 'NLFSI'                   ! psi, nlfsi_stiff

c     global comunication in nekton
      integer NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
      common /nekmpi/ NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
!     local variables
      integer ifield_tmp

      integer vecptr
!-----------------------------------------------------------------------

#ifdef MPI
      call pdnaupd(NEKCOMM,IDOARP,BMATARP,NVECAS,WHICHARP,ARNEGV,
     $     TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,IPARP,IPNTARP,WORKDA,
     $     WORKLA,NWLARP,INFARP,NPARP,RNMARP,NCARP)
#else
      call dnaupd(IDOARP,BMATARP,NVECAS,WHICHARP,ARNEGV,
     $     TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,IPARP,IPNTARP,WORKDA,
     $     WORKLA,NWLARP,INFARP)
#endif

c     error check
      if (INFARP.lt.0) then
         if (NIO.eq.0) then
            write(6,*) 'ARNOLDI_ARPACK: '
            write(6,*) ' Error with _naupd, info = ', INFARP
            write(6,*) ' Check the documentation of _naupd.'
         endif
         call exitt
      endif

      if (IDOARP.eq.2) then
         do
c           A*x = lambda*M*x
c           multiply by weights and masks

            if (IFFSI) then
!             FSI weights
              call fsi_weight_fcn
            elseif (IFNLFSI) then
!             NLFSI weights
              call nlfsi_weight_fcn
            else
              call arn_weight_fcn
            endif

#ifdef MPI
            call pdnaupd(NEKCOMM,IDOARP,BMATARP,NVECAS,WHICHARP,
     $           ARNEGV,TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,IPARP,
     $           IPNTARP,WORKDA,WORKLA,NWLARP,INFARP,NPARP,RNMARP,
     $           NCARP)
#else
            call dnaupd(IDOARP,BMATARP,NVECAS,WHICHARP,ARNEGV,
     $           TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,IPARP,IPNTARP,
     $           WORKDA,WORKLA,NWLARP,INFARP)
#endif

!           error check
            if (INFARP.lt.0) then
               if (NIO.eq.0) then
                  write(6,*) 'ARNOLDI_ARPACK: '
                  write(6,*) ' Error with _naupd, info = ', INFARP
                  write(6,*) ' Check the documentation of _naupd.'
               endif
               call exitt
            endif
            if (IDOARP.ne.2) then
              exit
            endif
         enddo
      endif                     ! IDOARP.eq.2

!     restart stepper
      if (IDOARP.eq.-1.or.IDOARP.eq.1) then

         if (NIO.eq.0) write(6,*) 'ARNOLDI_ARPACK: restarting stepper'

c        move renormed data back to nek arrays
         if (IFFSI) then
           call fsi_workda2flds
         elseif(IFNLFSI) then
           call nlfsi_workda2flds
         else
           call arn_workda2flds
         endif  

c        make sure the velocity and temperature fields are continuous at 
c        element faces and edges
         ifield_tmp = IFIELD
         IFIELD = 1
         call opdssum(VXP,VYP,VZP)
         call opcolv (VXP,VYP,VZP,VMULT)
         if(IFHEAT) then
            IFIELD = 2
            CALL dssum(TP,NX1,NY1,NZ1)
            call col2 (TP,TMULT,NVECAT)
         endif
         IFIELD = ifield_tmp
      endif                     ! IDOARP.eq.-1.or.IDOARP.eq.1

      return
      end subroutine arn_naupd
c***********************************************************************
c     subroutine to take care of savings for restart 
      subroutine arn_rst_save
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'TSTEP_DEF'
      include 'TSTEP'           ! LASTEP
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD' !
      include 'FSI'             ! save eta,etav
      include 'NLFSI'           ! save psi,psiv
     
!-----------------------------------------------------------------------
c     save checkpoint for IDOARP=-2
      if (IDOARP.eq.-2) then

         if (NIO.eq.0) then
            write(6,*) 'ARNOLDI_ARPACK: '
            write(6,*) 'IDOARP  = ', IDOARP
            write(6,*) 'Writing checkpoint'
         endif

!        save parameters and WORKLA; independent on processor; 
!        serial output
         call arn_write_par('ARP')

c        save big arrays; parallel output
         call mfo_arnv('ARV')

!        FSI
         if (IFFSI) then 
           call fsi_arnoldi_rstsave
         elseif  (IFNLFSI) then
           call nlfsi_arnoldi_rstsave
         endif  

c        this is the last step
         LASTEP=1

      endif

      return
      end subroutine arn_rst_save
c***********************************************************************
c     subroutine to read from checkpoints
      subroutine arn_rst_read
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !
      include 'FSI'                     ! IFFSI
      include 'NLFSI'                   ! IFNLFSI
!-----------------------------------------------------------------------
      if (NIO.eq.0) then
         write(6,*) 'ARNOLDI_ARPACK: '
         write(6,*) 'Reading checkpoint'
      endif
c     read parameters and WORKLA; independent on processor; serial input
      call arn_read_par('ARP')
      
c     read big arrays; parallel input
      call mfi_arnv('ARV')

!     FSI
      if (IFFSI) then 
        call fsi_arnoldi_rstread
      elseif  (IFNLFSI) then
        call nlfsi_arnoldi_rstread
      endif 

      return
      end subroutine arn_rst_read
c***********************************************************************
c     subroutine to finalise arnoldi
      subroutine arn_end
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NID, NDIM, NPERT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP, NSTEPS, LASTEP
      include 'ARNOLDI_ARPACKD'
!-----------------------------------------------------------------------
      if ((LASTEP.eq.1)) then
         if (NIO.eq.0) then
            write(6,*) ''
            write(6,*) 'ARPACK: finalize'
            write(6,*) '   Time spent in ARPACK      ',ARPTIMET
         endif
      endif

      return
      end subroutine arn_end
c-----------------------------------------------------------------------
c     Subroutines to manipulate files for arnoldi restart
c-----------------------------------------------------------------------
c     Subroutines to write/read procesor independent variables;
c     big arrays are saved/read by mfo_out_arnv/
c***********************************************************************
      subroutine arn_write_par(prefix)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'

!     argument list
      character*3 prefix

!     local variables
      character(LEN=132) fname

      character(LEN=6)  str

      integer lwdsizo, k, ierr

      logical lifreguo
!-----------------------------------------------------------------------
      call nekgsync()

c     copy and set output parameters
      lwdsizo= WDSIZO
      WDSIZO = 8

      lifreguo= IFREGUO
      IFREGUO = .false.

c     this is done by master node only
      if (NID.eq.0) then
c     create file name
         call IO_mfo_fname(prefix,fname,SESSION,k)
         write(str,540) ARNISTOP
 540     format(i5.5)
         fname(k:k+4)=str(1:5)
        
c     open file; only serial
         call IO_mbyte_open_srl(fname,NID,ierr)

         call mfo_arnp

c     close the file; only serial
         call byte_close(ierr)
      endif

c     put output variables back
      WDSIZO = lwdsizo

      IFREGUO = lifreguo

      return
      end
c***********************************************************************
      subroutine mfo_arnp
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'

!     local variables
      character*16 hdr
      integer ahdsize
      parameter (ahdsize=16)

      integer ibsw_out, i, itmp(33), ierr

      real*4 test_pattern

      real*4 rtmp4(6), workla4(2*WLDIMA)
      real*8 rtmp8(3), workla8(WLDIMA)
      equivalence (rtmp4,rtmp8)
      equivalence (workla4,workla8)
!-----------------------------------------------------------------------
c     write IDOARP and character varialbes 
      call blank(hdr,ahdsize)

      write(hdr,1) IDOARP,BMATARP,WHICHARP,TSTMODE! 14
 1    format('#arp',1x,i2,1x,a1,1x,a2,1x,i1)

      ! if we want to switch the bytes for output
      ! switch it again because the hdr is in ASCII
      call get_bytesw_write(ibsw_out)
      if (ibsw_out.ne.0) call set_bytesw_write(0)  

      call byte_write(hdr,ahdsize/4,ierr)

c     write test pattern for byte swap
      test_pattern = 6.54321

      call byte_write(test_pattern,1,ierr)

c     collect and write integer varialbes
      itmp(1) = NVECAS
      itmp(2) = ARNEGV
      itmp(3) = ARNKRYLOV
      itmp(4) = NWLARP
      itmp(5) = INFARP
      itmp(6) = NPARP
      itmp(7) = NCARP
      itmp(8) = TSTSTEP
      do i=1,11
         itmp(8+i) = IPARP(i)
      enddo
      do i=1,14
         itmp(19+i) = IPNTARP(i)
      enddo

      call byte_write(itmp,33,ierr)

c     collect and write real variables
      rtmp8(1) = TSTTOL
      rtmp8(2) = RNMARP
      rtmp8(3) = DT

      call byte_write(rtmp4,6,ierr)

c     write WORKLA
      call copy(workla8,WORKLA,NWLARP)

      call byte_write(workla4,2*NWLARP,ierr)

      return
      end
c***********************************************************************
      subroutine arn_read_par(prefix)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'

!     argument list
      character*3 prefix

!     local variables
      character(LEN=132) fname

      character(LEN=6)  str

      integer lwdsizo, i, k, ierr, error

      logical lifreguo
!-----------------------------------------------------------------------
      call nekgsync()

c     copy and set output parameters
      lwdsizo= WDSIZO
      WDSIZO = 8

      lifreguo= IFREGUO
      IFREGUO = .false.

      error=0

c     this is done by master node only
      if (NID.eq.0) then
c     create file name
         call IO_mfo_fname(prefix,fname,SESSION,k)
         write(str,540) ARNISTART
 540     format(i5.5)
         fname(k:k+4)=str(1:5)
        
c     open file; only serial
         call IO_mbyte_open_srl(fname,NID, ierr)

c     read parameters
         call mfi_arnp

c     close the file; only serial
         call byte_close(ierr)

c     check and copy parameters
c     is it correct restart
         if (IDOARP0.ne.-2) then
            write(6,*) 'ERROR: arn_read_par, wrong IDOARP0; abort.'
            write(6,*) 'IDOARP0 = ', IDOARP0
            error=1
         endif

c     is it the same ARPACK mode
         if (BMATARP0.ne.BMATARP) then
            write(6,*) 'ERROR: arn_read_par, different ARPACK modes.'
            write(6,*) 'BMATARP0 = ', BMATARP0
            write(6,*) 'BMATARP  = ', BMATARP
            error=1
         endif

c     do we look for the same eigenvectors
         if (WHICHARP0.ne.WHICHARP) then
            write(6,*) 'ERROR: arn_read_par, different mode selection.'
            write(6,*) 'WHICHARP0 = ', WHICHARP0
            write(6,*) 'WHICHARP  = ', WHICHARP
            error=1
         endif

c     is it the same integration mode
         if (TSTMODE0.ne.TSTMODE) then
            write(6,*) 'ERROR: arn_read_par, wrong simulation mode.'
            write(6,*) 'TSTMODE0 = ', TSTMODE0
            write(6,*) 'TSTMODE  = ', TSTMODE
            error=1
         endif

!     this should be removed later as it does not allow to change processor number
c     is the length of the vector the same
         if (NVECAS0.ne.NVECAS) then
            write(6,*) 'ERROR: arn_read_par, different vector lenght;
     $abort. Check IFHEAT!'
            write(6,*) 'NVECAS0 = ', NVECAS0
            write(6,*) 'NVECAS  = ', NVECAS
            error=1
         endif

c     what is the size of Krylov space
c     would it be possible to change this?; related NPARP, NWLARP,
c     IPNTARP
         if (ARNKRYLOV0.ne.ARNKRYLOV) then
            write(6,*) 'ERROR: arn_read_par, different Krylov space
     $size.'
            write(6,*) 'ARNKRYLOV0 = ', ARNKRYLOV0
            write(6,*) 'ARNKRYLOV  = ', ARNKRYLOV
            error=1
         endif

         if (NWLARP0.ne.NWLARP) then
            write(6,*) 'ERROR: arn_read_par, different size of work
     $array'
            write(6,*) 'NWLARP0 = ', NWLARP0
            write(6,*) 'NWLARP  = ', NWLARP
            error=1
         endif

c     stopping criterion 
         if (TSTTOL0.ne.TSTTOL) then
            write(6,*) 'WARNING: arn_read_par, different stopping
     $criterion'
            write(6,*) 'TSTTOL0 = ', TSTTOL0
            write(6,*) 'TSTTOL  = ', TSTTOL
         endif

c     number of eigenvalues
         if (ARNEGV0.ne.ARNEGV) then
            write(6,*) 'WARNING: arn_read_par, different number of 
     $eigenvalues'
            write(6,*) 'ARNEGV0 = ', ARNEGV0
            write(6,*) 'ARNEGV  = ', ARNEGV
         endif

c     stepper phase length
         if (DTARP0.ne.DT) then
            write(6,*) 'WARNING: arn_read_par, different timestep'
            write(6,*) 'DTARP0 = ', DTARP0
            write(6,*) 'DT     = ', DT
         endif

         if (TSTSTEP0.ne.TSTSTEP) then
            write(6,*) 'WARNING: arn_read_par, different number of 
     $steps in stepper phase'
            write(6,*) 'TSTSTEP0 = ', TSTSTEP0
            write(6,*) 'TSTSTEP  = ', TSTSTEP
         endif

c     check IPARP
         if (IPARP0(1).ne.IPARP(1)) then
            write(6,*) 'ERROR: arn_read_par, different shift in ARPACK'
            write(6,*) 'IPARP0(1) = ', IPARP0(1)
            write(6,*) 'IPARP(1)  = ', IPARP(1)
            error=1
         endif

         if (IPARP0(3).ne.IPARP(3)) then
            write(6,*) 'WARNING: arn_read_par, different cycle number'
            write(6,*) 'IPARP0(3) = ', IPARP0(3)
            write(6,*) 'IPARP(3)  = ', IPARP(3)
         endif

         if (IPARP0(7).ne.IPARP(7)) then
            write(6,*) 'ERROR: arn_read_par, different ARPACK modes'
            write(6,*) 'IPARP0(7) = ', IPARP0(7)
            write(6,*) 'IPARP(7)  = ', IPARP(7)
            error=1
         endif

c        copy rest of parameters
         NPARP = NPARP0
         NCARP = NCARP0
         INFARP= INFARP0
         RNMARP= RNMARP0
         do i=4,11
            IPARP(i) = IPARP0(i)
         enddo
         IPARP(2) = IPARP0(2)
         do i=1,14
            IPNTARP(i) = IPNTARP0(i)
         enddo
      endif                     ! NID

      call err_chk(error,'Error restarting arnoldi.$')

c      call nekgsync()


      IDOARP = -2
c     broadcast
      call bcast(NPARP,ISIZE)
      call bcast(NCARP,ISIZE)
      call bcast(INFARP,ISIZE)
      call bcast(IPARP,11*ISIZE)
      call bcast(IPNTARP,14*ISIZE)
      call bcast(RNMARP,WDSIZE)

      call bcast(WORKLA,NWLARP*WDSIZE)


c     put output variables back
      WDSIZO = lwdsizo

      IFREGUO = lifreguo

      return
      end
c***********************************************************************
      subroutine mfi_arnp
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'TIME_STEPPERD'
      INCLUDE 'ARNOLDI_ARPACKD'
!     local variables
      character*16 hdr
      character*4 dummy
      integer ahdsize
      parameter (ahdsize=16)

      integer ibsw_out, i, itmp(33), ierr

      real*4 test_pattern

      real*4 rtmp4(6), workla4(2*WLDIMA)
      real*8 rtmp8(3), workla8(WLDIMA)
      equivalence (rtmp4,rtmp8)
      equivalence (workla4,workla8)

      logical if_byte_swap_test, if_byte_sw_loc
!     functions
      integer indx2
!-----------------------------------------------------------------------
c     read IDOARP and character varialbes 
      call blank(hdr,ahdsize)

      call byte_read(hdr,ahdsize/4,ierr)

      if (indx2(hdr,132,'#arp',4).eq.1) then
         read(hdr,*) dummy,IDOARP0,BMATARP0,WHICHARP0,TSTMODE0! 14
      else
         write(6,*) 'ERROR: mfi_arnp, wrong header; abort.'
         call exitt
      endif

c     read test pattern for byte swap
      call byte_read(test_pattern,1,ierr)
c     determine endianess
      if_byte_sw_loc = if_byte_swap_test(test_pattern,ierr) 

c     read integer varialbes
      call byte_read(itmp,33,ierr)
      if (if_byte_sw) call byte_reverse(itmp,33,ierr)

      NVECAS0 = itmp(1)
      ARNEGV0  = itmp(2)
      ARNKRYLOV0  = itmp(3)
      NWLARP0 = itmp(4)
      INFARP0 = itmp(5)
      NPARP0  = itmp(6)
      NCARP0  = itmp(7)
      TSTSTEP0 = itmp(8)
      do i=1,11
         IPARP0(i) = itmp(8+i)
      enddo
      do i=1,14
         IPNTARP0(i) = itmp(19+i)
      enddo

c     read real variables
      call byte_read(rtmp4,6,ierr)
      if (if_byte_sw) call byte_reverse(rtmp4,6,ierr)

      TSTTOL0 = rtmp8(1)
      RNMARP0 = rtmp8(2)
      DTARP0  = rtmp8(3)

c     read WORKLA
      if (NWLARP0.le.WLDIMA) then
         call byte_read(workla4,2*NWLARP0,ierr)
         if (if_byte_sw) call byte_reverse(workla4,2*NWLARP0,ierr)

         call copy(WORKLA,workla8,NWLARP0)
      else
         write(6,*) 'ARNOLDI_ARPACKD: restart error'
         write(6,*) 'NWLARP0 = ',NWLARP0
         write(6,*) 'WLDIMA  = ',WLDIMA
         call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
c     following subroutines are modiffications of 
c     
c     mfo_outfld
c     from prepost.f;
c     mfi
c     from ic.f
c***********************************************************************
      subroutine mfo_arnv(prefix)  ! muti-file output
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'

!     argument list
      character*3 prefix

!     local variables
      integer*8 offs0,offs,nbyte,stride,strideB,nxyzo8
      integer lwdsizo, i, k, ierr
      integer ioflds, nout
      real dnbyte, tio, tiostart

      character(LEN=132)  fname

      character(LEN=6)  str

      logical lifxyo, lifpo, lifvo, lifto, lifreguo, lifpso(LDIMT1)

!     functions
      real dnekclock_sync, glsum

!     scratch space
      real UR1(LXO*LXO*LXO*LELT), UR2(LXO*LXO*LXO*LELT),
     $     UR3(LXO*LXO*LXO*LELT)
      common /SCRUZ/  UR1, UR2, UR3
!-----------------------------------------------------------------------
      tiostart=dnekclock_sync()

c     copy and set output parameters
      lwdsizo= WDSIZO
      WDSIZO = 8

      lifreguo= IFREGUO
      IFREGUO = .false.
      lifxyo= IFXYO
      IFXYO = .false.
      lifpo= IFPO
      IFPO = .false.
      lifvo= IFVO
      IFVO = .true.
      lifto= IFTO
      IFTO = .false.
      do i=1,LDIMT1
         lifpso(i)= IFPSO(i)
         IFPSO(i) = .false.
      enddo

      nout = NELT
      NXO  = NX1
      NYO  = NY1
      NZO  = NZ1

c     set offset
      offs0 = iHeaderSize + 4 + isize*nelgt

      if (nid.eq.pid0) then         ! open files on i/o nodes
!        call mfo_open_files(prefix)
c        create file name
         call IO_mfo_fname(prefix,fname,SESSION,k)
         write(str,5) ARNISTOP
 5       format(i5.5)
         fname(k:k+4)=str(1:5)
         call mbyte_open(fname,fid0,ierr)
c         write(6,*) nid,fid0,' FILE:',fname
      endif

      call mfo_write_hdr                     ! create element mapping +
                                             ! write hdr
      nxyzo8  = NXO*NYO*NZO
      strideB = nelB * nxyzo8*WDSIZO
      stride  = nelgt* nxyzo8*WDSIZO

      ioflds = 0
      ! dump all fields based on the t-mesh to avoid different
      ! topologies in the post-processor

c     resid array
      call  mfo_singlev(ioflds,nout,offs0,stride,strideB,
     $     UR1,UR2,UR3,RESIDA)

c     workd array
      do i=0,2
         call  mfo_singlev(ioflds,nout,offs0,stride,strideB,
     $     UR1,UR2,UR3,WORKDA(1+NVECAS*i))
      enddo

c     krylov space
      do i=1,ARNKRYLOV
         call  mfo_singlev(ioflds,nout,offs0,stride,strideB,
     $     UR1,UR2,UR3,VBASEA(1,i))
      enddo

      dnbyte = 1.*ioflds*nout*WDSIZO*NXO*NYO*NZO

c     put output variables back
      WDSIZO = lwdsizo

      IFREGUO = lifreguo
      IFXYO = lifxyo
      IFPO = lifpo
      IFVO = lifvo
      IFTO = lifto
      do i=1,LDIMT1
         IFPSO(i) = lifpso(i)
      enddo

      if (NID.eq.PID0) 
#ifdef MPIIO
     &   call byte_close_mpi(ifh_mbyte,ierr)
#else
     &   call byte_close(ierr)
#endif

      tio = dnekclock_sync()-tiostart
      dnbyte = glsum(dnbyte,1)
      dnbyte = dnbyte + iHeaderSize + 4. +ISIZE*NELGT
      dnbyte = dnbyte/1024/1024
      if(nid.eq.0) write(6,7)  ISTEP,TIME,dnbyte,dnbyte/tio,
     &             NFILEO
    7 format(/,i9,1pe12.4,' done :: Write checkpoint',/,
     &       30X,'file size = ',3pG12.2,'MB',/,
     &       30X,'avg data-throughput = ',0pf7.1,'MB/s',/,
     &       30X,'io-nodes = ',i5,/)

      return
      end
c-----------------------------------------------------------------------
      subroutine mfi_arnv(prefix)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TIME_STEPPERD'
      INCLUDE 'ARNOLDI_ARPACKD'

!     argument list
      character*3 prefix

!     local variables
      character(LEN=132)  fname

      character(LEN=6)  str

      character*132 hdr

      integer e, i, k, iofldsr, ierr, iofldr
      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8
      real dnbyte, tio, tiostart

!     functions
      real dnekclock, glsum

!     scratch space
      integer lwk
      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      real wk(lwk)
      common /scrns/ wk

      real UR1(LX1,LY1,LZ1,LELT), UR2(LX1,LY1,LZ1,LELT),
     $     UR3 (LX1,LY1,LZ1,LELT)
      COMMON /scruz/ UR1, UR2, UR3
c-----------------------------------------------------------------------
      tiostart=dnekclock()

c     create file name
      call IO_mfo_fname(prefix,fname,SESSION,k)
      write(str,5) ARNISTART
 5    format(i5.5)
      fname(k:k+4)=str(1:5)

      call mfi_prepare(fname)       ! determine reader nodes +
                                    ! read hdr + element mapping 

      offs0   = iHeadersize + 4 + ISIZE*NELGR
      nxyzr8  = NXR*NYR*NZR
      strideB = nelBr* nxyzr8*WDSIZR
      stride  = nelgr* nxyzr8*WDSIZR

c     read arrays
      iofldsr = 0
c     resid array
      call mfi_singlev(iofldr,offs0,stride,strideB,
     $     UR1,UR2,UR3,RESIDA(1))

c     workd array
      do i=0,2
         call mfi_singlev(iofldr,offs0,stride,strideB,
     $        UR1,UR2,UR3,WORKDA(1+NVECAS*i))
      enddo

c     krylov space
      do i=1,ARNKRYLOV
         call mfi_singlev(iofldr,offs0,stride,strideB,
     $        UR1,UR2,UR3,VBASEA(1,i))
      enddo

      nbyte = 0
      if(nid.eq.pid0r) nbyte = iofldsr*nelr*wdsizr*nxr*nyr*nzr

c     close files
#ifdef MPIIO
      if (nid.eq.pid0r) call byte_close_mpi(ifh_mbyte,ierr)
#else
      if (nid.eq.pid0r) call byte_close(ierr)
#endif
      call nekgsync
      tio = dnekclock()-tiostart

      dnbyte = nbyte
      nbyte = glsum(dnbyte,1)
      nbyte = nbyte + iHeaderSize + 4 + isize*nelgr

      if(nid.eq.0) write(6,7) istep,time,
     &             nbyte/tio/1024/1024/10,
     &             nfiler
    7 format(/,i9,1pe12.4,' done :: Read checkpoint data',/,
     &       30X,'avg data-throughput = ',f7.1,'MBps',/,
     &       30X,'io-nodes = ',i5,/)

      return
      end
c***********************************************************************
      subroutine mfo_singlev(ioflds,nout,offs0,stride,strideB,
     $     UR1,UR2,UR3,VECT)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFHEAT
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TIME_STEPPERD'
      INCLUDE 'ARNOLDI_ARPACKD'

!     argument list
      integer ioflds,nout
      integer*8 offs0,stride,strideB
      real UR1(LXO*LXO*LXO*LELT), UR2(LXO*LXO*LXO*LELT),
     $     UR3(LXO*LXO*LXO*LELT)
      real VECT(LVAS)

!     local variables
      integer*8 offs
!-----------------------------------------------------------------------
      offs = offs0 + ioflds*stride + NDIM*strideB
      call byte_set_view(offs,ifh_mbyte)

      call copy(UR1,VECT(1),NVECAV)
      call copy(UR2,VECT(1+NVECAV),NVECAV)
      if (IF3D) call copy(UR3,VECT(1+2*NVECAV),NVECAV)

      call mfo_outv(UR1,UR2,UR3,nout,NXO,NYO,NZO)
      ioflds = ioflds + NDIM

      if (IFHEAT) then
         offs = offs0 + ioflds*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         call copy(UR1,VECT(1+NDIM*NVECAV),NVECAT)

         call mfo_outs(UR1,nout,NXO,NYO,NZO)
         ioflds = ioflds + 1
      endif

      return
      end
c***********************************************************************
      subroutine mfi_singlev(iofldr,offs0,stride,strideB,
     $     UR1,UR2,UR3,VECT)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFHEAT
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TIME_STEPPERD'
      INCLUDE 'ARNOLDI_ARPACKD'
!     argument list
      integer iofldr
      integer*8 offs0,stride,strideB
      real UR1(LX1*LX1*LX1*LELT), UR2(LX1*LX1*LX1*LELT),
     $     UR3(LX1*LX1*LX1*LELT)
      real VECT(LVAS)

      integer lwk
      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      real wk(lwk)
      common /scrns/ wk

!     local variables
      integer*8 offs
!-----------------------------------------------------------------------
      offs = offs0 + iofldr*stride + NDIM*strideB
      call byte_set_view(offs,ifh_mbyte)
      call mfi_getv(UR1,UR2,UR3,wk,lwk,.false.)

      call copy(VECT(1),UR1,NVECAV)
      call copy(VECT(1+NVECAV),UR2,NVECAV)
      if (IF3D) call copy(VECT(1+2*NVECAV),UR3,NVECAV)
      iofldr = iofldr + NDIM

      if (IFHEAT) then
         offs = offs0 + iofldr*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         call mfi_gets(UR1,wk,lwk,.false.)
         call copy(VECT(1+NDIM*NVECAV),UR1,NVECAT)
         iofldr = iofldr + 1
      endif

      return
      end
c***********************************************************************
      subroutine arn_weight_fcn

!     Called when IDOARP.eq.2
!     if y=Ax is the Arnoldi vector
!     This subroutine calculates w = My
!     and stores it in WORDA(IPNTARP(2))
!     Needed by Arnoldi to perform M-weighted
!     orthogonalization

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'
      include 'MASS_DEF'
      include 'MASS'                ! BM1
      include 'INPUT_DEF'
      include 'INPUT'

      integer vecptr,aptr

      vecptr=0
      aptr  = IPNTARP(2)

!     velocity
      call copy(WORKDA(aptr+vecptr),BM1,NVECAV)
      vecptr=vecptr+NVECAV
      call copy(WORKDA(aptr+vecptr),BM1,NVECAV)
      vecptr=vecptr+NVECAV
      if (IF3D) then
        call copy(WORKDA(aptr+vecptr),
     $     BM1,NVECAV)
        vecptr=vecptr+NVECAV
      endif

!     Pressure. Zero weight. 
!     If we have included pressure in the arnoldi vector
!     and are using a Semi-norm
      if (TSTIFPR) then
        call rzero(WORKDA(aptr+vecptr),NVECAP)
        vecptr=vecptr+NVECAP
      endif  

!     Temperature. 
      if(IFHEAT) then
        call copy(WORKDA(aptr+vecptr),
     $       BM1,NVECAT)
             vecptr=vecptr+NVECAT
      endif

      call col2(WORKDA(IPNTARP(2)),WORKDA(IPNTARP(1)),NVECAS)

      return
      end subroutine arn_weight_fcn

!----------------------------------------------------------------------             

      subroutine arn_workda2flds

!     Copy fields from Arnoldi work array to Nek arrays            

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'
      include 'SOLN_DEF'
      include 'SOLN'

      integer vecptr
      integer aptr

!     Array pointers      
      vecptr=0
      aptr=IPNTARP(1)

!     Velocity
      call copy(VXP,WORKDA(aptr),NVECAV)
      vecptr=vecptr+NVECAV
      call copy(VYP,WORKDA(aptr+vecptr),NVECAV)
      vecptr=vecptr+NVECAV
      if (IF3D) then
        call copy(VZP,WORKDA(aptr+vecptr),NVECAV)
        vecptr=vecptr+NVECAV
      endif

!     Pressure 
      if (TSTIFPR) then
        call copy(PRP,WORKDA(aptr+vecptr),NVECAP)
        vecptr=vecptr+NVECAP
      endif  

!     temperature
      if (IFHEAT) then
        call copy(TP,WORKDA(aptr+vecptr),NVECAT)
        vecptr=vecptr+NVECAT
      endif  

      return
      end subroutine arn_workda2flds
!---------------------------------------------------------------------- 
      subroutine arn_flds2workda

!     Copy the vector resulting from
!     A*x to the arnoldi work array            

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'INPUT_DEF'
      include 'INPUT'           ! IFHEAT, IF3D
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, TP, V?MASK
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !

      integer vecptr
      integer aptr

      vecptr= 0
      aptr  = IPNTARP(2)      ! pointer for Arpack work array

      call col3(WORKDA(aptr),VXP,V1MASK,NVECAV)
      vecptr=vecptr+NVECAV
      call col3(WORKDA(aptr+vecptr),VYP,V2MASK,NVECAV)
      vecptr=vecptr+NVECAV
      if (IF3D) then
        call col3(WORKDA(aptr+vecptr),VZP,V3MASK,NVECAV)
        vecptr=vecptr+NVECAV
      endif  

!     if we include pressure in the arnoldi vector
!     No masks here
      if (TSTIFPR) then
        call copy(WORKDA(aptr+vecptr),PRP,NVECAP)
        vecptr=vecptr+NVECAP
      endif  

!     temperature
      if (IFHEAT) then
!       Not sure if we need a mask here. Needs to be tested.
        call col3(WORKDA(aptr+vecptr),TP,TMASK,NVECAT)
        vecptr=vecptr+NVECAT
      endif


      return
      end subroutine arn_flds2workda
!----------------------------------------------------------------------       
      subroutine arn_flds2resida

!     Copy the vector to the arnoldi resid array
!     Only done once during initialization
!     (Not entirely sure why its needed) 

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'INPUT_DEF'
      include 'INPUT'           ! IFHEAT, IF3D
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, TP, V?MASK
      include 'MASS_DEF'
      include 'MASS'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !
      include 'NLFSI'                   ! psi,psiv

      integer vecptr

      integer ntotv,ntott,ntotp

      real arn_glsc2_wt             ! function
      real beta


!     Apply masks      
      call opcol2(vxp,vyp,vzp,v1mask,v2mask,v3mask)

!     no masks for pressure

      if (ifheat) then
        ntott = nx1*ny1*nz1*nelt
        call col2(tp,tmask,ntott)
      endif  
   
!     Normalize starting vector
!     Its not really used so the normalization step is not necessary      
      beta = arn_glsc2_wt()
      beta = sqrt(beta)
      call opcmult(vxp,vyp,vzp,1./beta)

      if (TSTIFPR) then
        ntotp=nx2*ny2*nz2*nelv
        call cmult(prp,1./beta,ntotp)
      endif  

      if (ifheat) then
!       Not sure if we need a mask here. Needs to be tested.
        ntott = nx1*ny1*nz1*nelt
        call cmult(tp,1./beta,ntott)
      endif 
 
      vecptr=1
      call copy(RESIDA(vecptr),VXP,NVECAV)
      vecptr=vecptr+NVECAV
      call copy(RESIDA(vecptr),VYP,NVECAV)
      vecptr=vecptr+NVECAV
      if (IF3D) then
        call copy(RESIDA(vecptr),VZP,NVECAV)
        vecptr=vecptr+NVECAV
      endif  

!     if we include pressure in the arnoldi vector
      if (TSTIFPR) then
        call copy(RESIDA(vecptr),PRP,NVECAP)
        vecptr=vecptr+NVECAP
      endif 

!     temperature
      if (IFHEAT) then
        call copy(RESIDA(vecptr),TP,NVECAT)
        vecptr=vecptr+NVECAT
      endif


      return
      end subroutine arn_flds2resida
!----------------------------------------------------------------------       
      real function arn_glsc2_wt()

!     global weighted scalar product

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
      include 'SOLN_DEF'
      include 'SOLN'
      include 'MASS_DEF'
      include 'MASS'            ! VOLVM1, VOLTM1
!      include 'CONHT'           ! CHCST_SC,CHCFF_V,CHCFF_T

!     local variables
      integer ntotv, ntott
      real scal
      real op_glsc2_wt
      real glsc3  

      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      scal = 0.
      scal = op_glsc2_wt(vxp,vyp,vzp,vxp,vyp,vzp,bm1)

!     Zero weight for pressure

!     for conjugate heat transfer
      if (IFHEAT) then
        scal = scal + glsc3(tp,tp,bm1,ntott)
      endif  

      arn_glsc2_wt = scal 
      
      return
      end function arn_glsc2_wt

!---------------------------------------------------------------------- 

