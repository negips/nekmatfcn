c     Adam Peplinski; 20013.05.15
c     This file includes subroutines to support selective frequency
c     damping (SFD) in nekton
c
c WARNING: Mattias modified these vars when setting up the roughness case,
c          2014-03-13
c     Parameters used by this set of subroutines:
c     UPARAM(1) - /= 0 if restart
c     UPARAM(2) - checkpiont dump frequency (number of time steps)
c     UPARAM(4) - frequency of saving convergence rate and velocity probes
c     UPARAM(33) - 0 if no SFD
c     UPARAM(34) - period of disturbance
c     UPARAM(35) - forcing control parameter
c     UPARAM(36) - number of SFD checkpoint file to restart
c
!---------------------------------------------------------------------- 
!     read parameters selective frequency damping  
      subroutine sfd_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'SFD'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /SFD/ IFSFD,SFDD,SFDCHI,SFDFCONV,SFDIRST


!     default values
      IFSFD = .FALSE.         ! if sfd
      SFDD = 1.05000          ! sfd time-period.
      SFDCHI = 0.5            ! sfd strength
      SFDFCONV = 50           ! sfd convergence
      SFDIRST  = 1            ! SFD file no to restart from

!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=SFD,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading SFD parameters.$')

!     broadcast data
      call bcast(SFDD,     WDSIZE)
      call bcast(SFDCHI,   WDSIZE)
      call bcast(SFDFCONV, ISIZE)
      call bcast(IFSFD,    LSIZE)
      call bcast(SFDIRST,  ISIZE)

      return
      end
!---------------------------------------------------------------------- 
!     write parameters selective frequency damping 
      subroutine sfd_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'SFD'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /SFD/ IFSFD,SFDD,SFDCHI,SFDFCONV,SFDIRST

!     read the file
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=SFD,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing SFD parameters.$')

      return
      end
!-----------------------------------------------------------------------

      subroutine sfd_main
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP
      include 'SFD'

      if (ISTEP.eq.0) then      ! before first step
!     initialisation of SFD
         call sfd_init
        
      else
!     SFD evolution
         if (IFSFD) then
            call sfd_solve
            call sfd_rst_save
            call sfd_end
         endif
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine sfd_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'       !
      include 'INPUT'           ! IF3D
      include 'PARALLEL_DEF'
      include 'PARALLEL'        ! GLLEL
      include 'SFD'

!     argument list
      real ffx, ffy, ffz
      integer ix,iy,iz,ieg

!     local variables
      integer iel

      if (IFSFD) then
         iel=GLLEL(ieg)
         ffx = ffx - SFDCHI*BFSX(ix,iy,iz,iel)
         ffy = ffy - SFDCHI*BFSY(ix,iy,iz,iel)
         if (IF3D) ffz = ffz - SFDCHI*BFSZ(ix,iy,iz,iel)
      endif

      return
      end
!---------------------------------------------------------------------- 

      subroutine sfd_init
c     initialise all SFD variables

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TSTEP_DEF'
      include 'TSTEP'           ! NSTEP
      include 'SFD'
      include 'CHKPOINT'        ! ifchkptrst,chkptstep

      integer ntot1, ilag
      real dnekclock            ! function

c     Mattias bug fix :)
      ntot1 = nx1*ny1*nz1*nelv

c     should we perform SFD
      if(IFSFD) then

c     timing
         SFDTIME1=dnekclock()

c     check nekton parameters
         if (.not.IFTRAN) then
            if (NID.eq.0) write(6,*)
     $           'ERROR: SFD requres transient equations'
            call exitt
         endif

         if (NSTEPS.eq.0) then
            if (NID.eq.0) write(6,*)
     $           'ERROR: SFD requires NSTEPS>0'
            call exitt
         endif

         if (IFPERT) then
            if (NID.eq.0) write(6,*)
     $           "ERROR: SFD shouldn't be run in perturbation mode"
            call exitt
         endif

c     set SFD parameters
c     do we restart
         IFSFDRST = .FALSE.
         if (ifchkptrst) IFSFDRST = .TRUE.

c     checkpointing frequency
!         SFDSAVE = int(UPARAM(2))
         SFDSAVE = chkptstep
c     frequency for saving convegence frequency
!         SFDFCONV = int(UPARAM(4))
c     open file for saving convergence history
         if (NID.eq.0) open(150,file='SFDconv.out',status='unknown')

c     delta
         if (SFDD.gt.0.0) then
            SFDD = 8.0*ATAN(1.0)/SFDD
         else
            if(NID.eq.0) then
               write(6,*) 'SFD: ERROR'
               write(6,*) 'SFDD = ', SFDD
               write(6,*) 'Disturbance period must be positive.'
            endif
            call exitt
         endif

c     chi
         if (SFDCHI.eq.0.0) then
            if(NID.eq.0) then
               write(6,*) 'SFD: ERROR'
               write(6,*) 'SFDCHI = ', SFDCHI
               write(6,*) 'Forcing control must be positive.'
            endif
            call exitt
         endif

c     initialise arrays
c     PLACE FOR RESTART
         if (IFSFDRST) then
c     read checkpoint
            call sfd_rst_read
         else
            do ilag=1,3
               call rzero(VSXLAG(1,1,1,1,ilag),ntot1)
               call rzero(VSYLAG(1,1,1,1,ilag),ntot1)
               call rzero(VSZLAG(1,1,1,1,ilag),ntot1)
            enddo

            call opcopy (VSX,VSY,VSZ,VX,VY,VZ)
         endif

c     find the difference between V? and VS?
         call opsub3(BFSX,BFSY,BFSZ,VX,VY,VZ,VSX,VSY,VSZ)

c     print info
         if (NID.eq.0) then
            write(6,*)
            write(6,*) 'SFD initialised'
            write(6,*) 'Parameters:'
            write(6,'(A15,G13.5)') 'DELTA = ',1.0/SFDD
            write(6,'(A15,G13.5)') 'CHI = ',SFDCHI
         endif

c     timing
         SFDTIME2=dnekclock()
         SFDTIME = SFDTIME2 -SFDTIME1
      endif

      return
      end
c-----------------------------------------------------------------------
c     Following subroutine are based on
c     makeabf
c     makebdf
c     from navier1.f
c-----------------------------------------------------------------------
      subroutine sfd_solve
C-----------------------------------------------------------------------
C     
C     Sum up contributions to kth order extrapolation scheme and
C     get new filtered velocity field.
c     This subroutine is based on makeabf and makebdf from navier1.f
C     
C-----------------------------------------------------------------------

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D
c     user include files
      include 'SFD'

C     temporary storage
      real TA1,TA2,TA3
      COMMON /SCRUZ/ TA1 (LX1,LY1,LZ1,LELV)
     $     ,             TA2 (LX1,LY1,LZ1,LELV)
     $     ,             TA3 (LX1,LY1,LZ1,LELV)

c     storage of rhs
      real absx1(LX1,LY1,LZ1,LELV),absx2(LX1,LY1,LZ1,LELV),
     $     absy1(LX1,LY1,LZ1,LELV),absy2(LX1,LY1,LZ1,LELV),
     $     absz1(LX1,LY1,LZ1,LELV),absz2(LX1,LY1,LZ1,LELV)
      save absx1,absx2,absy1,absy2,absz1,absz2

      integer ntot1, ilag

      real ab0, ab1, ab2

      integer icalld
      save    icalld
      data    icalld  /0/

      real dnekclock                ! function
      real gl2norm

c     should we perform SFD
      if(IFSFD) then

c     timing
         SFDTIME1=dnekclock()

         ntot1 = NX1*NY1*NZ1*NELV

c     this is done only once
         if (icalld.eq.0) then
            icalld = icalld + 1

c     initialise arrays
            call rzero(absx1,ntot1)
            call rzero(absx2,ntot1)
            call rzero(absy1,ntot1)
            call rzero(absy2,ntot1)
            call rzero(absz1,ntot1)
            call rzero(absz2,ntot1)
         endif

c     A-B part
c     current rhs
c     I use BFS? vectors generated during the convergence tests
c     so skip it
c         call opsub3(BFSX,BFSY,BFSZ,VXLAG,VYLAG,VZLAG,VSX,VSY,VSZ)
c     finish rhs
         call opcmult(BFSX,BFSY,BFSZ,SFDD)

C     old timesteps
         ab0 = AB(1)
         ab1 = AB(2)
         ab2 = AB(3)
         call add3s2 (TA1,absx1,absx2,ab1,ab2,ntot1)
         call add3s2 (TA2,absy1,absy2,ab1,ab2,ntot1)
c     save rhs
         call copy   (absx2,absx1,ntot1)
         call copy   (absy2,absy1,ntot1)
         call copy   (absx1,BFSX,ntot1)
         call copy   (absy1,BFSY,ntot1)
c     current
         call add2s1 (BFSX,TA1,ab0,ntot1)
         call add2s1 (BFSY,TA2,ab0,ntot1)
         if (IF3D) then
            call add3s2 (TA3,absz1,absz2,ab1,ab2,ntot1)
            call copy   (absz2,absz1,ntot1)
            call copy   (absz1,BFSZ,ntot1)
            call add2s1 (BFSZ,TA3,ab0,ntot1)
         endif
c     multiplication by timestep
         call opcmult(BFSX,BFSY,BFSZ,DT)

c     BD part
         ab0 = BD(2)
         call opadd2cm(BFSX,BFSY,BFSZ,VSX,VSY,VSZ,ab0)
C     
         do ilag=2,NBD
            ab0 = BD(ilag+1)
            call opadd2cm(BFSX,BFSY,BFSZ,VSXLAG (1,1,1,1,ILAG-1),
     $           VSYLAG (1,1,1,1,ILAG-1),VSZLAG (1,1,1,1,ILAG-1),ab0)
         enddo
C     take into account restart option
         if (IFSFDRST.and.(ISTEP.lt.SFDNRSF))
     $        call opcopy (TA1,TA2,TA3
     $        ,VSXLAG(1,1,1,1,3),VSYLAG(1,1,1,1,3),VSZLAG(1,1,1,1,3))

C     Keep old filtered velocity fields
         do ilag=3,2,-1
            call opcopy(VSXLAG(1,1,1,1,ilag),VSYLAG(1,1,1,1,ilag),
     $           VSZLAG(1,1,1,1,ilag),
     $           VSXLAG(1,1,1,1,ilag-1),VSYLAG(1,1,1,1,ilag-1),
     $           VSZLAG(1,1,1,1,ilag-1))
         enddo

         call opcopy (VSXLAG,VSYLAG,VSZLAG,VSX,VSY,VSZ)

c     calculate new filtered velocity field
c     take into account restart option
         if (IFSFDRST.and.(ISTEP.lt.SFDNRSF)) then
            call opcopy (VSX,VSY,VSZ,TA1,TA2,TA3)
         else
            ab0 = 1.0/BD(1)
            call opcopy (VSX,VSY,VSZ,BFSX,BFSY,BFSZ)
            call opcmult(VSX,VSY,VSZ,ab0)
         endif

c     convergence test
c     find the difference between V? and VS?
         call opsub3(BFSX,BFSY,BFSZ,VX,VY,VZ,VSX,VSY,VSZ)

c     calculate L2 norms
         ab0 = gl2norm(BFSX,ntot1)
         ab1 = gl2norm(BFSY,ntot1)
         if (IF3D) ab2 = gl2norm(BFSZ,ntot1)
c     for tracking convergence
         if (NID.eq.0.and.mod(ISTEP,SFDFCONV).eq.0) then 
            if (IF3D) then
               write(150,'(4(1pE13.5))') TIME, ab0, ab1, ab2
            else
               write(150,'(3(1pE13.5))') TIME, ab0, ab1
            endif
            flush(150)
c     stamp logs
            write(6,*) 'SFD: Convergence (L2 norm per grid point):'
            write(6,'(A15,G13.5)') 'DVX = ',ab0
            write(6,'(A15,G13.5)') 'DVY = ',ab1
            if (IF3D) write(6,'(A15,G13.5)') 'DVZ = ',ab2
         endif

c     timing
         SFDTIME2=dnekclock()
         SFDTIME = SFDTIME + SFDTIME2 -SFDTIME1

      endif                     ! IFSFD

      return
      end
C
c-----------------------------------------------------------------------
c     create checkpoint
      subroutine sfd_rst_save

      implicit none
         
      include 'SIZE_DEF'
      include 'SIZE'            ! NID, NDIM, NPERT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP, NSTEPS, LASTEP
      include 'SFD'             !

      real VSX0(LX1,LY1,LZ1,LELV), VSY0(LX1,LY1,LZ1,LELV),
     $     VSZ0(LX1,LY1,LZ1,LELV)

      real VSXLAG0(LX1,LY1,LZ1,LELV,3), VSYLAG0(LX1,LY1,LZ1,LELV,3),
     $     VSZLAG0(LX1,LY1,LZ1,LELV,3)

      real dnekclock

c     save checkpoint
      if (IFSFD.and.mod(ISTEP,SFDSAVE).eq.(SFDNRSF-1)) then

c     timing
         SFDTIME1=dnekclock()

         if (NID.eq.0) write(6,*) 'SFD: writing checkpoint'

c     save filtered valocity field
         call mfo_sfd('SFD')

c     timing
         SFDTIME2=dnekclock()
         SFDTIME = SFDTIME + SFDTIME2 -SFDTIME1

      endif

      return
      end
c-----------------------------------------------------------------------
c     subroutine to read from checkpoints
      subroutine sfd_rst_read

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NID, NDIM, NPERT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP, NSTEPS, LASTEP
      include 'SFD'             !

C     temporary storage
      real TA1,TA2,TA3
      COMMON /SCRUZ/ TA1 (LX1,LY1,LZ1,LELV)
     $     ,             TA2 (LX1,LY1,LZ1,LELV)
     $     ,             TA3 (LX1,LY1,LZ1,LELV)

      integer ilag


      if (IFSFD) then

         if (NID.eq.0) write(6,*) 'SFD: reading checkpoint'

c     read filtered velocity field
         call mfi_sfd('SFD')

c     move velcity fields to sotre oldest one in VS?
         call opcopy (TA1,TA2,TA3
     $        ,VSXLAG(1,1,1,1,3),VSYLAG(1,1,1,1,3),VSZLAG(1,1,1,1,3))

         do ilag=3,2,-1
            call opcopy(VSXLAG(1,1,1,1,ilag),VSYLAG(1,1,1,1,ilag),
     $           VSZLAG(1,1,1,1,ilag),
     $           VSXLAG(1,1,1,1,ilag-1),VSYLAG(1,1,1,1,ilag-1),
     $           VSZLAG(1,1,1,1,ilag-1))
         enddo

         call opcopy (VSXLAG,VSYLAG,VSZLAG,VSX,VSY,VSZ)

         call opcopy (VSX,VSY,VSZ,TA1,TA2,TA3)


      endif

      return
      end
c-----------------------------------------------------------------------
c     subroutine to finalise SFD
      subroutine sfd_end

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NID, NDIM, NPERT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP, NSTEPS, LASTEP
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D
      include 'SFD'             ! BFS?
      include 'RESTART_DEF'
      include 'RESTART'
      include 'SOLN_DEF'
      include 'SOLN'

      integer ntot1
      real ab0, ab1, ab2

      logical lifxyo, lifpo, lifvo, lifto, lifreguo, lifpso(LDIMT1)

      integer i
      real gl2norm

      if (IFSFD.and.(ISTEP.eq.NSTEPS.or.LASTEP.eq.1)) then

c     final convergence
         ntot1 = NX1*NY1*NZ1*NELV
c     calculate L2 norms
         ab0 = gl2norm(BFSX,ntot1)
         ab1 = gl2norm(BFSY,ntot1)
         if (IF3D) ab2 = gl2norm(BFSZ,ntot1)

         if (NID.eq.0) then
            write(6,*) ''
            write(6,*) 'SFD: finalize'
            write(6,*) '   Time spent in SFD  ',SFDTIME
            write(6,*) '   Convergence (L2 norm per grid point):'
            write(6,'(A15,G13.5)') 'DVX = ',ab0
            write(6,'(A15,G13.5)') 'DVY = ',ab1
            if (IF3D) write(6,'(A15,G13.5)') 'DVZ = ',ab2
            write(6,*) ''
            write(6,*) 'SFD: saving velocity difference'
         endif

c     save the velocity difference for checking 
         lifxyo= IFXYO
         IFXYO = .TRUE.
         lifpo= IFPO
         IFPO = .FALSE.
         lifvo= IFVO
         IFVO = .TRUE.
         lifto= IFTO
         IFTO = .FALSE.
         do i=1,LDIMT1
            lifpso(i)= IFPSO(i)
            IFPSO(i) = .FALSE.
         enddo

         call outpost2(BFSX,BFSY,BFSZ,PR,T,0,'VDF')

         IFXYO = lifxyo
         IFPO = lifpo
         IFVO = lifvo
         IFTO = lifto
         do i=1,LDIMT1
            IFPSO(i) = lifpso(i)
         enddo

c     close file with convergence history
         if (NID.eq.0) close(150)

      endif

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     following subroutines are baseds on 
c     
c     mfo_outfld
c     from prepost.f;
c     mfi
c     from ic.f
c-----------------------------------------------------------------------
      subroutine mfo_sfd(prefix)  ! muti-file output

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
!      include 'TOTAL'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'SFD'

      integer i,ioflds,ierr,nout
      integer*8 offs0,offs,nbyte,stride,strideB,nxyzo8, lwdsizo

      character*3 prefix

      logical lifxyo, lifpo, lifvo, lifto, lifreguo, lifpso(LDIMT1)

      integer dnbyte
      real tiostart,tio
      real dnekclock_sync

      real glsum

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

      if (nid.eq.pid0) then
         call mfo_open_files(prefix,ierr)         ! open files on i/o nodes
      endif

      call mfo_write_hdr                     ! create element mapping +
                                             ! write hdr
      nxyzo8  = NXO*NYO*NZO
      strideB = nelB * nxyzo8*WDSIZO
      stride  = nelgt* nxyzo8*WDSIZO

      ioflds = 0
      ! dump all fields based on the t-mesh to avoid different
      ! topologies in the post-processor

c     current filtered velocity field
      offs = offs0 + NDIM*strideB
      call byte_set_view(offs,ifh_mbyte)

      call mfo_outv(VSX,VSY,VSZ,nout,NXO,NYO,NZO)
      ioflds = ioflds + NDIM

c     history
      do i=1,3
         offs = offs0 + ioflds*stride + NDIM*strideB
         call byte_set_view(offs,ifh_mbyte)

         call mfo_outv(VSXLAG(1,1,1,1,i),VSYLAG(1,1,1,1,i),
     $           VSZLAG(1,1,1,1,i),nout,NXO,NYO,NZO)
         ioflds = ioflds + NDIM
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
      call nekgsync
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
      subroutine mfi_sfd(prefix)

      implicit none

      include 'SIZE_DEF'            
      include 'SIZE'
!      include 'TOTAL'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'SFD'
      include 'CHKPOINT'
!      INCLUDE 'USER_PAR'        ! UPARAM

      character*3 prefix

      character*132  fname, bname
      character*1    fname1(132)
      equivalence   (fname1,fname)

      character*6  str

      character*132 hdr

      integer lwk
      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      real wk
      common /scrns/ wk(lwk)
      integer i, e, k, iofldsr, ierr

      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8
      integer len

      real dnbyte
      real tio, tiostart
      real dnekclock                ! function
      integer ltrunc                ! function
      real glsum                    ! function
      

      tiostart=dnekclock()

c     create file name
      call blank(fname,132)
      call blank(bname,132)

      len = ltrunc(SESSION,132) !  Add SESSION
      call chcopy(bname,SESSION,len)

!      call mfo_filename(prefix,fname,bname,k)
      call IO_mfo_fname(prefix,fname,bname,k)
      write(str,5) SFDIRST
 5    format(i5.5)
      call chcopy(fname1(k),str,5)

      call mfi_prepare(fname)       ! determine reader nodes +
                                    ! read hdr + element mapping 

      offs0   = iHeadersize + 4 + ISIZE*NELGR
      nxyzr8  = NXR*NYR*NZR
      strideB = nelBr* nxyzr8*WDSIZR
      stride  = nelgr* nxyzr8*WDSIZR


c     read arrays
      iofldsr = 0
c     filtered velocity
      offs = offs0 + NDIM*strideB
      call byte_set_view(offs,ifh_mbyte)
      call mfi_getv(VSX,VSY,VSZ,wk,lwk,.false.)

      iofldsr = iofldsr + NDIM

c     history
      do i=1,3
         offs = offs0 + iofldsr*stride + NDIM*strideB
         call byte_set_view(offs,ifh_mbyte)
         call mfi_getv(VSXLAG(1,1,1,1,i),VSYLAG(1,1,1,1,i),
     $           VSZLAG(1,1,1,1,i),wk,lwk,.false.)

         iofldsr = iofldsr + NDIM
      enddo

      nbyte = 0
      if(nid.eq.pid0r) nbyte = iofldsr*nelr*wdsizr*nxr*nyr*nzr

c     close files
#ifdef MPIIO
      if (nid.eq.pid0r) call byte_close_mpi(ifh_mbyte, ierr)
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
c-----------------------------------------------------------------------
