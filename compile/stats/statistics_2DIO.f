!=======================================================================
! Name        : statistics_2DIO
! Author      : Adam Peplinski
! Version     : last modification 2015.05.20
! Copyright   : GPL
! Description : This is a set of routines to write 2D statistics to 
!     the file. They are modiffication of the existing nek5000 routines
!=======================================================================
!     this is just modification of mfo_outfld
!     muti-file output
      subroutine stat_mfo_outfld2D(iftimedep)

      implicit none

      include 'SIZE_DEF'        ! missing definitions in include files
      include 'SIZE'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'STATS'           ! 2D statistics speciffic variables
      include 'TSTEP_DEF'
      include 'TSTEP'

!     local variables
      logical iftimedep         ! if time dependent    

!     temporary variables to overwrite global values
      logical ifreguol          ! uniform mesh
      logical ifxyol            ! write down mesh
      integer wdsizol           ! store global wdsizo
      integer nel2DB            ! running sum for owned 2D elements
      integer nelBl             ! store global nelB

      integer il, jl, kl        ! loop index
      integer itmp              ! dummy integer
      integer ierr              ! error mark
      integer nxyzo             ! element size

      character*3 prefix        ! file prefix

      integer*8 offs0, offs     ! offset      
      integer*8 stride,strideB  ! stride

      integer ioflds            ! fields count

      real dnbyte               ! byte sum
      real tiostart, tio        ! simple timing

!     dummy arrays
      real ur1(STAT_LM1,STAT_LM1,2*LELT)
      common /SCRUZ/  ur1

!     functions
      integer igl_running_sum
      real dnekclock_sync, glsum

!     simple timing
      tiostart=dnekclock_sync()

!     save and set global IO variables
!     no uniform mesh
      ifreguol = IFREGUO
      IFREGUO = .FALSE.

!     save mesh
      ifxyol = IFXYO
      IFXYO = .TRUE.

!     force double precission
      wdsizol = WDSIZO
!     for testing
      WDSIZO = WDSIZE

!     get number of 2D elements owned by proceesor with smaller nid
      itmp = STAT_LOWN
      nel2DB = igl_running_sum(itmp)
      nel2DB = nel2DB - STAT_LOWN
!     replace value
      nelBl = NELB
      NELB = nel2DB

!     set element size
      if (if3d) then     ! #2D
           NXO   = STAT_NM2
           NYO   = STAT_NM3
           NZO   = 1
           nxyzo = NXO*NYO*NZO
      else
          NXO    = STAT_NM1
          NYO    = STAT_NM2
          NZO    = 1
          nxyzo  = NXO*NYO*NZO
      endif

!     open files on i/o nodes
      if (iftimedep) then
        prefix = 'stt'
      else
        prefix='sts'
      endif
      ierr=0
      if (NID.eq.PID0) call mfo_open_files(prefix,ierr)

      call err_chk(ierr,'Error; opening file in stat_mfo_outfld2D. $')

!     write header, byte key, global ordering
      call stat_mfo_write_hdr2D(iftimedep)

!     initial offset: header, test pattern, global ordering
      offs0 = iHeaderSize + 4 + ISIZE*STAT_GNUM
      offs = offs0

!     stride
      strideB =      NELB * nxyzo * WDSIZO
      stride  = STAT_GNUM * nxyzo * WDSIZO

!     count fields
      ioflds = 0

!     write coordinates
      kl = 0
!     copy vector
      do il=1,STAT_LNUM
         if(STAT_OWN(il).eq.NID) then
            if (iftimedep) then
              call copy(ur1(1,1,2*kl+1),STAT_XMT(1,1,il),nxyzo)
              call copy(ur1(1,1,2*kl+2),STAT_YMT(1,1,il),nxyzo)
            else  
              call copy(ur1(1,1,2*kl+1),STAT_XM1(1,1,il),nxyzo)
              call copy(ur1(1,1,2*kl+2),STAT_YM1(1,1,il),nxyzo)
            endif
            kl = kl +1
         endif
      enddo
!     check consistency
      ierr = 0
      if (kl.ne.STAT_LOWN) ierr=1
      call err_chk(ierr,'Error stat; inconsistent STAT_LOWN.1 $')
!     offset
      kl = 2*kl
      offs = offs0 + stride*ioflds + 2*strideB         ! 2 flds: x,y
      call byte_set_view(offs,IFH_MBYTE)
      call mfo_outs(ur1,kl,NXO,NYO,NZO)
      ioflds = ioflds + 2

!     write fields
      do jl=1,STAT_NVAR
        kl = 0
!       copy vector
        do il=1,STAT_LNUM
          if(STAT_OWN(il).eq.NID) then
            kl = kl +1
            if (iftimedep) then
              call copy(ur1(1,1,kl),STAT_RAVGT(1,1,il,jl),
     $                nxyzo)
            else
              call copy(ur1(1,1,kl),STAT_RUAVG(1,1,il,jl),
     $                nxyzo)
            endif  
          endif
        enddo
!       check consistency
        ierr = 0
        if (kl.ne.STAT_LOWN) ierr=1
        call err_chk(ierr,'Error stat; inconsistent STAT_LOWN.2 $')
!       offset
        offs = offs0 + stride*ioflds + strideB
        call byte_set_view(offs,IFH_MBYTE)
        call mfo_outs(ur1,kl,NXO,NYO,NZO)
        ioflds = ioflds + 1
      enddo

!     write averaging data
      call stat_mfo_write_stat2D

!     count bytes
      dnbyte = 1.*ioflds*STAT_LOWN*WDSIZO*nxyzo

      ierr = 0
      if (NID.eq.PID0) 
#ifdef MPIIO
     &     call byte_close_mpi(IFH_MBYTE,ierr)
#else
     &     call byte_close(ierr)
#endif
      call err_chk(ierr,
     $     'Error closing file in stat_mfo_outfld2D. Abort. $')

      tio = dnekclock_sync()-tiostart
      if (tio.le.0) tio=1.

      dnbyte = glsum(dnbyte,1)
      dnbyte = dnbyte + iHeaderSize + 4 + ISIZE*STAT_GNUM
      dnbyte = dnbyte/1024/1024
      if(NIO.eq.0) write(6,7) ISTEP,TIME,dnbyte,dnbyte/tio,
     &     NFILEO
    7 format(/,i9,1pe12.4,' done :: Write checkpoint',/,
     &     30X,'file size = ',3pG12.2,'MB',/,
     &     30X,'avg data-throughput = ',0pf7.1,'MB/s',/,
     &     30X,'io-nodes = ',i5,/)

!     set global IO variables back
      IFREGUO = ifreguol
      IFXYO = ifxyol
      WDSIZO = wdsizol
      NELB = nelBl

!     clean up variables array
      il = STAT_LM1*STAT_LM1*LELT*STAT_NVAR
      if (iftimedep) then
        call rzero(STAT_RAVGT,il)
        il = STAT_LM1*STAT_LM1*LELT
        call rzero(STAT_XMT,il)
        call rzero(STAT_YMT,il)
      else
        call rzero(STAT_RUAVG,il)
      endif

!     reset averagign parameters
!     to be added

!     update timing and counters
      STAT_TIO = STAT_TIO + tio
      STAT_ION = STAT_ION + 1

      return
      end
c-----------------------------------------------------------------------
!     based on mfo_write_hdr
!     write hdr, byte key, global ordering
      subroutine stat_mfo_write_hdr2D(iftimedep)

      implicit none

      include 'SIZE_DEF'        ! missing definitions in include files
      include 'SIZE'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'STATS'           ! 2D statistics speciffic variables


!     local variables
      logical iftimedep         ! if time dependent stats
      real tstart_tmp           ! To specify which start time
      integer evrn_tmp          ! No of samples 

      real*4 test_pattern       ! byte key
      integer lglist(0:LELT)    ! dummy array
      common /ctmp0/ lglist
      integer idum, inelp
      integer nelo              ! number of elements to write
      integer nfileoo           ! number of files to create
      
      integer il, jl, kl        ! loop index
      integer mtype             ! tag

      integer ierr              ! error mark
      integer ibsw_out, len
      integer*8 ioff            ! offset

      character*132 hdr         ! header

#ifdef MPIIO
      nfileoo = 1   ! all data into one file
      nelo = STAT_GNUM
#else
      nfileoo = NFILEO
      if(NID.eq.PID0) then                ! how many elements to dump
        nelo = STAT_LOWN
        do jl = PID0+1,PID1
           mtype = jl
           call csend(mtype,idum,ISIZE,jl,0)   ! handshake
           call crecv(mtype,inelp,ISIZE)
           nelo = nelo + inelp
        enddo
      else
        mtype = NID
        call crecv(mtype,idum,ISIZE)          ! hand-shake
        call csend(mtype,STAT_LOWN,ISIZE,PID0,0)   ! u4 :=: u8
      endif 
#endif

!     write header
      ierr = 0
      if(NID.eq.PID0) then
         call blank(hdr,132)

!     varialbe set
         call blank(RDCODE1,10)

!     we save coordinates
         RDCODE1(1)='X'
!     and set of fields marked as passive scalars
         RDCODE1(2) = 'S'
         write(RDCODE1(3),'(I1)') STAT_NVAR/10
         write(RDCODE1(4),'(I1)') STAT_NVAR-(STAT_NVAR/10)*10

!     two last variables added for statistics
         if (iftimedep) then
           tstart_tmp = T_STAT_TSTART
           evrn_tmp = T_STAT_AVGITER 
         else
           tstart_tmp = STAT_TSTART
           evrn_tmp = STAT_EVRN 
         endif
         write(hdr,1) WDSIZO,NXO,NYO,NZO,nelo,STAT_GNUM,TIME,evrn_tmp,
     $        FID0, nfileoo, (rdcode1(il),il=1,10),
     $        tstart_tmp, STAT_TACC
 1       format('#std',1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,
     $        e20.13,1x,i9,1x,i6,1x,i6,1x,10a,1x,e17.10,1x,e17.10)

!     if we want to switch the bytes for output
!     switch it again because the hdr is in ASCII
         call get_bytesw_write(ibsw_out)
c      if (ibsw_out.ne.0) call set_bytesw_write(ibsw_out)
         if (ibsw_out.ne.0) call set_bytesw_write(0)  

!     write test pattern for byte swap
         test_pattern = 6.54321 

#ifdef MPIIO
! only rank0 (pid00) will write hdr + test_pattern + time list
         call byte_write_mpi(hdr,iHeaderSize/4,PID00,IFH_MBYTE,ierr)
         call byte_write_mpi(test_pattern,1,PID00,IFH_MBYTE,ierr)
#else
         call byte_write(hdr,iHeaderSize/4,ierr)
         call byte_write(test_pattern,1,ierr)
#endif

      endif

      call err_chk(ierr,
     $     'Error writing header in stat_mfo_write_hdr2D. $')

!     write global 2D elements numbering for this group
!     copy data
      lglist(0) = STAT_LOWN
      kl = 0
      do il=1,STAT_LNUM
         if(STAT_OWN(il).eq.NID) then
            kl = kl +1
            lglist(kl) = STAT_GMAP(il)
         endif
      enddo
!     check consistency
      ierr = 0
      if (kl.ne.STAT_LOWN) ierr=1
      call err_chk(ierr,'Error stat; inconsistent STAT_LOWN.3 $')

      if(NID.eq.PID0) then
#ifdef MPIIO
         ioff = iHeaderSize + 4 + NELB*ISIZE
         call byte_set_view (ioff,IFH_MBYTE)
         call byte_write_mpi (lglist(1),lglist(0),-1,IFH_MBYTE,ierr)
#else
         call byte_write(lglist(1),lglist(0),ierr)
#endif
         do jl = PID0+1,PID1
            mtype = jl
            call csend(mtype,idum,ISIZE,jl,0) ! handshake
            len = ISIZE*(LELT+1)
            call crecv(mtype,lglist,len)
            if(ierr.eq.0) then
#ifdef MPIIO
               call byte_write_mpi
     $              (lglist(1),lglist(0),-1,IFH_MBYTE,ierr)
#else
               call byte_write(lglist(1),lglist(0),ierr)
#endif
            endif
         enddo
      else
         mtype = NID
         call crecv(mtype,idum,ISIZE) ! hand-shake

         len = ISIZE*(STAT_LOWN+1)
         call csend(mtype,lglist,len,PID0,0)  
      endif 

      call err_chk(ierr,
     $     'Error writing global nums in stat_mfo_write_hdr2D. $')

      return
      end
c-----------------------------------------------------------------------
!     write additional data at the end of the file
      subroutine stat_mfo_write_stat2D

      implicit none

      include 'SIZE_DEF'        ! missing definitions in include files
      include 'SIZE'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'STATS'           ! 2D statistics specific variables

!     to be added
      if(NID.eq.PID0) then
#ifdef MPIIO
! only rank0 (pid00) will write hdr + test_pattern + time list
!         call byte_write_mpi(hdr,iHeaderSize/4,PID00,IFH_MBYTE,ierr)
!         call byte_write_mpi(test_pattern,1,PID00,IFH_MBYTE,ierr)
#else
!         call byte_write(hdr,iHeaderSize/4,ierr)
!         call byte_write(test_pattern,1,ierr)
#endif

      endif


      return
      end
c-----------------------------------------------------------------------
!     this is just modification of mfi
!     muti-file input
      subroutine stat_mfi_dat(nfile)

      implicit none

      include 'SIZE_DEF'        ! missing definitions in include files
      include 'SIZE'
      include 'RESTART_DEF'
      include 'RESTART'         ! IHEADERSIZE,WDSIZR,NXR,NYR,NZR,NELR,
                                ! NELGR,TIMER,ISTPR,IFILER,NFILER,RDCODE,NELBR,LELR,ER
      include 'PARALLEL_DEF'
      include 'PARALLEL'        ! ISIZE, WDSIZE
      include 'INPUT_DEF'
      include 'INPUT'           ! SESSION
      include 'STATS'           ! 2D statistics speciffic variables

!     argument list
      integer nfile             ! file number

!     local variables
      character*132 fname, bname ! file name and base name
      character*1 fnam1(132)
      equivalence (fnam1,fname)
      character*3 prefix        ! file name prefix
      character*6 str           ! string for file number
      character*132 hdr         ! file header
      character*4 dummy         ! dummy string for reading header

      integer len               ! string length
      integer ierr              ! error flag

      real deltat, dacct        !

      real*4 bytetest           ! prabal 

!     functions
      integer ltrunc

!     get the file name
!     clean strings
      call blank(fname,132)
      call blank(bname,132)

!     add SESSION
      len = ltrunc(SESSION,132) 
      call chcopy(bname,SESSION,len)

!     prefix and name
      prefix = 'sts'
      call IO_mfo_fname(prefix,fname,bname,len)
!     is fname too long?
      if ((len+5).gt.80) then
         if(NIO.eq.0) write(6,*)
     $        'Error: stat_mfi_2D; too long file name'
         call exitt
      endif

!     add file number
      write(str,"(i5.5)") nfile
      call chcopy(fnam1(len),str,5)

!     read data form the header that will not be extracted by
!     mfi_parse_hdr
!     mater node only
      ierr = 0
      if (NID.eq.0) then
         call mbyte_open(fname,0,ierr) ! open  blah000.fldnn
         if(ierr.ne.0) goto 10
         call blank     (hdr,IHEADERSIZE)
         call byte_read (hdr, IHEADERSIZE/4,ierr)
         if(ierr.ne.0) goto 10
         call byte_read (bytetest,1,ierr)
         if(ierr.ne.0) goto 10

         read(hdr,*,err=5) dummy
     $         ,  WDSIZR,NXR,NYR,NZR,NELR,NELGR,TIMER,STAT_EVRN
     $         ,  IFILER,NFILER
     $         ,  RDCODE, STAT_TSTART, STAT_TACC

!     close the file
         call byte_close(ierr)
         goto 10
!     mark reading error
 5       ierr = 1
      endif
 
 10   continue
!     broadcast statistics data
      call bcast(STAT_EVRN,ISIZE)
      call bcast(STAT_TSTART,WDSIZE)
      call bcast(STAT_TACC,WDSIZE)

      return
      end
