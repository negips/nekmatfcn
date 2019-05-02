!----------------------------------------------------------------------
!
!     Name: Save 2D slices time history
!     Author: Prabal Negi
!     Description: Testing to save 2D slices to get time history
!     Last Modified: 21/12/2016
!
!---------------------------------------------------------------------- 
!----------------------------------------------------------------------

      subroutine SLICE_MAKE

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'

      integer icalld
      save icalld
      data icalld /0/


      if (icalld.eq.0) then
        icalld=icalld+1
        call SLICE_INIT
        return
      endif              

      call SLICE_SAVE

      return
      end subroutine SLICE_MAKE
!---------------------------------------------------------------------- 
      subroutine SLICE_MFO_OUTFLD

      implicit none

      include 'SIZE_DEF'        ! missing definitions in include files
      include 'SIZE'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SLICE'           ! 2D slice specific variables
      include 'TSTEP_DEF'
      include 'TSTEP'

!     local variables
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
      integer nxyo              ! Slice size    

      character*3 prefix        ! file prefix

      integer*8 offs0, offs     ! offset      
      integer*8 stride,strideB  ! stride

      integer ioflds            ! fields count

      real dnbyte               ! byte sum
      real tiostart, tio        ! simple timing

!     dummy arrays
      real ur1(lx1,ly1,SLICE_MAXSAVES,3*LELT)
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

!     get number of elements owned by proceesor with smaller nid
      itmp = SLICE_LOWN
      nel2DB = igl_running_sum(itmp)
      nel2DB = nel2DB - SLICE_LOWN
!     replace value
      nelBl = NELB
      NELB = nel2DB

!     set element size
      NXO   = lx1
      NYO   = ly1
      NZO   = SLICE_NSAVES
      nxyo  = NXO*NYO
      nxyzo = NXO*NYO*NZO

!     open files on i/o nodes
      prefix='slc'
      ierr=0
      if (NID.eq.PID0) call mfo_open_files(prefix,ierr)

      call err_chk(ierr,'Error; opening file in slice_mfo_outfld. $')

!     write header, byte key, global ordering
      call slice_mfo_write_hdr

!     initial offset: header, test pattern, global ordering
      offs0 = iHeaderSize + 4 + ISIZE*SLICE_GNUM
      offs = offs0

!     stride
      strideB =      NELB * nxyzo * WDSIZO
      stride  = SLICE_GNUM * nxyzo * WDSIZO

!     count fields
      ioflds = 0

!     write coordinates
      kl = 0
!     copy vector
      do il=1,NELV
        if (SLICE_IFSLICE(il)) then
          kl = kl +1
          do jl = 1,NZO
            call copy(ur1(1,1,jl,3*kl+1),SLICE_XM1(1,1,il),nxyo)
            call copy(ur1(1,1,jl,3*kl+2),SLICE_YM1(1,1,il),nxyo)
            call cfill(ur1(1,1,jl,3*kl+3),SLICE_TIME(jl),nxyo) 
          enddo
        endif 
      enddo
!     check consistency
      ierr = 0
      if (kl.ne.SLICE_LOWN) ierr=1
      call err_chk(ierr,'Error stat; inconsistent SLICE_LOWN.1 $')
!     offset
      kl = 3*kl
      offs = offs0 + stride*ioflds + 3*strideB         ! 3 flds: x,y,z
      call byte_set_view(offs,IFH_MBYTE)
      call mfo_outs(ur1,kl,NXO,NYO,NZO)
      ioflds = ioflds + 3

!     write fields
      do jl=1,SLICE_NVAR
         kl = 0
!        copy vector
         do il=1,NELV
            if(SLICE_IFSLICE(il)) then
               kl = kl +1
               call copy(ur1(1,1,1,kl),SLICE_VAR(1,1,1,il,jl),nxyzo)
            endif
         enddo
!        check consistency
         ierr = 0
         if (kl.ne.SLICE_LOWN) ierr=1
         call err_chk(ierr,'Error stat; inconsistent SLICE_LOWN.2 $')
!        offset
         offs = offs0 + stride*ioflds + strideB
         call byte_set_view(offs,IFH_MBYTE)
         call mfo_outs(ur1,kl,NXO,NYO,NZO)
         ioflds = ioflds + 1
      enddo

!     write averaging data
!      call stat_mfo_write_stat2D

!     count bytes
      dnbyte = 1.*ioflds*SLICE_LOWN*WDSIZO*nxyzo

      ierr = 0
      if (NID.eq.PID0) 
#ifdef MPIIO
     &     call byte_close_mpi(IFH_MBYTE,ierr)
#else
     &     call byte_close(ierr)
#endif
      call err_chk(ierr,
     $     'Error closing file in slice_mfo_outfld2D. Abort. $')

      tio = dnekclock_sync()-tiostart
      if (tio.le.0) tio=1.

      dnbyte = glsum(dnbyte,1)
      dnbyte = dnbyte + iHeaderSize + 4 + ISIZE*SLICE_GNUM
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

!     clean up array
!      il = STAT_LM1*STAT_LM1*LELT*STAT_NVAR
!     reset averagign parameters
!     to be added

!     update timing and counters
!      STAT_TIO = STAT_TIO + tio
!      STAT_ION = STAT_ION + 1

      return
      end subroutine SLICE_MFO_OUTFLD

!---------------------------------------------------------------------- 

!     based on stat_mfo_write_hdr2D
!     write hdr, byte key, global ordering
      subroutine slice_mfo_write_hdr

      implicit none

      include 'SIZE_DEF'        ! missing definitions in include files
      include 'SIZE'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'SLICE'           ! Slice speciffic variables


!     local variables
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
      nelo = SLICE_GNUM
#else
      nfileoo = NFILEO
      if(NID.eq.PID0) then                ! how many elements to dump
        nelo = SLICE_LOWN
        do jl = PID0+1,PID1
           mtype = jl
           call csend(mtype,idum,ISIZE,jl,0)   ! handshake
           call crecv(mtype,inelp,ISIZE)
           nelo = nelo + inelp
        enddo
      else
        mtype = NID
        call crecv(mtype,idum,ISIZE)          ! hand-shake
        call csend(mtype,SLICE_LOWN,ISIZE,PID0,0)   ! u4 :=: u8
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
         write(RDCODE1(3),'(I1)') SLICE_NVAR/10
         write(RDCODE1(4),'(I1)') SLICE_NVAR-(SLICE_NVAR/10)*10

!     two last variables added for statistics
         write(hdr,1) WDSIZO,NXO,NYO,NZO,nelo,SLICE_GNUM,TIME,999,
     $        FID0, nfileoo, (rdcode1(il),il=1,10)
 1       format('#std',1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,
     $        e20.13,1x,i9,1x,i6,1x,i6,1x,10a)

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
     $     'Error writing header in slice_mfo_write_hdr. $')

!     write global 2D elements numbering for this group
!     copy data
      lglist(0) = SLICE_LOWN
      kl = 0
      do il=1,nelv
        if (SLICE_IFSLICE(il)) then
           kl = kl +1
           lglist(kl) = SLICE_GMAP(il)
        endif
      enddo
!     check consistency
      ierr = 0
      if (kl.ne.SLICE_LOWN) ierr=1
      call err_chk(ierr,'Error stat; inconsistent SLICE_LOWN.3 $')

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

         len = ISIZE*(SLICE_LOWN+1)
         call csend(mtype,lglist,len,PID0,0)  
      endif 

      call err_chk(ierr,
     $     'Error writing global nums in slice_mfo_write_hdr. $')

      return
      end

!-----------------------------------------------------------------------

      subroutine SLICE_INIT

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'SLICE' 

      integer el
      integer ix,iy,iz
      integer nxyz
      integer cnt,itmp

      integer iglsum

      real lzmin,lzmax
      real vlmin,vlmax

      nxyz=nx1*ny1*nz1

      SLICE_Z0=0.05
      SLICE_COMP = 2

      cnt = 0
      do el=1,nelv
        lzmin=vlmin(zm1(1,1,1,el),nxyz)
        lzmax=vlmax(zm1(1,1,1,el),nxyz)
        if ((SLICE_Z0.gt.lzmin).and.(SLICE_Z0.lt.lzmax)) then
          SLICE_IFSLICE(el)=.TRUE.
          cnt=cnt+1  
        else
          SLICE_IFSLICE(el)=.FALSE.
        endif
      enddo 
      
      SLICE_LOWN=cnt
      
      itmp = iglsum(cnt,1)
      SLICE_GNUM = itmp
      itmp = itmp-SLICE_LOWN

      cnt=0
      do el=1,nelv
        if (SLICE_IFSLICE(el)) then
          cnt=cnt+1  
          SLICE_GMAP(el)=itmp+cnt              ! element mapping

!         X/Y            
          do ix=1,nx1
            do iy=1,ny1
              do iz=1,1
                SLICE_XM1(ix,iy,el)=XM1(ix,iy,iz,el)  
                SLICE_YM1(ix,iy,el)=YM1(ix,iy,iz,el)
              enddo
            enddo
          enddo    
        endif
      enddo 

      return
      end subroutine SLICE_INIT

!---------------------------------------------------------------------- 

      subroutine SLICE_SAVE

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SLICE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'SOLN_DEF'
      include 'SOLN'      

      integer npos
      
      if ((mod(ISTEP,SLICE_COMP).eq.0).and.(ISTEP.gt.0)) then

        if (nid.eq.0) write(6,*) 'Saving fld slice'

        SLICE_NSAVES=SLICE_NSAVES+1    

        npos=1    
        call SLICE_SAVE_FLD(SLICE_VAR(1,1,SLICE_NSAVES,1,npos),vx,
     $  SLICE_IFSLICE)

        npos=npos+1  
        call SLICE_SAVE_FLD(SLICE_VAR(1,1,SLICE_NSAVES,1,npos),vy,
     $  SLICE_IFSLICE)

        npos=npos+1    
        call SLICE_SAVE_FLD(SLICE_VAR(1,1,SLICE_NSAVES,1,npos),vz,
     $  SLICE_IFSLICE)

        SLICE_TIME(SLICE_NSAVES)=Time    

        if ((NSTEPS-ISTEP.lt.SLICE_COMP).or.
     $         (SLICE_NSAVES.EQ.SLICE_MAXSAVES)) then
          call SLICE_MFO_OUTFLD
          SLICE_NSAVES=0
        endif 
 
      endif



      return
      end subroutine SLICE_SAVE

!---------------------------------------------------------------------- 

      subroutine SLICE_SAVE_FLD(svar,fld,ifslice)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'

      integer cnt,el
      integer ix,iy,iz
      real svar(lx1,ly1,lelt)
      real fld(lx1,ly1,lz1,lelt)
      logical ifslice(lelt)

      cnt=0
      do el=1,nelv
        if (ifslice(el)) then
          cnt=cnt+1
          do ix=1,nx1
            do iy=1,ny1
              do iz=1,1
                svar(ix,iy,el)=fld(ix,iy,iz,el)
              enddo
            enddo
          enddo
        endif
      enddo            

      return
      end subroutine SLICE_SAVE_FLD

!---------------------------------------------------------------------- 










     
