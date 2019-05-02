!====================================================================== 
!    Name           : statistics_pdf.f
!    Author         : Prabal S. Negi
!    Last Modified  : Apr 04, 2017 
!    Description    : PDF for defined surfaces 
!    Notes          : Requires the statistics_surf.f developed by
!                   : Prabal S. Negi
!======================================================================  
!----------------------------------------------------------------------

      subroutine surf_save_fld

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SURF_STATS'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'TOPOL_DEF'
      include 'TOPOL'
      include 'TSTEP_DEF'
      include 'TSTEP'

      real tmp_fld(lx1,ly1,lz1,lelt,7)
         
      common /scrns/ tmp_fld

      integer isf,sel,iel,iface
      integer KX1,KX2,KY1,KY2,KZ1,KZ2
      integer ix,iy,iz
      integer cnt

      real alpha,beta

      integer icalld
      save icalld
      data icalld /0/

      real snx,sny,snz        ! surface normals
      real stx,sty,stz        ! surface tangents
      integer f

      real vtmp


      if (icalld.eq.0) then
        surf_nsamples=0
        icalld=icalld+1
        call rzero(surf_fld,lx1*ly1*lelt*MAXSURFS*SURF_NFLDS)
      endif

!      SURF_FLDOUT = 5000

      call gradm1(tmp_fld(1,1,1,1,1),tmp_fld(1,1,1,1,2),
     $             tmp_fld(1,1,1,1,3),vx,nelt)

      call gradm1(tmp_fld(1,1,1,1,4),tmp_fld(1,1,1,1,5),
     $             tmp_fld(1,1,1,1,6),vy,nelt)
     
      surf_nsamples = surf_nsamples + 1

!     pdf fields
      do isf=1,nsurfs
        cnt=0
        do sel=1,SURF_MEMCOUNT(isf)
          iel=SURF_OBJ(1,sel,isf)
          iface=SURF_OBJ(2,sel,isf)

          CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $                  ,iface)
           
          do ix=kx1,kx2
            do iy=ky1,ky2
              do iz=kz1,kz2
                cnt=cnt+1

!               all are negative because we want inward normals                  
                f = eface1(iface)
                if (f.eq.1.or.f.eq.2) then         ! r face
                  snx = -unx(iy,iz,iface,iel)
                  sny = -uny(iy,iz,iface,iel)
                  snz = -unz(iy,iz,iface,iel)
                elseif (f.eq.3.or.f.eq.4) then     ! s face
                  snx = -unx(ix,iz,iface,iel)
                  sny = -uny(ix,iz,iface,iel)
                  snz = -unz(ix,iz,iface,iel)
                elseif (f.eq.5.or.f.eq.6) then     ! t face
                  snx = -unx(ix,iy,iface,iel)
                  sny = -uny(ix,iy,iface,iel)
                  snz = -unz(ix,iy,iface,iel)
                endif

                stx = -sny
                sty = snx
                stz = 0.                           ! we know it is 2D.
                if (stx.lt.0) then
                  stx = -stx
                  sty = -sty
                endif  

!               This should be normal gradient of tangential velocity                  
                vtmp = tmp_fld(ix,iy,iz,iel,1)*stx*snx + 
     $                  tmp_fld(ix,iy,iz,iel,2)*stx*sny
                vtmp = vtmp +  tmp_fld(ix,iy,iz,iel,4)*sty*snx + 
     $                  tmp_fld(ix,iy,iz,iel,5)*sty*sny

!               usual fields  
                if ((mod(ISTEP,SURF_FLDOUT)).eq.0.and.
     $             (ISTEP.gt.0)) then  
                  surf_fldx(cnt,isf)=XM1(ix,iy,iz,iel)
                  surf_fldy(cnt,isf)=YM1(ix,iy,iz,iel)
                  surf_fldz(cnt,isf)=ZM1(ix,iy,iz,iel)

                  surf_fld (cnt,1,isf)=snx
                  surf_fld (cnt,2,isf)=sny
                
                  surf_fld (cnt,3,isf)=vtmp

                endif  

                if (vtmp.le.0) then
                  vtmp = 1.
                else
                  vtmp = 0.  
                endif 
                surf_fld (cnt,4,isf) = surf_fld(cnt,4,isf)+vtmp

!               normalize pdf at the last step 
                if ((mod(ISTEP,SURF_FLDOUT)).eq.0.and.
     $             (ISTEP.gt.0)) then
                  surf_fld(cnt,4,isf)=
     $                   surf_fld(cnt,4,isf)/(surf_nsamples+0.)
                endif   
                  
              enddo     ! iz
            enddo       ! iy
          enddo         ! ix

        enddo           ! sel
      enddo             ! isf
       

      if ((mod(ISTEP,SURF_FLDOUT)).eq.0.and.
     $             (ISTEP.gt.0)) then  
        do isf = 1,NSURFS    
          call surf_mfo_outvtk(isf)
        enddo
        surf_nsamples = 0
        call rzero(surf_fld,lx1*ly1*lelt*MAXSURFS*SURF_NFLDS)
      endif    

      return
      end subroutine surf_save_fld 

!----------------------------------------------------------------------

      subroutine surf_mfo_outvtk(isf)

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
      include 'SURF_STATS'

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
      integer*8 offs_sc         ! scalar headers
      integer*8 stride_cr,stride_crB      ! strides for carriage return 

      integer ioflds            ! fields count

      real dnbyte               ! byte sum
      real tiostart, tio        ! simple timing

!     dummy arrays
      real ur1(lx1*ly1*LELT*3)
      common /SCRUZ/  ur1

!     functions
      integer igl_running_sum
      real dnekclock_sync, glsum
      integer ivlsum
      integer iglsum

      integer glitmp
      integer j
      integer offs_cr
      integer isf

!      integer vtkHeaderSize

!      real tmp_surf_fld(lx1*ly1*LELT,MAXSURFS)
!      real tmp_surf_fldx(lx1*ly1*LELT,MAXSURFS)
!      real tmp_surf_fldy(lx1*ly1*LELT,MAXSURFS)
!      real tmp_surf_fldz(lx1*ly1*LELT,MAXSURFS)
!
!      COMMON /SRFXYZU/ tmp_surf_fld,tmp_surf_fldx,
!     $       tmp_surf_fldy,tmp_surf_fldz 

      character*32 vtkpts
      character ch


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
      itmp = SURF_MEMCOUNT(isf)                  ! prabal. Only doing 1 surface for now
      glitmp = iglsum(itmp,1)
      nel2DB = igl_running_sum(itmp)
      nel2DB = nel2DB - itmp
!     replace value
      nelBl = NELB
      NELB = nel2DB

!     set element size
      NXO   = lx1
      NYO   = ly1
      if (ndim.eq.2) NYO = 1 
      NZO   = 1
      nxyo  = NXO*NYO
      nxyzo = NXO*NYO*NZO

!     open files on i/o nodes
      prefix='vtk'
      ierr=0
      if (NID.eq.PID0) call mfo_open_files(prefix,ierr)

      call err_chk(ierr,'Error; opening file in surf_mfo_outvtk. $')

!     write header, byte key, global ordering
      call surf_mfo_write_hdr(isf)            ! prabal. needs to be written

!     initial offset: header, test pattern, global ordering
!      vtkHeaderSize = 27+17+7+26+18+16+8+1+14+4+1+21       ! see surf_mfo_vtk_hdr routine
      offs0 = iHeaderSize + 4 + ISIZE*glitmp 
      offs = offs0
      offs_cr = 0

!     stride
      strideB =      NELB * nxyzo * WDSIZO
      stride  = glitmp * nxyzo * WDSIZO
!      stride_cr  = glitmp * nxyzo * 1
!      stride_crB  = NELB * nxyzo * 1

!     count fields
      ioflds = 0

!     write coordinates
      kl = 0
!     copy vector
      do kl=1,SURF_MEMCOUNT(isf)
        do il = 1,nxyzo
          ur1(3*(kl-1)*nxyzo+(il-1)*3+1) =
     $                   surf_fldx((kl-1)*nxyzo+il,isf)
          ur1(3*(kl-1)*nxyzo+(il-1)*3+2) =
     $                   surf_fldy((kl-1)*nxyzo+il,isf)
          ur1(3*(kl-1)*nxyzo+(il-1)*3+3) = 
     $                   surf_fldz((kl-1)*nxyzo+il,isf)
        enddo 
      enddo
      kl = 3*SURF_MEMCOUNT(isf)
!     check consistency
      ierr = 0
!      if (kl.ne.SLICE_LOWN) ierr=1
!      call err_chk(ierr,'Error stat; inconsistent SLICE_LOWN.1 $')
!     offset
      offs = offs0 + stride*ioflds + 3*strideB  ! + stride_cr + stride_crB 
      call byte_set_view(offs,IFH_MBYTE)
      call mfo_outs(ur1,kl,NXO,NYO,NZO)
      ioflds = ioflds + 3

!     write fields
      do jl=1,SURF_NFLDS
        kl = 0
!       copy vector
        do il=1,SURF_MEMCOUNT(isf)
          kl = kl +1
          do j=1,nxyzo 
            ur1((kl-1)*nxyzo+j) =
     $            surf_fld((kl-1)*nxyzo+j,jl,isf)
          enddo 
        enddo
!       check consistency
        ierr = 0
!        if (kl.ne.SLICE_LOWN) ierr=1
!        call err_chk(ierr,'Error stat; inconsistent SLICE_LOWN.2 $')
!       offset
        offs=offs0+stride*ioflds+strideB !+offs_sc+(ioflds-2)*stride_cr
!     $             +stride_crB
        call byte_set_view(offs,IFH_MBYTE)
        call mfo_outs(ur1,kl,NXO,NYO,NZO)
        ioflds = ioflds + 1
      enddo

!     write averaging data
!      call stat_mfo_write_stat2D

!     count bytes
      dnbyte = 1.*ioflds*glitmp*WDSIZO*nxyzo

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

!      dnbyte = glsum(dnbyte,1)
      dnbyte = dnbyte + iHeaderSize + 4 + ISIZE*glitmp
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
      end subroutine surf_mfo_outvtk

!---------------------------------------------------------------------- 

!     based on stat_mfo_write_hdr2D
!     write hdr, byte key, global ordering
      subroutine surf_mfo_write_hdr(isf)

      implicit none

      include 'SIZE_DEF'        ! missing definitions in include files
      include 'SIZE'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'SURF_STATS'           ! Slice speciffic variables


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

!     function      
      integer iglsum

      integer itmp,isf,glno

      itmp = iglsum(SURF_MEMCOUNT(isf),1)

#ifdef MPIIO
      nfileoo = 1   ! all data into one file
      nelo = iglsum(SURF_MEMCOUNT(isf),1)
#else
      nfileoo = NFILEO
      if(NID.eq.PID0) then                ! how many elements to dump
        nelo = SURF_MEMCOUNT(isf)
        do jl = PID0+1,PID1
           mtype = jl
           call csend(mtype,idum,ISIZE,jl,0)   ! handshake
           call crecv(mtype,inelp,ISIZE)
           nelo = nelo + inelp
        enddo
      else
        mtype = NID
        call crecv(mtype,idum,ISIZE)          ! hand-shake
        call csend(mtype,SURF_MEMCOUNT(isf),ISIZE,PID0,0)   ! u4 :=: u8
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
         write(RDCODE1(3),'(I1)') SURF_NFLDS/10
         write(RDCODE1(4),'(I1)') SURF_NFLDS-(SURF_NFLDS/10)*10

!     two last variables added for statistics
         write(hdr,1) WDSIZO,NXO,NYO,NZO,nelo,itmp,TIME,SURF_NSAMPLES,
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
     $     'Error writing header in surf_mfo_write_hdr. $')

!     write global 2D elements numbering for this group
!     copy data
      lglist(0) = SURF_MEMCOUNT(isf)
      kl = 0
      do il=1,SURF_MEMCOUNT(isf)
        kl = kl +1
        glno = lglel(SURF_OBJ(1,il,isf))
        lglist(kl) = glno
      enddo
!     check consistency
      ierr = 0
      if (kl.ne.SURF_MEMCOUNT(isf)) ierr=1
      call err_chk(ierr,'Error stat; inconsistent SURF_MEMCOUNT.3 $')

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

         len = ISIZE*(SURF_MEMCOUNT(isf)+1)
         call csend(mtype,lglist,len,PID0,0)  
      endif 

      call err_chk(ierr,
     $     'Error writing global nums in slice_mfo_write_hdr. $')

      return
      end subroutine surf_mfo_write_hdr

!-----------------------------------------------------------------------


