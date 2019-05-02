!=======================================================================
! Name        : statistics_2D
! Author      : Prabal Singh Negi, Adam Peplinski
! Version     : last modification 2015.05.22
! Copyright   : GPL
! Description : This is a set of routines to calculate 2D statistics
!=======================================================================
c----------------------------------------------------------------------
!     read parameters Statistics 
      subroutine t_stat_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'STATS'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /TSTATS/ stat_comp,stat_outp,stat_idir,
     $                 stat_ifstpts,statpts_comp

!     default values
      stat_comp     = 10            ! compute interval 
      stat_outp     = 1.00E+4       ! saving interval
      stat_idir     = 3             ! Integration direction. 3=>Z    
      stat_ifstpts  = .false.       ! turn off time series
      statpts_comp  = 100           ! Time history saving interval    

!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=STATS,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading STATS parameters.$')

!     broadcast data
      call bcast(stat_comp,         ISIZE)
      call bcast(stat_outp,         ISIZE)
      call bcast(stat_idir,         ISIZE)
      call bcast(stat_ifstpts,      LSIZE)
      call bcast(statpts_comp,      ISIZE)

      return
      end
!-----------------------------------------------------------------------
!     write parameters Statistics 
      subroutine t_stat_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'STATS'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /STATS/ stat_comp,stat_outp,stat_idir,
     $                 stat_ifstpts,statpts_comp

!     read the file
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=STATS,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing STATS parameters.$')

      return
      end
!-----------------------------------------------------------------------
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!     Driver for statistics computation
      subroutine T_STAT_AVG_ALL

      implicit none 

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TSTATS'

      integer COMPS
      save COMPS
      data COMPS /0/


!     Time dependent stats calculation
      if (mod(ISTEP,T_STAT_COMP).LT.T_STAT_AVGITER) then
        if (COMPS.eq.0) then
          T_STAT_TSTART=TIME-DT
          dtime=TIME-T_STAT_TSTART
          T_STAT_ATIME=TIME-T_STAT_TSTART
          tbeta=dtime/T_STAT_ATIME
          talpha=1.0 - tbeta
        else
          dtime=TIME-T_STAT_ATIME-T_STAT_TSTART
          T_STAT_ATIME=TIME-T_STAT_TSTART
          tbeta=dtime/T_STAT_ATIME
          talpha=1.0 - tbeta
        endif
 
        COMPS=COMPS+1
        call T_STAT_COMPUTE(talpha,tbeta,comps)

        if (COMPS.eq.T_STAT_AVGITER) then
          call T_STATS_SAVE
          COMPS=0
        endif
      endif

      if (T_STAT_NSAVES.eq.T_STAT_MAXTSAVES) then
        call T_STAT_COLLATE
        call T_STATS_OUT
!        call exitt                ! prabal  
        continue
      endif

      if (ISTEP.eq.NSTEPS) call stat_end        ! finalize

      return
      end subroutine T_STAT_AVG_ALL

!====================================================================== 


