!=======================================================================
!     Adam Peplinski; 2015.11.30
!     Set of subroutines to add sponge inthe domain
!
!     Parameters used by this set of subroutines:
!     SPONGE:
!
!
!=======================================================================
!***********************************************************************
!     read parameters SPONGE
      subroutine spng_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'SPONGE'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr, len, il, ip, lsent
      parameter (lsent = (1+3*LDIM))
      real rtmp(lsent)

!     namelists
      namelist /SPONGE/ SPNG_STR, SPNG_W, SPNG_WL, SPNG_WR
!-----------------------------------------------------------------------
!     default values
      SPNG_STR = 1.0
      SPNG_W(1) = 5.0
      SPNG_WL(1) = 2.0
      SPNG_WR(1) = 2.0
      do il=2,NDIM
         SPNG_W(il) = 0.0
         SPNG_WL(il) = 0.0
         SPNG_WR(il) = 0.0
      enddo
!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=SPONGE,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading SPONGE parameters.$')

!     broadcast data
      if (NID.eq.0) then
         ip = 1
         rtmp(ip) = SPNG_STR
         do il=1,NDIM
            ip = ip +1
            rtmp(ip) = SPNG_W(il)
         enddo
         do il=1,NDIM
            ip = ip +1
            rtmp(ip) = SPNG_WL(il)
         enddo
         do il=1,NDIM
            ip = ip +1
            rtmp(ip) = SPNG_WR(il)
         enddo
      endif
      len = lsent*WDSIZE
      call bcast(rtmp,len)
      if (NID.ne.0) then
         ip = 1
         SPNG_STR = rtmp(ip)
         do il=1,NDIM
            ip = ip +1
            SPNG_W(il) = rtmp(ip)
         enddo
         do il=1,NDIM
            ip = ip +1
            SPNG_WL(il) = rtmp(ip)
         enddo
         do il=1,NDIM
            ip = ip +1
            SPNG_WR(il) = rtmp(ip)
         enddo
      endif

      return
      end
!***********************************************************************
!     write parameters SPONGE
      subroutine spng_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'SPONGE'          !

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /SPONGE/ SPNG_STR, SPNG_W, SPNG_WL, SPNG_WR
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=SPONGE,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing SPONGE parameters.$')

      return
      end
!***********************************************************************
!     calcualte spnge forcing
      subroutine spng_box_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D
      include 'PARALLEL_DEF'
      include 'PARALLEL'        ! GLLEL
      include 'SOLN_DEF'
      include 'SOLN'            ! JP
      include 'SPONGE'

!     argument list
      real ffx, ffy, ffz
      integer ix,iy,iz,ieg

!     local variables
      integer e, ip
      integer rtmp
!-----------------------------------------------------------------------
      if (SPNG_STR.gt.0.0) then
         e=GLLEL(ieg)
         ip=ix+NX1*(iy-1+NY1*(iz-1+NZ1*(e-1)))

!     dns
         if (JP.eq.0) then
            ffx = ffx + SPNG_FUN(ip)*(SPNG_VR(ip,1) - VX(ix,iy,iz,e))
            ffy = ffy + SPNG_FUN(ip)*(SPNG_VR(ip,2) - VY(ix,iy,iz,e))
            if (IF3D) ffz = ffz + SPNG_FUN(ip)*
     $           (SPNG_VR(ip,NDIM) - VZ(ix,iy,iz,e))
         else
            ffx = ffx - SPNG_FUN(ip)*VXP(ip,JP)
            ffy = ffy - SPNG_FUN(ip)*VYP(ip,JP)
            if(IF3D) ffz = ffz - SPNG_FUN(ip)*VZP(ip,JP)
         endif

      endif

      return
      end
!***********************************************************************
!     init sponge in the box domain
      subroutine spng_box_init(lvx,lvy,lvz)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'            ! [XYZ]M1
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D
      include 'SPONGE'          !

!     argument list
!     reference field
      real lvx(LX1*LY1*LZ1*LELV),lvy(LX1*LY1*LZ1*LELV),
     $     lvz(LX1*LY1*LZ1*LELV)

!     local variables
!     do this only once
      integer icalld
      save    icalld
      data    icalld  /0/

      integer ntot, il, jl
      real bmin(LDIM), bmax(LDIM)

      real rtmp, xxmax, xxmax_c, xxmin, xxmin_c, arg
      real lcoord(LX1*LY1*LZ1*LELV)


!     functions
      real glmin, glmax, spng_step
!-----------------------------------------------------------------------
      if (icalld.eq.0) then  ! do it only once
         icalld = icalld + 1

!     get box size
         ntot = NX1*NY1*NZ1*NELV
         bmin(1) = glmin(XM1,ntot)
         bmax(1) = glmax(XM1,ntot)
         bmin(2) = glmin(YM1,ntot)
         bmax(2) = glmax(YM1,ntot)
         if(IF3D) then
            bmin(NDIM) = glmin(ZM1,ntot)
            bmax(NDIM) = glmax(ZM1,ntot)
         endif

!     zero SPNG_FUN
         call rzero(SPNG_FUN,ntot)
            
         if(SPNG_STR.gt.0.0) then

!     stamp the file
            if (NIO.eq.0) then 
               write(6,*) 'SPONGE TURNED ON'
               write(6,*) 'sponge strenght = ' , SPNG_STR
               write(6,*) 'sponge width = ', SPNG_W
               write(6,*) 'sponge drop width = ', SPNG_WL
               write(6,*) 'sponge rise width = ', SPNG_WR
            endif

!     save reference field
            call copy(SPNG_VR(1,1),lvx, ntot)
            call copy(SPNG_VR(1,2),lvy, ntot)
            if (IF3D) call copy(SPNG_VR(1,NDIM),lvz, ntot)

!     for every dimension
            do il=1,NDIM

               if (SPNG_W(il).gt.0.0) then
                  if (SPNG_W(il).lt.(SPNG_WL(il)+SPNG_WR(il)).or.
     $                 SPNG_WL(il)+SPNG_WR(il).eq.0.0) then
                     if (NIO.eq.0) then
                        write(6,*) 'ERROR; wrong sponge parameters'
                     endif
                     call exitt
                  endif
     
                  rtmp = SPNG_W(il)/(SPNG_WL(il)+SPNG_WR(il))
!     sponge beginning (rise at xmax; right)
                  xxmax = bmax(il) - SPNG_WR(il)*rtmp
!     end (drop at xmin; left)
                  xxmin = bmin(il) + SPNG_WL(il)*rtmp
!     beginnign of constant part (right)
                  xxmax_c = xxmax + SPNG_WR(il)
!     beginnign of constant part (left)
                  xxmin_c = xxmin - SPNG_WL(il)

!     get SPNG_FUN
                  if (xxmax.le.xxmin) then
                     if (NIO.eq.0) write(6,*) 'ERROR; sponge to wide'
                     call exitt
                  else
!     this should be done by pointers, but for now I avoid it
                     if (il.eq.1) then
                        call copy(lcoord,XM1, ntot)
                     elseif (il.eq.2) then
                        call copy(lcoord,YM1, ntot)
                     elseif (il.eq.3) then
                        call copy(lcoord,ZM1, ntot)
                     endif

                     do jl=1,ntot
                        rtmp = lcoord(jl)
                        if(rtmp.le.xxmin_c) then ! constant; xmin
                           rtmp=SPNG_STR
                        elseif(rtmp.lt.xxmin) then ! fall; xmin
                           arg = (xxmin-rtmp)/SPNG_WL(il)
                           rtmp = spng_step(arg)
                        elseif (rtmp.le.xxmax) then ! zero
                           rtmp = 0.0
                        elseif (rtmp.lt.xxmax_c) then ! rise
                           arg = (rtmp-xxmax)/SPNG_WR(il)
                           rtmp = spng_step(arg)
                        else    ! constant
                           rtmp = SPNG_STR
                        endif
                        SPNG_FUN(jl)=max(SPNG_FUN(jl),rtmp)
                     enddo

                  endif         ! xxmax.le.xxmin

               endif            ! SPNG_W(il).gt.0.0
            enddo

         endif                  ! SPNG_STR.gt.0.0
      endif                     ! icalld.eq.0

!      call outpost2(SPNG_VR,SPNG_VR(1,2),SPNG_VR(1,3),SPNG_FUN,
!     &              SPNG_FUN,0,'   ')
!      call outpost2(SPNG_FUN,SPNG_VR(1,2),SPNG_VR(1,3),SPNG_FUN,
!     &              SPNG_FUN,0,'   ')


      return
      end
!***********************************************************************
      real function spng_step(x)
!
!     Smooth step function:
!     x<=0 : step(x) = 0
!     x>=1 : step(x) = 1
!     Non-continuous derivatives at x=0.02 and x=0.98
!
      implicit none
!     argument list
      real x
!     local variables
!     this should be done in a general way, but for now
      real eps
      parameter (eps=0.01)
!-----------------------------------------------------------------------
      if (x.le.eps) then
         spng_step = 0.0
      else
         if (x.le.(1.0d0-eps)) then
            spng_step = 1./( 1. + exp(1./(x - 1.) + 1./x) )
         else
            spng_step = 1.
         end if
      end if

      end function spng_step
!***********************************************************************
