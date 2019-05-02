!    Template to extract wall values.
!---------------------------------------------------------------------- 
      subroutine extract_wall()  
c     Set the GLL points as defined in the CASENAME.grid file

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
!      include 'TOTAL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'PARALLEL_DEF'
      include 'PARALLEL'      
      include 'mpif.h'
      include 'TOPOL_DEF'
      include 'TOPOL'
      include 'SOLN_DEF'
      include 'SOLN'

      real dum(3)

      CHARACTER CB*3
      integer iel,counter,i,cnt_rem,j
      integer ntot1,ntott

      integer wallpars
      parameter (wallpars=10)
      real wall(lx1,ly1,lelt,wallpars),wall_tmp(lx1*ly1*lelt,wallpars)
      common /wallval/ wall

      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer iface,nfaces
      integer ix,iy,iz
      integer ifield
      integer gllcnt,elcnt

      real iobj

      character outfmt*32

      integer ierr,length,ip
      
      integer idir
      real ARC_LEN1D(LX1,LY1,LZ1,LELT)
      COMMON /ARCLEN/ ARC_LEN1D

      integer ppf
      real snx,sny,snz

      call blank(outfmt,32)
      write (outfmt, "(A1,I2,A13)") '(', wallpars, '(E15.8E2,2X))'
!      if (nid.eq.0) write(*,*) 'outfmt', outfmt

      ntott =  nx1*ny1*nz1*nelt
      ntot1 = nx1*ny1*nz1*nelv

      NFACES=2*NDIM

      ifxyo = .true.
      
      IFIELD = 1
      counter = 0
      elcnt = 0
      do  iel=1,NELV            !do ieg=1,nelgt
        idir=1
        call calc_arc_len(idir)
c       Read GLL points (2D) from file
        do IFACE = 1,NFACES
           CB = CBC(IFACE,iel,IFIELD)
           if (CB.EQ.'W  '.or.CB.eq.'mv ') then
              elcnt = elcnt + 1   
              CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $             ,IFACE)
              gllcnt = 0   
              do 202 iz=KZ1,KZ2
              do 202 iy=KY1,KY2
              do 202 ix=KX1,KX2
                counter = counter +1
                gllcnt = gllcnt+1
                wall(counter,1,1,1) = xm1(ix,iy,iz,iel)
                wall(counter,1,1,2) = ym1(ix,iy,iz,iel)
                wall(counter,1,1,3) = ARC_LEN1D(ix,iy,iz,iel)

                ppf = eface1(iface)  
              
                if (ppf.eq.1.or.ppf.eq.2) then      ! "r face"  
                   snx = unx(iy,iz,iface,iel)                 ! Note:  iy,iz  
                   sny = uny(iy,iz,iface,iel)  
                   snz = unz(iy,iz,iface,iel)  
                elseif (ppf.eq.3.or.ppf.eq.4)  then ! "s face"  
                   snx = unx(ix,iz,iface,iel)                 !        ix,iz  
                   sny = uny(ix,iz,iface,iel)  
                   snz = unz(ix,iz,iface,iel)  
                elseif (ppf.eq.5.or.ppf.eq.6)  then ! "t face"  
                   snx = unx(ix,iy,iface,iel)                 !        ix,iy  
                   sny = uny(ix,iy,iface,iel)  
                   snz = unz(ix,iy,iface,iel)  
                endif

                wall(counter,1,1,4) = snx
                wall(counter,1,1,5) = sny
                wall(counter,1,1,6) = snz
                wall(counter,1,1,7) = lglel(iel)+0.
                wall(counter,1,1,8) = iface+0.
                wall(counter,1,1,9) = v1mask(ix,iy,iz,iel)
                wall(counter,1,1,10)= v2mask(ix,iy,iz,iel)

!                if (IFACE.eq.2) then
!                      iobj = 1.0  ! upper wall
!                elseif (IFACE.eq.4) then 
!                      iobj = -1.0  ! lower wall
!                else
!                      iobj = 0.0
!                endif
!                wall(4,counter)    =  iobj
 202          continue   
           endif        ! cb.eq.'W  '
         enddo          ! iface
      enddo             ! iel

      if (nid.eq.0) then
         open(unit=74, FILE='surf_data.dat')
         write(*,*) 'Proc/NPTS',nid, counter
         write(74,'(10(A15,2x))') 'x','y','dx','snx','sny','snz','iel',
     $            'face','v1mask','v2mask'
         do i = 1,counter
            write(74,outfmt)  (wall(i,1,1,j), j = 1,wallpars)
         end do
         
         do ip=1,np-1
!     hand sahking
            length = 1*isize
            ierr = ip
            call csend(ip,ierr,length,ip,0)
!     count number
            call crecv(ip,cnt_rem,length)
            write(6,*) nid, 'Receiving: ', ip,cnt_rem
!     get data
            if (cnt_rem.gt.0) then
               length = cnt_rem*wdsize*wallpars*2
               call crecv(ip,wall_tmp,length)

               do i=1,wallpars
                  call copy(wall(1,1,1,i),wall_tmp((i-1)*cnt_rem+1,1),
     $                                                    cnt_rem)
               enddo

               do i = 1,cnt_rem
                  write(74,outfmt)  (wall(i,1,1,j),j = 1,wallpars)
               end do
            endif
            counter = counter+cnt_rem
         enddo
         close(74)

      else
!     hand shaking
         length = 1*isize
         call crecv(nid,ierr,length)
!     number of points
         call csend(nid,counter,length,0,0)

!     send data
         if(counter.gt.0) then
            do i=1,wallpars
               call copy(wall_tmp(counter*(i-1)+1,1),wall(1,1,1,i),
     $                                                 counter)
            enddo
            length = counter*wdsize*wallpars*2
            call csend(nid,wall_tmp,length,0,0)
         endif

      endif
      call nekgsync()

      if(nid.eq.0)then
         write(*,*)  'Total PTS:', counter
      end if

      return
      end

!----------------------------------------------------------------------

      subroutine calc_arc_len(arc_dir)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'MASS_DEF'
      include 'MASS'
      include 'WZ_DEF'
      include 'WZ'
      include 'DXYZ_DEF'
      include 'DXYZ'

!     scratch space
      real lxyzd(LX1,LY1,LZ1,LELT,3)
      common /SCRSF/ lxyzd      ! coordinate derivatives

!     local variables
      integer i, j, k, e        ! loop index
      integer el                ! index of 2D element
      real lwm1(LX1)       ! wieghts for 1D integration

      real ARC_LEN1D(LX1,LY1,LZ1,LELT)
      COMMON /ARCLEN/ ARC_LEN1D

      integer arc_dir

!      arc_dir = 1

!     copy wieghts depending on the uniform direction
      if (arc_dir.eq.1) then
         call copy(lwm1,WXM1,NX1)
!     get coordinates derivatives d[XYZ]/dr
         i = NY1*NZ1
         do e = 1, NELT
            call mxm(DXM1,NX1,XM1(1,1,1,e),NX1,lxyzd(1,1,1,e,1),i)
            call mxm(DXM1,NX1,YM1(1,1,1,e),NX1,lxyzd(1,1,1,e,2),i)
            call mxm(DXM1,NX1,ZM1(1,1,1,e),NX1,lxyzd(1,1,1,e,3),i)
         enddo
      elseif (arc_dir.eq.2) then
         call copy(lwm1,WYM1,NY1)
!     get coordinates derivatives d[XYZ]/ds
         do e = 1, NELT
           do i=1, NZ1
             call mxm(XM1(1,1,i,e),NX1,DYTM1,NY1,
     $            lxyzd(1,1,i,e,1),NY1)
             call mxm(YM1(1,1,i,e),NX1,DYTM1,NY1,
     $            lxyzd(1,1,i,e,2),NY1)
             call mxm(ZM1(1,1,i,e),NX1,DYTM1,NY1,
     $            lxyzd(1,1,i,e,3),NY1)
           enddo
         enddo
      else
         call copy(lwm1,WZM1,NZ1)
!     get coordinates derivatives d[XYZ]/dt
         i = NX1*NY1
         do e = 1, NELT
            call mxm(XM1(1,1,1,e),i,DZTM1,NZ1,lxyzd(1,1,1,e,1),NZ1)
            call mxm(YM1(1,1,1,e),i,DZTM1,NZ1,lxyzd(1,1,1,e,2),NZ1)
            call mxm(ZM1(1,1,1,e),i,DZTM1,NZ1,lxyzd(1,1,1,e,3),NZ1)
         enddo

      endif

!     get 1D mass matrix ordering directions in such a way that 
!     the uniform direction corresponds to the the first index
      i = NX1*NY1*NZ1
!     get arc length
      do e = 1, NELT
        call vsq(lxyzd(1,1,1,e,1),i)
        call vsq(lxyzd(1,1,1,e,2),i)
        call vsq(lxyzd(1,1,1,e,3),i)
      
        call add2(lxyzd(1,1,1,e,1),lxyzd(1,1,1,e,2),i)
        call add2(lxyzd(1,1,1,e,1),lxyzd(1,1,1,e,3),i)

        call vsqrt(lxyzd(1,1,1,e,1),i)
      enddo

      i=i*NELT
      call rzero(ARC_LEN1D,i)
      call copy (ARC_LEN1D,lxyzd,i)

      do e=1, NELT
        do k=1,NZ1
           do j=1, NY1
              do i=1, NX1
                 ARC_LEN1D(i,j,k,e) = lwm1(i)*ARC_LEN1D(i,j,k,e)
              enddo
           enddo
        enddo
      enddo


      return
      end subroutine calc_arc_len

!----------------------------------------------------------------------       
      
