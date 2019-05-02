!!======================================================================
!!
!!   Bunch of miscelineous routines.
!!   misc.f made by Prabal Negi
!!   Routines are not mine...
!!
!!======================================================================
 
      real function step_simson(x)
!     belongs to bla version 2.2
!     for more info see the bla.f file

      implicit none
     
      real x
      if(x.le.0.02) then
        step_simson=0.
      else
        if(x.le.0.98) then
           step_simson=1./(1.+exp(1./(x-1.)+1./x))
        else
           step_simson=1.
        endif
      endif

      return
      end
!-----------------------------------------------------------------------

      real function fringe_simson(x,xstart,xrise,xend,xfall,lambda_max)

!     Fringe function from Simson

      implicit none

      real step_simson        ! function
      real x                  ! position
      real xstart             ! start of fringe
      real xrise              ! rise distance at the start
      real xend               ! end of fringe
      real xfall              ! fall distance at the end
      real lambda_max         ! maximum fringe strength


      fringe_simson = lambda_max*(step_simson((x-xstart)/xrise) -   
     $                            step_simson((x-xend)/xfall + 1.) )

      return
      end     
!-----------------------------------------------------------------------
      subroutine set_obj  ! define objects for surface integrals
c
      include 'SIZE'
      include 'TOTAL'
!      include 'SURF_STATS'

      integer e,f,eg,isf

      nobj = 1
      iobj = 0
      do ii=nhis+1,nhis+nobj
         iobj = iobj+1
         hcode(10,ii) = 'I'
         hcode( 1,ii) = 'F'
         hcode( 2,ii) = 'F'
         hcode( 3,ii) = 'F'
         lochis(1,ii) = iobj
      enddo
      nhis = nhis + nobj

      if (maxobj.lt.nobj) call exitti('increase maxobj in SIZE$',nobj)

      nxyz  = nx1*ny1*nz1
      nface = 2*ndim

      do e=1,nelv
        do f=1,nface
          do isf = 1,1 !NSURFS 
!            if (cbc(f,e,1).eq.SURF_DEF(isf)) then
            if (cbc(f,e,1).eq.'mv ') then
              iobj  = isf
              if (iobj.gt.0) then
                 nmember(iobj) = nmember(iobj) + 1
                 mem = nmember(iobj)
                 eg  = lglel(e)
                 object(iobj,mem,1) = eg
                 object(iobj,mem,2) = f
!                write(6,1) iobj,mem,f,eg,e,nid,' OBJ'
!   1            format(6i9,a4)
              endif
            endif
          enddo
        enddo
      enddo

c     write(6,*) 'number',(nmember(k),k=1,4)
c
      return
      end
!-----------------------------------------------------------------------

      subroutine torque_calc_axis(scale,x0,ifdout,iftout,velx,vely,velz,
     $                            press)
c
c     Compute torque about point x0
c
c     Scale is a user-supplied multiplier so that results may be
c     scaled to any convenient non-dimensionalization.
c
c
      implicit none

      INCLUDE 'SIZE_DEF'
      INCLUDE 'SIZE'
      INCLUDE 'INPUT_DEF'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM_DEF'
      INCLUDE 'GEOM'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'
      INCLUDE 'SOLN_DEF'
      INCLUDE 'SOLN'          ! vdiff
      INCLUDE 'PARALLEL_DEF'
      INCLUDE 'PARALLEL'      ! gllnid

      real flow_rate,base_flow,domain_length,xsec,scale_vf
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                ,scale_vf(3)

      real scale
      real x0(3),w1(0:maxobj)
      logical ifdout,iftout

      real sij,pm1,xm0,ym0,zm0
      common /scrns/         sij (lx1*ly1*lz1*6*lelv)
      common /scrcg/         pm1 (lx1,ly1,lz1,lelv)
      common /scrsf/         xm0(lx1,ly1,lz1,lelt)
     $,                      ym0(lx1,ly1,lz1,lelt)
     $,                      zm0(lx1,ly1,lz1,lelt)

      real ur,us,ut,vr,vs,vt,wr,ws,wt
      integer lr
      parameter (lr=lx1*ly1*lz1)
      common /scruz/         ur(lr),us(lr),ut(lr)
     $                     , vr(lr),vs(lr),vt(lr)
     $                     , wr(lr),ws(lr),wt(lr)
c

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

      real xmx,xmn,ymx,ymn,zmx,zmn
      COMMON /XYZRANGE/ xmx,xmn,ymx,ymn,zmx,zmn

      integer lt,lt2
      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real velx(lt),vely(lt),velz(lt),press(lt2)

      real glmax,glmin        ! functions
      integer iglmax
      real dnekclock          ! function
      
      real torq_timer
      real sa                 ! local area
      real sarea(0:maxobj)    ! total area

      integer icalld
      save icalld
      data icalld /0/

      integer i0,ie,ieg,ifc,iobj,mem,memtot
      integer i,ii,n,nij


      torq_timer = dnekclock()      
c
      n = nx1*ny1*nz1*nelv
c
!      call mappr(pm1,pr,xm0,ym0) ! map pressure onto Mesh 1
      call mappr(pm1,press,xm0,ym0)       ! xm0,ym0 used as work arrays
c
c    Add mean_pressure_gradient.X to p:

      dpdx_mean=0.
      dpdy_mean=0.
      dpdz_mean=0.

      if (param(55).ne.0) then
        dpdx_mean = -scale_vf(1)
        dpdy_mean = -scale_vf(2)
        dpdz_mean = -scale_vf(3)
      endif

      call add2s2(pm1,xm1,dpdx_mean,n)  ! Doesn't work if object is cut by 
      call add2s2(pm1,ym1,dpdy_mean,n)  ! periodicboundary.  In this case,
      call add2s2(pm1,zm1,dpdz_mean,n)  ! set ._mean=0 and compensate in
c
c    Compute sij
c
      nij = 3
      if (if3d.or.ifaxis) nij=6
      call comp_sij(sij,nij,velx,vely,velz,ur,us,ut,vr,vs,vt,wr,ws,wt)
c
c
c     Fill up viscous array w/ default
c
      if (icalld.eq.0) then
        call cfill(vdiff,param(2),n)
        xmx = glmax(xm1,n)
        xmn = glmin(xm1,n)
        ymx = glmax(ym1,n)
        ymn = glmin(ym1,n)
        zmx = glmax(zm1,n)
        zmn = glmin(zm1,n)

        icalld = icalld+1
      endif

!     If point of torque calculation is located outside the domain
!     Then calculate about that axis instead of the point.
      if (x0(1).lt.xmn.or.x0(1).gt.xmx) then
          call rzero(xm0,n)
      else
          call cadd2(xm0,xm1,-x0(1),n)
      endif

      if (x0(2).lt.ymn.or.x0(2).gt.ymx) then
          call rzero(ym0,n)
      else
          call cadd2(ym0,ym1,-x0(2),n)
      endif

      if ((x0(3).lt.zmn.or.x0(3).gt.zmx).or.(ndim.eq.2)) then
          call rzero(zm0,n)
      else
          call cadd2(zm0,zm1,-x0(3),n)
      endif

!     Don't think we need these.      
!      x1min=glmin(xm0(1,1,1,1),n)
!      x2min=glmin(ym0(1,1,1,1),n)
!      x3min=glmin(zm0(1,1,1,1),n)
!
!      x1max=glmax(xm0(1,1,1,1),n)
!      x2max=glmax(ym0(1,1,1,1),n)
!      x3max=glmax(zm0(1,1,1,1),n)

!      x1min=xmn
!      x2min=ymn
!      x3min=zmn
!
!      x1max=xmx
!      x2max=ymx
!      x3max=zmx

      do i=0,maxobj
         dragpx(i) = 0   ! BIG CODE  :}
         dragvx(i) = 0
         dragx (i) = 0
         dragpy(i) = 0
         dragvy(i) = 0
         dragy (i) = 0
         dragpz(i) = 0
         dragvz(i) = 0
         dragz (i) = 0
         torqpx(i) = 0
         torqvx(i) = 0
         torqx (i) = 0
         torqpy(i) = 0
         torqvy(i) = 0
         torqy (i) = 0
         torqpz(i) = 0
         torqvz(i) = 0
         torqz (i) = 0
         sarea(i) = 0   
      enddo
c
c
      nobj = 0
      do ii=1,nhis
        if (hcode(10,ii).EQ.'I') then
          iobj   = lochis(1,ii)
          memtot = nmember(iobj)
          nobj   = max(iobj,nobj)
c
          if (hcode(1,ii).ne.' ' .or. hcode(2,ii).ne.' ' .or.
     $      hcode(3,ii).ne.' ' ) then
            ifield = 1
c
c           Compute drag for this object
c
            do mem=1,memtot
               ieg   = object(iobj,mem,1)
               ifc   = object(iobj,mem,2)
               if (gllnid(ieg).eq.nid) then ! this processor has a contribution
                  ie = gllel(ieg)
                  call clcdcm(dgtq,sa,xm0,ym0,zm0,sij,pm1,vdiff,ifc,ie)
c
                  call cmult(dgtq,scale,12)
c
                  dragpx(iobj) = dragpx(iobj) + dgtq(1,1)  ! pressure 
                  dragpy(iobj) = dragpy(iobj) + dgtq(2,1)
                  dragpz(iobj) = dragpz(iobj) + dgtq(3,1)
c
                  dragvx(iobj) = dragvx(iobj) + dgtq(1,2)  ! viscous
                  dragvy(iobj) = dragvy(iobj) + dgtq(2,2)
                  dragvz(iobj) = dragvz(iobj) + dgtq(3,2)
c
                  torqpx(iobj) = torqpx(iobj) + dgtq(1,3)  ! pressure 
                  torqpy(iobj) = torqpy(iobj) + dgtq(2,3)
                  torqpz(iobj) = torqpz(iobj) + dgtq(3,3)
c
                  torqvx(iobj) = torqvx(iobj) + dgtq(1,4)  ! viscous
                  torqvy(iobj) = torqvy(iobj) + dgtq(2,4)
                  torqvz(iobj) = torqvz(iobj) + dgtq(3,4)
                  sarea(iobj) = sarea(iobj) + sa
c
               endif
            enddo
          endif
        endif
      enddo

c
c     Sum contributions from all processors
c
      call gop(dragpx,w1,'+  ',maxobj+1)
      call gop(dragpy,w1,'+  ',maxobj+1)
      call gop(dragpz,w1,'+  ',maxobj+1)
      call gop(dragvx,w1,'+  ',maxobj+1)
      call gop(dragvy,w1,'+  ',maxobj+1)
      call gop(dragvz,w1,'+  ',maxobj+1)
c
      call gop(torqpx,w1,'+  ',maxobj+1)
      call gop(torqpy,w1,'+  ',maxobj+1)
      call gop(torqpz,w1,'+  ',maxobj+1)
      call gop(torqvx,w1,'+  ',maxobj+1)
      call gop(torqvy,w1,'+  ',maxobj+1)
      call gop(torqvz,w1,'+  ',maxobj+1)
      call gop(sarea,w1,'+  ',maxobj+1)

!      write(6,*) nid, 'Area :', sarea
c
      nobj = iglmax(nobj,1)
c
      do i=1,nobj
         dragx(i) = dragpx(i) + dragvx(i)
         dragy(i) = dragpy(i) + dragvy(i)
         dragz(i) = dragpz(i) + dragvz(i)
c
         torqx(i) = torqpx(i) + torqvx(i)
         torqy(i) = torqpy(i) + torqvy(i)
         torqz(i) = torqpz(i) + torqvz(i)
c
         dragpx(0) = dragpx (0) + dragpx (i)
         dragvx(0) = dragvx (0) + dragvx (i)
         dragx (0) = dragx  (0) + dragx  (i)
c
         dragpy(0) = dragpy (0) + dragpy (i)
         dragvy(0) = dragvy (0) + dragvy (i)
         dragy (0) = dragy  (0) + dragy  (i)
c
         dragpz(0) = dragpz (0) + dragpz (i)
         dragvz(0) = dragvz (0) + dragvz (i)
         dragz (0) = dragz  (0) + dragz  (i)
c
         torqpx(0) = torqpx (0) + torqpx (i)
         torqvx(0) = torqvx (0) + torqvx (i)
         torqx (0) = torqx  (0) + torqx  (i)
c
         torqpy(0) = torqpy (0) + torqpy (i)
         torqvy(0) = torqvy (0) + torqvy (i)
         torqy (0) = torqy  (0) + torqy  (i)
c
         torqpz(0) = torqpz (0) + torqpz (i)
         torqvz(0) = torqvz (0) + torqvz (i)
         torqz (0) = torqz  (0) + torqz  (i)
c
      enddo
c
      i0 = 0
      if (nobj.le.1) i0 = 1  ! one output for single-object case

      torq_timer = dnekclock() - torq_timer
c
      do i=i0,nobj
        if (nio.eq.0) then
!          write(6,*) 'Drag/Torque calculations'
          if (if3d.or.ifaxis) then
           if (ifdout) then
            write(6,6) istep,time,torq_timer,
     $              dragx(i),dragpx(i),dragvx(i),sarea(i),i,'dragx'
            write(6,6) istep,time,torq_timer,
     $              dragy(i),dragpy(i),dragvy(i),sarea(i),i,'dragy'
            write(6,6) istep,time,torq_timer,
     $              dragz(i),dragpz(i),dragvz(i),sarea(i),i,'dragz'
           endif
           if (iftout) then
            write(6,6) istep,time,torq_timer,
     $              torqx(i),torqpx(i),torqvx(i),sarea(i),i,'torqx'
            write(6,6) istep,time,torq_timer,
     $              torqy(i),torqpy(i),torqvy(i),sarea(i),i,'torqy'
            write(6,6) istep,time,torq_timer,
     $              torqz(i),torqpz(i),torqvz(i),sarea(i),i,'torqz'
           endif
          else
           if (ifdout) then
            write(6,6) istep,time,torq_timer,
     $              dragx(i),dragpx(i),dragvx(i),sarea(i),i,'dragx'
            write(6,6) istep,time,torq_timer,
     $              dragy(i),dragpy(i),dragvy(i),sarea(i),i,'dragy'
           endif
           if (iftout) then
            write(6,6) istep,time,torq_timer,
     $              torqz(i),torqpz(i),torqvz(i),sarea(i),i,'torqz'
           endif
          endif
        endif
    6   format(i8,1p6e16.8,1x,i3.1,a6)
      enddo
c
      return
      end subroutine torque_calc_axis
!-----------------------------------------------------------------------
      subroutine clcdcm(dgtq,a,xm0,ym0,zm0,sij,pm1,visc,f,e)
c
      implicit none

      INCLUDE 'SIZE_DEF'
      INCLUDE 'SIZE'
      INCLUDE 'GEOM_DEF'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT_DEF'
      INCLUDE 'INPUT'
      INCLUDE 'TOPOL_DEF'
      INCLUDE 'TOPOL'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'
c
      real dgtq(3,4)
      real xm0 (lx1,ly1,lz1,lelt)
      real ym0 (lx1,ly1,lz1,lelt)
      real zm0 (lx1,ly1,lz1,lelt)
      real sij (lx1,ly1,lz1,3*ldim-3,lelv)
      real pm1 (lx1,ly1,lz1,lelv)
      real visc(lx1,ly1,lz1,lelv)
c
      real dg(3,2)
c
      integer f,e
      real    n1,n2,n3

      real a                  ! total (local) area

      integer i,l,k
      integer j1,j2,js1,js2,jf1,jf2,jskip1,jskip2
      integer pf
      real s11,s21,s31,s12,s22,s32,s13,s23,s33
      real v
      real r1,r2,r3

c
      call dsset(nx1,ny1,nz1)    ! set up counters
      pf     = eface1(f)         ! convert from preproc. notation
      js1    = skpdat(1,pf)
      jf1    = skpdat(2,pf)
      jskip1 = skpdat(3,pf)
      js2    = skpdat(4,pf)
      jf2    = skpdat(5,pf)
      jskip2 = skpdat(6,pf)
C
      call rzero(dgtq,12)
c
      if (if3d.or.ifaxis) then
       i = 0
       a = 0
       do j2=js2,jf2,jskip2
       do j1=js1,jf1,jskip1
         i = i+1
         n1 = unx(i,1,f,e)*area(i,1,f,e)
         n2 = uny(i,1,f,e)*area(i,1,f,e)
         n3 = unz(i,1,f,e)*area(i,1,f,e)
         a  = a +          area(i,1,f,e)
c
         v  = visc(j1,j2,1,e)
c
         s11 = sij(j1,j2,1,1,e)
         s21 = sij(j1,j2,1,4,e)
         s31 = sij(j1,j2,1,6,e)
c
         s12 = sij(j1,j2,1,4,e)
         s22 = sij(j1,j2,1,2,e)
         s32 = sij(j1,j2,1,5,e)
c
         s13 = sij(j1,j2,1,6,e)
         s23 = sij(j1,j2,1,5,e)
         s33 = sij(j1,j2,1,3,e)
c
         dg(1,1) = pm1(j1,j2,1,e)*n1     ! pressure drag
         dg(2,1) = pm1(j1,j2,1,e)*n2
         dg(3,1) = pm1(j1,j2,1,e)*n3
c
         dg(1,2) = -v*(s11*n1 + s12*n2 + s13*n3) ! viscous drag
         dg(2,2) = -v*(s21*n1 + s22*n2 + s23*n3)
         dg(3,2) = -v*(s31*n1 + s32*n2 + s33*n3)
c
         r1 = xm0(j1,j2,1,e)
         r2 = ym0(j1,j2,1,e)
         r3 = zm0(j1,j2,1,e)
c
         do l=1,2
         do k=1,3
            dgtq(k,l) = dgtq(k,l) + dg(k,l)
         enddo
         enddo
c
         dgtq(1,3) = dgtq(1,3) + (r2*dg(3,1)-r3*dg(2,1)) ! pressure
         dgtq(2,3) = dgtq(2,3) + (r3*dg(1,1)-r1*dg(3,1)) ! torque
         dgtq(3,3) = dgtq(3,3) + (r1*dg(2,1)-r2*dg(1,1))
c
         dgtq(1,4) = dgtq(1,4) + (r2*dg(3,2)-r3*dg(2,2)) ! viscous
         dgtq(2,4) = dgtq(2,4) + (r3*dg(1,2)-r1*dg(3,2)) ! torque
         dgtq(3,4) = dgtq(3,4) + (r1*dg(2,2)-r2*dg(1,2))
       enddo
       enddo

      else ! 2D

       i = 0
       a = 0
       do j2=js2,jf2,jskip2
       do j1=js1,jf1,jskip1
         i = i+1
         n1 = unx(i,1,f,e)*area(i,1,f,e)
         n2 = uny(i,1,f,e)*area(i,1,f,e)
         a  = a +          area(i,1,f,e)
         v  = visc(j1,j2,1,e)

         s11 = sij(j1,j2,1,1,e)
         s12 = sij(j1,j2,1,3,e)
         s21 = sij(j1,j2,1,3,e)
         s22 = sij(j1,j2,1,2,e)

         dg(1,1) = pm1(j1,j2,1,e)*n1     ! pressure drag
         dg(2,1) = pm1(j1,j2,1,e)*n2
         dg(3,1) = 0

         dg(1,2) = -v*(s11*n1 + s12*n2) ! viscous drag
         dg(2,2) = -v*(s21*n1 + s22*n2)
         dg(3,2) = 0.

         r1 = xm0(j1,j2,1,e)
         r2 = ym0(j1,j2,1,e)
         r3 = 0.

         do l=1,2
         do k=1,3
            dgtq(k,l) = dgtq(k,l) + dg(k,l)
         enddo
         enddo

         dgtq(1,3) = 0! dgtq(1,3) + (r2*dg(3,1)-r3*dg(2,1)) ! pressure
         dgtq(2,3) = 0! dgtq(2,3) + (r3*dg(1,1)-r1*dg(3,1)) ! torque
         dgtq(3,3) = dgtq(3,3) + (r1*dg(2,1)-r2*dg(1,1))

         dgtq(1,4) = 0! dgtq(1,4) + (r2*dg(3,2)-r3*dg(2,2)) ! viscous
         dgtq(2,4) = 0! dgtq(2,4) + (r3*dg(1,2)-r1*dg(3,2)) ! torque
         dgtq(3,4) = dgtq(3,4) + (r1*dg(2,2)-r2*dg(1,2))
       enddo
       enddo
      endif

      return
      end

!-----------------------------------------------------------------------
!      subroutine get_spectra(coeff,f)
!     
!      implicit none
!
!      include 'SIZE_DEF'
!      include 'SIZE'
!
!      integer i
!      real coeff(lx1,ly1,lz1,lelt)
!      real f(lx1,ly1,lz1,lelt)
!
!      do i=1,nelt
!          call err_est_el_lget(coeff(1,1,1,i),f(1,1,1,i))
!      enddo
!
!      return
!      end subroutine get_spectra

c---------------------------------------------------------------------- 

      subroutine write_hessenberg(workl,ncv,ptr)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'

      real workl(1)
      integer ncv,ptr
      integer i,j,k,l
      integer iunit,ierror
      character outfmt*132
        

      call blank(outfmt,132)
      write(outfmt,'(A1,I3,A14)') '(', ncv,'(E16.07E3,1x))'

      call IO_freeid(iunit, ierror)
      if (nid.eq.0) then
        open (unit=iunit,file='arpack_hessenberg',
     $      action='write',status='unknown',iostat=ierror)
        do i=0,ncv-1
          j=ptr+i
          k=ptr+i+(ncv-1)*ncv
          write(iunit,outfmt) (workl(l), l=j,k,ncv)
        enddo
        flush(iunit)
        close(iunit)
      endif 

      return
      end subroutine write_hessenberg
!-----------------------------------------------------------------------









