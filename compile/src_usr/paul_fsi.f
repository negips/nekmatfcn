!====================================================================== 
!     FSI subroutine from Paul
!     Author: Paul Fischer
!          
!
!====================================================================== 

      subroutine paul_fsi_main

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /mybase/ bvx(lt,3),bpr(lt),f_lift_base
      common /mycomm/ v_cyl(0:10),cyl_mass,stiffness

      n =nx1*ny1*nz1*nelt
      n2=nx2*ny2*nz2*nelt

!      nio = -99    ! Turn off verbose logfile
!      if (istep.le.100.or.mod(istep,iostep).eq.0) nio=nid

      if (istep.eq.0) then
       rad       = 0.5
       cyl_vol   = pi*rad*rad
       cyl_rho   = param(33)
       cyl_mass  = cyl_vol*cyl_rho
       fn        = param(35)           ! From Arne's paper
       stiffness = cyl_mass*(2*pi*fn)**2
       call rzero(v_cyl(0),11)
       v_cyl(0)=param(34)
       call cyl_center_mass(y_del)
       if (nio.eq.0) write(6,2) cyl_mass,stiffness,y_del,' cyl mass'
   2   format(1p3e16.8,a9)

       if (nio.eq.0) write(6,3) 'istep','time','umx','y_del',
     $       'f_lift_ns','v_cyl(1)','f_lift_ns','f_spring'
   3   format(A5,1x,7(A14,1x),' umx')

       if (nio.eq.0) then
         write(6,'(A6,1x,A3,1x,11(A10,1x))') 'Perms:','ist','Fs','Fg',
     $  'Fk_ext','alpha','bd_etav','etav_s','etav','eta','rhs',
     $  'lhs'

         write(6,'(A6,1x,A5,1x,4(A13,1x))') 'Alpha:','istep',
     $             'alpha','eta','etav','etav_s'
       endif

      else
!       Done in userchk
!       call set_obj              ! Define object for lift calculations

        call ns_lift  (f_lift_ns) ! Lift due to NS update

        call base_lift ! base_lift  ! Lift from unit vertical shift

        rhs = f_lift_ns + cyl_mass
     $                * ( bd(2)*v_cyl(1)
     $                  + bd(3)*v_cyl(2)
     $                  + bd(4)*v_cyl(3)
     $                  - bd(1)*v_cyl(0) ) / dt

        call cyl_center_mass(y_del)
        f_spring = -stiffness*y_del
        rhs = rhs + f_spring

        den = cyl_mass*bd(1)/dt - f_lift_base  ! Linear in base lift if m=0
        v_del = rhs / den                      ! Denominator should never be 0

        call add2s2(vx,bvx(1,1),v_del,n)   ! Update velocity field with
        call add2s2(vy,bvx(1,2),v_del,n)   ! implicit contribution
        call add2s2(vz,bvx(1,3),v_del,n)
        call add2s2(pr,bpr     ,v_del,n2)

        v_cyl(3)=v_cyl(2)
        v_cyl(2)=v_cyl(1)
        v_cyl(1)=v_cyl(0) + v_del  ! Preset, plus delta

        v_cyl(0)=(bd(2)*v_cyl(1)       ! This is just a projected preset
     $           +bd(3)*v_cyl(2)       ! estimate of Omega for next
     $           +bd(4)*v_cyl(3))/bd(1)! timestep that will generally
      endif                           ! give zero acceleration.  It is
                                      ! an arbitrary choice because we
                                      ! subsequently correct with delta.

      umx = glmax (vx,n)
      vmx = glamax(vy,n)

      if (istep.gt.0) then
       if (nio.eq.0) write(6,1) istep,time,umx,y_del,f_lift_ns,v_cyl(1)
     $              ,f_lift_ns,f_spring
    1  format(i5,7(e14.6,1x),' umx')

       if (nio.eq.0) write(6,4) 'Perms:',istep,f_lift_ns,
     $        f_lift_base,f_spring,v_del,0.,v_cyl(0),v_cyl(1),y_del,
     $        rhs,den
    4   format(A6,1x,I5,1x,11(E10.3E2,1x))
!        call outpost2(bvx(1,1),bvx(1,2),bvx(1,3),bpr,bvx,1,'pfs')
       if (nio.eq.0) write(6,'(A6,1x,I5,1x,4(E13.6E2,1x))') 
     $  'Alpha:',istep,v_del,y_del,v_cyl(1),v_cyl(0)

      endif


      ifusermv = .true.
!      if (ifusermv) call my_mv_mesh    ! Compute our own mesh velocity
      if (ifusermv) call my_meshv_eta(v_cyl(0))

      return
      end

!---------------------------------------------------------------------- 
      subroutine ns_lift(f_lift_ns)
      include 'SIZE'
      include 'TOTAL'

      common /ctorq/ dragx(0:maxobj),dragpx(0:maxobj),dragvx(0:maxobj)
     $             , dragy(0:maxobj),dragpy(0:maxobj),dragvy(0:maxobj)
     $             , dragz(0:maxobj),dragpz(0:maxobj),dragvz(0:maxobj)

     $             , torqx(0:maxobj),torqpx(0:maxobj),torqvx(0:maxobj)
     $             , torqy(0:maxobj),torqpy(0:maxobj),torqvy(0:maxobj)
     $             , torqz(0:maxobj),torqpz(0:maxobj),torqvz(0:maxobj)

     $             , dpdx_mean,dpdy_mean,dpdz_mean
     $             , dgtq(3,4)


      common /mycntr/ x0(3)

      call rzero(x0,3)      ! Center of cylinder, for torque_calc

      scale = 2.  ! Cd = F/(.5 rho U^2 ) = 2*F
      call torque_calc(scale,x0,.true.,.false.)

      f_lift_ns = dragy(0)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_my_vel_bc(vxc,vyc,vzc)

      include 'SIZE'
      include 'TOTAL'

      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)

      integer e,f

      call opzero(vxc,vyc,vzc)

      nface = 2*ndim
      do e=1,nelv
      do f=1,nface
c        if (cbc(f,e,1).eq.'v  '.and.abs(xm1(1,1,1,e)).lt.1.0) then
         if (cbc(f,e,1).eq.'mv '.and.abs(xm1(1,1,1,e)).lt.1.0) then
c        write(6,*) e,f,' ',cbc(f,e,1),xm1(1,1,1,e),' cb2'

           l=0
           call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)
           do k=k0,k1
           do j=j0,j1
           do i=i0,i1
              l=l+1
              x = xm1(i,j,k,e)
              y = ym1(i,j,k,e)
              z = zm1(i,j,k,e)
              vxc(i,j,k,e) =  0.
              vyc(i,j,k,e) =  1.
              vzc(i,j,k,e) =  0.
c             write(6,*) i,j,e,f,vxc(i,j,k,e),' vxc'
           enddo
           enddo
           enddo
         endif
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine plan3_base_flow(vxc,vyc,vzc,prc)

c     Compute PnPn-2 pressure and velocity using fractional step method.

      include 'SIZE'
      include 'TOTAL'
c
      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)
C
      common /scrns/ rw1   (lx1,ly1,lz1,lelv)
     $ ,             rw2   (lx1,ly1,lz1,lelv)
     $ ,             rw3   (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             respr (lx2,ly2,lz2,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)
      common /cvflow_i/ icvflow,iavflow


c     Compute velocity, 1st part 

      ntot1  = nx1*ny1*nz1*nelv
      ntot2  = nx2*ny2*nz2*nelv
      ifield = 1

      call set_my_vel_bc(vxc,vyc,vzc)
c     call outpost (vxc,vyc,vzc,prc,t,'   ')
c     write(6,*) dt,bd(1),' this is bd1'
c     stop


      intype = -1
      call sethlm   (h1,h2,intype)

      call ophx     (rw1,rw2,rw3,vxc,vyc,vzc,h1,h2)
      call opchsgn  (rw1,rw2,rw3)  ! rhs = -H*U_bdry

      call ophinv   (dv1,dv2,dv3,rw1,rw2,rw3,h1,h2,tolhv,nmxh)
      call opadd2   (vxc,vyc,vzc,dv1,dv2,dv3)

c     Compute pressure  (from "incompr")

      intype = 1
      dtinv  = 1./dt

      ifield = 1
      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call cmult   (h2,dtinv,ntot1)
      call invers2 (h2inv,h2,ntot1)
      call opdiv   (respr,vxc,vyc,vzc)
      call chsign  (respr,ntot2)
      call ortho   (respr)


c     Set istep=0 so that h1/h2 will be re-initialized in eprec
      i_tmp = istep
      istep = 0
      call esolver (respr,h1,h2,h2inv,intype)
      istep = i_tmp

      call opgradt (rw1,rw2,rw3,respr)
      call opbinv  (dv1,dv2,dv3,rw1,rw2,rw3,h2inv)
      call opadd2  (vxc,vyc,vzc,dv1,dv2,dv3)

      call cmult2  (prc,respr,bd(1),ntot2)

c     call outpost (vxc,vyc,vzc,prc,t,'   ')
c     write(6,*) dt,bd(1),' this is bd1'
c     stop
c     call exitti('quit base flow$',istep)

      return
      end
c-----------------------------------------------------------------------
      subroutine base_lift ! Lift and base flow from unit rotation rate
      include 'SIZE'
      include 'TOTAL'

      real dtbd_last
      save dtbd_last
      data dtbd_last / -9.9 /

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrmu/ svx(lt,3),spr(lt)
      common /mybase/ bvx(lt,3),bpr(lt),f_lift_base

      logical if_done

      dtbd    = bd(1)/dt

      if_done = .true.
      if (dtbd.ne.dtbd_last) if_done = .false.
      if (ifmvbd)            if_done = .false.

      if (if_done)           return

      dtbd_last = dtbd

      n  = nx1*ny1*nz1*nelt
      n2 = nx2*ny2*nz2*nelt
      call opcopy(svx(1,1),svx(1,2),svx(1,3),vx,vy,vz) ! Save current flow
      call copy  (spr,pr,n2)                           ! and pressure

      call plan3_base_flow(vx,vy,vz,pr)                ! Compute base flow
      call ns_lift(f_lift_base)                        ! Compute base lift

      call opcopy(bvx(1,1),bvx(1,2),bvx(1,3),vx,vy,vz) ! Save base flow
      call copy  (bpr,pr,n2)

      call opcopy(vx,vy,vz,svx(1,1),svx(1,2),svx(1,3)) ! Restore current
      call copy  (pr,spr,n2)                           ! flow and pressure

      return
      end
c-----------------------------------------------------------------------
      subroutine cyl_center_mass(y_bar)
      include 'SIZE'
      include 'TOTAL'
      common /ctmp0/ one(lx1*ly1*lz1*lelt)
      integer e,f,eg

      n=nx1*ny1*nz1*nelv
      call rone(one,n)

      y_bar = 0.
      y_sum = 0.
      a_sum = 0.

      nface = 2*ndim

      do e=1,nelv
      do f=1,nface
         if (cbc(f,e,1).eq.'mv '.and.abs(xm1(1,1,1,e)).lt.1.0) then
            y_sum = y_sum + facint_v(ym1,area,f,e)
            a_sum = a_sum + facint_v(one,area,f,e)
         endif
      enddo
      enddo

      y_sum=glsum(y_sum,1)
      a_sum=glsum(a_sum,1)
      if (a_sum.gt.0) y_bar = y_sum/a_sum

      return
      end
c-----------------------------------------------------------------------


