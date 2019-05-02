!======================================================================
!     Routines for fluid structure interaction implicit solutions
!     using fixed-point method with dynamic relaxation.
!     Subroutines for adjoint FSI calculation      
!     Author: Prabal S. Negi
!======================================================================       
      subroutine nlfsi_main_adj

      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'MVGEOM_DEF'
      include 'MVGEOM'
      include 'NLFSI'
      include 'CHKPOINT'
      include 'ADJOINT_DEF'
      include 'ADJOINT'

      real scale
      real x0(3)
      logical ifdout,iftout
      integer ierr

      integer icalld
      save icalld
      data icalld /0/

!     Takes from Calc_torque_axis/pert 
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

      integer lvp,lvp2
      parameter(lvp=lpx1*lpy1*lpz1*lpelv)
      parameter(lvp2=lpx2*lpy2*lpz2*lpelv)
      real wk_vxp(lvp),wk_vyp(lvp),wk_vzp(lvp),wk_prp(lvp2)

      real wk_vx0(lvp),wk_vy0(lvp),wk_vz0(lvp),wk_pr0(lvp2)
      real tmp_pr(lvp)

      real rot_sn(3,3)
      real rot_cs(3,3)
      real dxyz(3),r0(3)
      real psi_sx,psi_sy,psi_sz

      real*8 dnekclock          ! function. timing
      real*8 steptime           ! timing for step

      real rxyz 
      common /scrsf/ rxyz(lx1*ly1*lz1*lelt,3)

      integer i


      if (icalld.eq.0.and..not.nlfsi_ifinit) then
!       Initialize timer 
        nlfsi_timea = 0.
        nlfsi_timei = 0.
           
        call nlfsi_init
        
        icalld=icalld+1
      endif  

      if (.not.IFNLFSI) then
        psiv_s = 0.
        psiv_m = 0.
        psiv   = 0.
        psi    = 0.
        psi_s  = 0.
        psia_s = 0.
        psia   = 0.    
        return
      endif  

!     Start timing for this step   
      steptime = dnekclock()

      if (ifchkptrst.and.istep.lt.chkptnrsf) then
!       Read saved values for restart
!       If the simulation has been restarted
        call nlfsi_rstread

        call nlfsi_output

        if (ifusermv.and.ifpert) then
          if (ifadj) then          
            call nlfsi_meshv_adj
          else
            call nlfsi_meshv_pert
          endif
        endif

        nlfsi_timei=dnekclock()-steptime
        nlfsi_timea=nlfsi_timea+nlfsi_timei
        return

      elseif (istep.eq.0) then
!       If there are no restarts              
!       At first time step the projected velocity is
!       just the initialized velocity
        psiv_s = psiv_ini
        psiv_m = 0.
        psiv   = psiv_ini
        psia   = 0.
        psia_s = 0.    
        psi    = psi_ini
        psi_s  = psi

        call nlfsi_output

        if (ifusermv.and.ifpert) then
          if (ifadj) then          
            call nlfsi_meshv_adj
          else
            call nlfsi_meshv_pert
          endif
        endif

        nlfsi_timei=dnekclock()-steptime
        nlfsi_timea=nlfsi_timea+nlfsi_timei
        return
      endif

      x0(1) = nlfsi_x0
      x0(2) = nlfsi_y0
      if (ndim.eq.3) x0(3) = 1000.         ! set outside the domain

      ifdout=.false.
      iftout=.false.

      scale=nlfsi_rescale
      
      call torque_calc_axis(scale,x0,ifdout,iftout,
     $             vxp(1,jp),vyp(1,jp),vzp(1,jp),prp(1,jp))

!     Who came up with this code? 
      call copy(nlfs_dragx (0),dragx (0),maxobj+1)
      call copy(nlfs_dragpx(0),dragpx(0),maxobj+1)
      call copy(nlfs_dragvx(0),dragpx(0),maxobj+1)
      call copy(nlfs_dragy (0),dragy (0),maxobj+1)
      call copy(nlfs_dragpy(0),dragpy(0),maxobj+1)
      call copy(nlfs_dragvy(0),dragvy(0),maxobj+1)
      call copy(nlfs_dragz (0),dragz (0),maxobj+1)
      call copy(nlfs_dragpz(0),dragpz(0),maxobj+1)
      call copy(nlfs_dragvz(0),dragvz(0),maxobj+1)
      call copy(nlfs_torqx (0),torqx (0),maxobj+1)
      call copy(nlfs_torqpx(0),torqpx(0),maxobj+1)
      call copy(nlfs_torqvx(0),torqpx(0),maxobj+1)
      call copy(nlfs_torqy (0),torqy (0),maxobj+1)
      call copy(nlfs_torqpy(0),torqpy(0),maxobj+1)
      call copy(nlfs_torqvy(0),torqvy(0),maxobj+1)
      call copy(nlfs_torqz (0),torqz (0),maxobj+1)
      call copy(nlfs_torqpz(0),torqpz(0),maxobj+1)
      call copy(nlfs_torqvz(0),torqvz(0),maxobj+1)

      
!     Adjoint forces for psi^+ equation
      call torque_calc_adj(scale,x0,ifdout,iftout,
     $         vxp(1,jp),vyp(1,jp),vzp(1,jp),prp(1,jp),
     $         nlfsi_grvx0,nlfsi_grvy0,nlfsi_grvz0)

      call copy(nlfadj_dragx (0),dragx (0),maxobj+1)
      call copy(nlfadj_dragpx(0),dragpx(0),maxobj+1)
      call copy(nlfadj_dragvx(0),dragpx(0),maxobj+1)
      call copy(nlfadj_dragy (0),dragy (0),maxobj+1)
      call copy(nlfadj_dragpy(0),dragpy(0),maxobj+1)
      call copy(nlfadj_dragvy(0),dragvy(0),maxobj+1)
      call copy(nlfadj_dragz (0),dragz (0),maxobj+1)
      call copy(nlfadj_dragpz(0),dragpz(0),maxobj+1)
      call copy(nlfadj_dragvz(0),dragvz(0),maxobj+1)
      call copy(nlfadj_torqx (0),torqx (0),maxobj+1)
      call copy(nlfadj_torqpx(0),torqpx(0),maxobj+1)
      call copy(nlfadj_torqvx(0),torqpx(0),maxobj+1)
      call copy(nlfadj_torqy (0),torqy (0),maxobj+1)
      call copy(nlfadj_torqpy(0),torqpy(0),maxobj+1)
      call copy(nlfadj_torqvy(0),torqvy(0),maxobj+1)
      call copy(nlfadj_torqz (0),torqz (0),maxobj+1)
      call copy(nlfadj_torqpz(0),torqpz(0),maxobj+1)
      call copy(nlfadj_torqvz(0),torqvz(0),maxobj+1)

      call nlfsi_struc_solve_adj

      nlfsi_timei=dnekclock()-steptime
      nlfsi_timea=nlfsi_timea+nlfsi_timei

      return
      end subroutine nlfsi_main_adj 
!---------------------------------------------------------------------- 
      subroutine nlfsi_struc_solve_adj

!     Calculate superposition of Solutions to get implicit update
!     In this routine we start with an assumed or extrapolated value of
!     psi (psi_s). Then solve the structural equations to obtain the new
!     value to psi and psiv.
!     Convergence is achieved if the new and extrapolated values are
!     consistent within a given tolerance.            

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'GEOM_DEF'      
      include 'GEOM'
      include 'NLFSI'

      integer icalld
      save icalld
      data icalld /0/

      real omegamin                 ! minimum omega
      parameter (omegamin=1.0e-10)  ! can't be negative

      real ext_k,bd_psi,bd_psiv,bd_psia
      real ab0,ab1,ab2
      integer iobj
      integer ilag
      real const
      real solid_area
      real fb

!      real alpha
      real Fs,Fd,Fk,Ft,Fadj

      integer istep_old

      integer it
      real rnorm,dr(2),dr_1(2)
      real omg,tmp1,tmp2
      logical ifbisection

      real inert,stiff,damp
      real lhs(2,2)          ! two-by-two matrix
      real inv_lhs(2,2)      ! inverse of the two-by-two matrix
      real rhs(2),psi_soln(2)
      real det, invdet       ! determinant and determinant inverse

      inert=nlfsi_inertia
      stiff=nlfsi_stiff
      damp=nlfsi_damp

      it = nlfsi_iter

!     Only doing for one object right now 
      iobj=0
      if (nlfsi_ifrot) then
!       Rotational coordinate system is opposite to the one used in Nek 
        Fs   = -nlfs_torqz(iobj)
!       Prabal. Must figure this out
        Fadj = -nlfadj_torqz(iobj) 
      else  
        Fs   = nlfs_dragy(iobj)
        Fadj = nlfadj_dragy(iobj)
      endif

!     Backward differentiation for psi
      bd_psi = psi*bd(2)
      do ilag=2,NBD
        bd_psi = bd_psi + bd(ilag+1)*psilag(ilag-1)
      enddo

!     Backward differentiation for psiv
      bd_psiv = psiv*bd(2)
      do ilag=2,NBD
        bd_psiv = bd_psiv + bd(ilag+1)*psivlag(ilag-1)
      enddo
      bd_psiv = bd_psiv

!     Backward differentiation for psia
      bd_psia = psia*bd(2)
      do ilag=2,NBD
        bd_psia = bd_psia + bd(ilag+1)*psialag(ilag-1)
      enddo

!     Nek makes a change of variable tau = -t
!     Therefore -d/dt --> d/dtau
!     Hence the negative signs from the time derivative are absorbed.      
      lhs(1,1) = bd(1)/DT*inert + damp
      lhs(1,2) = -1.
      lhs(2,1) = stiff
      lhs(2,2) = bd(1)/DT

      det    = lhs(1,1)*lhs(2,2) - lhs(1,2)*lhs(2,1)  ! determinant 
      invdet = 1./(det)

!     Direct inverse of a 2x2 matrix
!     Should probably use some iterative method for larger matrices      
      inv_lhs(1,1) = lhs(2,2)*invdet
      inv_lhs(2,2) = lhs(1,1)*invdet
      inv_lhs(1,2) = -lhs(1,2)*invdet
      inv_lhs(2,1) = -lhs(2,1)*invdet
   
      rhs(1) = inert*bd_psiv/DT + Fs
      rhs(2) = bd_psi/DT + Fadj

      call mxm(inv_lhs,2,rhs,2,psi_soln,1)

!     Velocity
      psiv_iter(it) = psi_soln(1)

!     Position
      psi_iter(it)  = psi_soln(2)

!     Acceleration
      psia_iter(it) = (bd(1)*psiv_iter(it)-bd_psiv)/DT

!     Calculate residual
      psi_res(it) = psi_iter(it) - psi_s 
      psiv_res(it) = psiv_iter(it) - psiv_s
      psia_res(it) = psia_iter(it) - psia_s

!     Check convergence
      rnorm = sqrt(psiv_res(it)**2 + psi_res(it)**2)

      nlfsi_rnorm(it) = rnorm 

      if (rnorm.le.nlfsi_tol) then
        nlfsi_ifconv = .true.
      endif  

      if (it.eq.nlfsi_miters.and..not.nlfsi_ifconv) then
        if (nio.eq.0) write(6,*) 'FSI unconverged', rnorm
        nlfsi_ifconv = .true.
      endif  

!     update omega if not converged
      ifbisection=.false.            ! use simple bisection
      if (nlfsi_ifconv) then
        
!       Copy lag array for psi.
        do ilag=lorder-1,2,-1
          psilag(ilag)=psilag(ilag-1)
        enddo
        psilag(1)=psi
        psi = psi_s

!       Copy lag array for psiv 
        do ilag=lorder-1,2,-1
          psivlag(ilag)=psivlag(ilag-1)
        enddo
        psivlag(1)=psiv
        psiv = psiv_s

!       Acceleration
!       Copy lag array for psia
        do ilag=lorder-1,2,-1
          psialag(ilag)=psialag(ilag-1)
        enddo
        psialag(1)=psia
        psia = (bd(1)*psiv - bd_psiv)/DT

!       Copy lag array for psivm
        do ilag=lorder-1,2,-1
          psivmlag(ilag)=psivmlag(ilag-1)
        enddo
        psivmlag(1)=psiv_m

!       Extrapolate velocity for next time step 
        psiv_s = ab(1)*psiv + ab(2)*psivlag(1) + ab(3)*psivlag(2)

!       Extrapolate position for next time step 
        psi_s = ab(1)*psi + ab(2)*psilag(1) + ab(3)*psilag(2)

        if (it.gt.1) then 
          nlfsi_omegai(it) = nlfsi_omegai(it-1)
        else
          nlfsi_omegai(it) = 0.
        endif  
        omega_prev = nlfsi_omegai(it-1)

      else        ! if not converged
        if (it.le.1) then
          omg = omega_prev          ! From last time step
          nlfsi_omegai(it) = omg
        elseif (ifbisection) then
!         using simple bisection for now  
          omg = 0.5
          nlfsi_omegai(it) = omg
        else

          tmp1 = psiv_res(it-1)*(psiv_res(it) - psiv_res(it-1))+
     $           psi_res(it-1)*(psi_res(it)   - psi_res(it-1))
          tmp2 = (psiv_res(it) - psiv_res(it-1))**2 +
     $           (psi_res(it)  - psi_res(it-1))**2
          omg = -nlfsi_omegai(it-1)*tmp1/tmp2

!         Can't be negative
          omg = max(omg,omegamin) 
          nlfsi_omegai(it) = omg  
        endif

!       Calculate new velocity projection
        psiv_s = psiv_s + omg*psiv_res(it)

!       Calculate position using backward difference
        psi_s = psi_s + omg*psi_res(it)
       
      endif 

      psiv_m = 0.       ! not needed

      Fk = -stiff*psiv
      Fd = damp*psiv
      Ft = Fs

!     Save forces for output 
      nlfsi_Ff   = Fs
      nlfsi_Fk   = Fk
      nlfsi_Fd   = Fd
      nlfsi_Fadj = Fadj
      nlfsi_Ft   = Ft

!     prabal
!      if (nio.eq.0) write(6,12) 'iters',it,psia_iter(it),
!     $         psia_res(it),tmp,omg
!   12 format(A5,1x,I3,1x,4(1pE14.7,1x))

      return
      end subroutine nlfsi_struc_solve_adj

!---------------------------------------------------------------------- 

      subroutine torque_calc_adj(scale,x0,ifdout,iftout,velx,vely,velz,
     $                            press,grvx,grvy,grvz)

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
!      INCLUDE 'TOTAL' 

      real flow_rate,base_flow,domain_length,xsec,scale_vf
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)


      real scale,x0(3),w1(0:maxobj)
      logical ifdout,iftout

!     Gradient of Base flow      
      real grvx(lx1,ly1,lz1,lelv,ldim)
      real grvy(lx1,ly1,lz1,lelv,ldim)
      real grvz(lx1,ly1,lz1,lelv,ldim)

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
      call mappr(pm1,press,xm0,ym0)       ! xm0,ym0 used as work arrays
c
c     Add mean_pressure_gradient.X to p:

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
                                        ! Compute sij
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
               ieg   = object(iobj,mem,1)       ! global no
               ifc   = object(iobj,mem,2)       ! face no
               if (gllnid(ieg).eq.nid) then ! this processor has a contribution
                  ie = gllel(ieg)
                  call clcdcm_adj(dgtq,sa,xm0,ym0,zm0,sij,pm1,vdiff,
     $                           grvx,grvy,grvz,ifc,ie)
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
      end subroutine torque_calc_adj
!-----------------------------------------------------------------------
      subroutine clcdcm_adj(dgtq,a,xm0,ym0,zm0,sij,pm1,visc,
     $                      grvx,grvy,grvz,f,e)
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

!     Gradient of Base flow      
      real grvx(lx1,ly1,lz1,lelv,ldim)
      real grvy(lx1,ly1,lz1,lelv,ldim)
      real grvz(lx1,ly1,lz1,lelv,ldim)

      real dudx,dvdx,dwdx
      real dudy,dvdy,dwdy
      real dudz,dvdz,dwdz

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

!        x-gradients         
         dudx = grvx(j1,j2,1,e,1)         
         dvdx = grvy(j1,j2,1,e,1)
         dwdx = grvz(j1,j2,1,e,1)
!        y-gradients         
         dudy = grvx(j1,j2,1,e,2)
         dvdy = grvy(j1,j2,1,e,2)
         dwdy = grvz(j1,j2,1,e,2)
!        z-gradients
         dudz = grvx(j1,j2,1,e,3)
         dvdz = grvy(j1,j2,1,e,3)
         dwdz = grvz(j1,j2,1,e,3)

!        Pressure drag            
         dg(1,1) = -dudx*pm1(j1,j2,1,e)*n1 
     $             -dvdx*pm1(j1,j2,1,e)*n2 
     $             -dwdx*pm1(j1,j2,1,e)*n3

         dg(2,1) = -dudy*pm1(j1,j2,1,e)*n1
     $             -dvdy*pm1(j1,j2,1,e)*n2
     $             -dwdy*pm1(j1,j2,1,e)*n3

         dg(3,1) = -dudz*pm1(j1,j2,1,e)*n1
     $             -dvdz*pm1(j1,j2,1,e)*n2
     $             -dwdz*pm1(j1,j2,1,e)*n3
     
!        Viscous drag
         dg(1,2) =  v*(s11*n1 + s12*n2 + s13*n3)*dudx
     $             +v*(s21*n1 + s22*n2 + s23*n3)*dvdx
     $             +v*(s31*n1 + s32*n2 + s33*n3)*dwdx

         dg(2,2) =  v*(s11*n1 + s12*n2 + s13*n3)*dudy
     $             +v*(s21*n1 + s22*n2 + s23*n3)*dvdy
     $             +v*(s31*n1 + s32*n2 + s33*n3)*dwdy

         dg(3,2) =  v*(s11*n1 + s12*n2 + s13*n3)*dudz
     $             +v*(s21*n1 + s22*n2 + s23*n3)*dvdz
     $             +v*(s31*n1 + s32*n2 + s33*n3)*dwdz

!        Prabal. Rotational forces not implemented yet.         
         r1 = 0. !xm0(j1,j2,1,e)
         r2 = 0. !ym0(j1,j2,1,e)
         r3 = 0. !zm0(j1,j2,1,e)

         do l=1,2
         do k=1,3
            dgtq(k,l) = dgtq(k,l) + dg(k,l)
         enddo
         enddo

         dgtq(1,3) = dgtq(1,3) + (r2*dg(3,1)-r3*dg(2,1)) ! pressure
         dgtq(2,3) = dgtq(2,3) + (r3*dg(1,1)-r1*dg(3,1)) ! torque
         dgtq(3,3) = dgtq(3,3) + (r1*dg(2,1)-r2*dg(1,1))

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

!         x-gradients         
          dudx = grvx(j1,j2,1,e,1)         
          dvdx = grvy(j1,j2,1,e,1)
          dwdx = 0.
!         y-gradients         
          dudy = grvx(j1,j2,1,e,2)
          dvdy = grvy(j1,j2,1,e,2)
          dwdy = 0.

!         Pressure drag         
          dg(1,1) = -dudx*pm1(j1,j2,1,e)*n1 
     $              -dvdx*pm1(j1,j2,1,e)*n2 

          dg(2,1) = -dudy*pm1(j1,j2,1,e)*n1
     $              -dvdy*pm1(j1,j2,1,e)*n2

          dg(3,1) = 0.
     
!         Viscous drag
          dg(1,2) =  v*(s11*n1 + s12*n2)*dudx
     $              +v*(s21*n1 + s22*n2)*dvdx

          dg(2,2) =  v*(s11*n1 + s12*n2)*dudy
     $              +v*(s21*n1 + s22*n2)*dvdy

          dg(3,2) =  0.                      

!         Prabal. Rotational forces not implemented yet.         
          r1 = 0. !xm0(j1,j2,1,e)
          r2 = 0. !ym0(j1,j2,1,e)
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
      end subroutine clcdcm_adj

!-----------------------------------------------------------------------
      subroutine nlfsi_meshv_adj

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'          ! xm1,ym1,zm1
      include 'MVGEOM_DEF'
      include 'MVGEOM'        ! wx,wy,wz
      include 'WING_MVMSH'
      include 'NLFSI'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'

      real ucx,ucy,ucz        ! mesh velocities
      real dx,dy,dz
      real psi_sx,psi_sy,psi_sz

      integer i,n
      real rot_sn(3,3)
      real rot_cs(3,3)
      real dxyz(3),r0(3)

      call rzero(rot_sn,9)            ! sine terms
      call rzero(rot_cs,9)            ! cosine terms
!     Matrix of sine terms
!     For a clockwise rotation for a positive psi
      if (nlfsi_ifrot) then
        rot_sn(1,2) =  psi_s
        rot_sn(2,1) = -psi_s
      else
        continue
      endif  

!     Matrix of cosine terms
!     In a linear case this is just Identity      
      rot_cs(1,1) = 1
      rot_cs(2,2) = 1
      rot_cs(3,3) = 1

      n = nx1*ny1*nz1*nelv

      do i=1,n                          ! Translational velocity
!       Current time step mesh velocity at the wall is the same
!       as the velocity at the wall
        if (.not.nlfsi_ifrot) then

!         Mesh motion                
          wx(i,1,1,1) = 0.
          wy(i,1,1,1) = 0.
          if (ndim.eq.3) wz(i,1,1,1) = 0.

!         Boundary condition            
          ucx =  0.
          ucy =  psiv_s
          ucz =  0.

          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz

        else
!         Rotational FSI

!         Distance from Rotational axis                
          dx = xm1(i,1,1,1) - nlfsi_x0
          dy = ym1(i,1,1,1) - nlfsi_y0
          dz = 0.

!         Mesh motion 
          wx(i,1,1,1) = 0.
          wy(i,1,1,1) = 0.
          if (ndim.eq.3) wz(i,1,1,1) = 0.

!         Calculate distance from axis   
          r0(1) = dx 
          r0(2) = dy
          r0(3) = dz

!         perturbation in position
          call mxm(rot_sn,3,r0,3,dxyz,1)
          psi_sx = dxyz(1)
          psi_sy = dxyz(2)
          psi_sz = dxyz(3) 

          ucx =  psiv_s*dy
          ucy =  -psiv_s*dx
          ucz =  0.
         
          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz
         
        endif
       
      enddo

      return
      end subroutine nlfsi_meshv_adj
!---------------------------------------------------------------------- 



