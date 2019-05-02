!======================================================================
!     Routines for fluid structure interaction implicit solutions
!     using fixed-point method with dynamic relaxation      
!     Author: Prabal S. Negi
!======================================================================       
      subroutine nlfsi_advancep

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'CTIMER_DEF'
      include 'CTIMER'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'MVGEOM_DEF'
      include 'MVGEOM'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'MASS_DEF'
      include 'MASS'
      include 'ADJOINT_DEF'
      include 'ADJOINT'       ! ifadj
      include 'NLFSI'

      integer igeom
      common /cgeom/ igeom    ! Apparently this common block is in NekNek

      integer icalld2
      save icalld2
      data icalld2 /0/

      integer lt
      parameter (lt=lx1*ly1*lz1*lelv)

      integer lt2
      parameter (lt2=lx2*ly2*lz2*lelv)

      character str*20
      integer ilag

      real wk1(lt),wk2(lt),wk3(lt)


      if (icalld2.eq.0) then
        if (ifadj) then
          if (nio.eq.0) write(6,*) 'FSI Adjoint perturbation mode'
        else
          if (nio.eq.0) write(6,*) 'FSI perturbation mode'
        endif
        icalld2=icalld+1
      endif  


      call nekgsync
      if (iftran) call settime
        
      call setsolv
      call comment
            
!     PN-2/PN-2 formulation
      call setprop

      nlfsi_ifconv = .false.
      nlfsi_iter = 0
      if (ifadj) then
        str='NL FSI Adj. Solve:  '
      else
        str='NL FSI Pert. Solve: '
      endif

      if (npert.gt.1) then
        if (nid.eq.0) write(6,*) 'NLFSI implemented only npert=1'
        call exitt
      endif  

      do jp=1,npert

!       Save current time-step velocity and pressure
        call opcopy(vx0,vy0,vz0,vxp(1,jp),vyp(1,jp),vzp(1,jp))
        call copy(pr0,prp(1,jp),lt2)
        call copy(bm0,bm1,lt)
        call opzero(nlfsi_wx,nlfsi_wy,nlfsi_wz)

        do while (.not.nlfsi_ifconv)
        
          nlfsi_iter = nlfsi_iter+1

          do igeom=1,ngeom

!           In the linear case the mesh does not move
!           Done in nlfsi_updcoor   
            if (ifgeom) then
              call nlfsi_gengeom (igeom)
              call geneig  (igeom)
            endif

            ifield = 1
            imesh = 1
            call unorm
            call settolv

            if (ifflow) call nlfsi_fluid_pert(igeom)

            if (igeom.ge.2) call vol_flow        ! check for fixed flow rate

          enddo

!         assuming jp=1 has the fsi perturbation fields
          if (jp.eq.1) then
            if (ifadj) then 
              call nlfsi_main_adj
            else
              call nlfsi_main
            endif
          endif   

!         Write to std output            
          if (nio.eq.0) write(6,11) str,istep,time,
     $               nlfsi_iter,nlfsi_timei,nlfsi_tol,
     $               nlfsi_rnorm(nlfsi_iter),psi_res(nlfsi_iter),
     $               psiv_res(nlfsi_iter),nlfsi_omegai(nlfsi_iter) 
   11   format(A20,I7,1x,1pE14.7E2,1x,I3,1x,1pE11.4E2,1x,1pE9.2E2,1x,
     $         4(1pE14.7E2,1x))


!         For debugging            
!          if (nlfsi_iter.eq.2) nlfsi_ifconv=.true.

          if (.not.nlfsi_ifconv) then       ! not converged
!           revert to old velocities                
            call opcopy(vxp(1,jp),vyp(1,jp),vzp(1,jp),vx0,vy0,vz0)
            call copy(prp(1,jp),pr0,lt2)
          endif 

          if (ifusermv.and.ifpert) then
            if (ifadj) then          
              call nlfsi_meshv_adj
            else
              call nlfsi_meshv_pert
            endif
          endif

        enddo       ! end of non-linear iterations

!       Assuming FSI is coupled only for jp=1
        if (jp.eq.1) then
          call nlfsi_rstsave
          call nlfsi_output
          
!         Update lag arrays            
          call nlfsi_lagmass
          call nlfsi_lagmshv(nelv)  ! mesh velocities
        endif 

        call nlfsi_lagfieldp
        call nlfsi_lagprp
        call nlfsi_lagbfp

      enddo       ! jp=1,npert 

      return
      end subroutine nlfsi_advancep
c-----------------------------------------------------------------------

      subroutine nlfsi_fluid_pert(igeom)
C
C     Compute pressure and velocity using consistent approximation spaces.     
C     Operator splitting technique.
C
C-----------------------------------------------------------------------

      INCLUDE 'SIZE_DEF'      
      INCLUDE 'SIZE'
      INCLUDE 'INPUT_DEF'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN_DEF'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'

!     prabal
      include 'FSI_DEBUG'

C
      COMMON /SCRNS/  RESV1 (LX1,LY1,LZ1,LELV)
     $ ,              RESV2 (LX1,LY1,LZ1,LELV)
     $ ,              RESV3 (LX1,LY1,LZ1,LELV)
     $ ,              DV1   (LX1,LY1,LZ1,LELV)
     $ ,              DV2   (LX1,LY1,LZ1,LELV)
     $ ,              DV3   (LX1,LY1,LZ1,LELV)
      COMMON /SCRVH/  H1    (LX1,LY1,LZ1,LELV)
     $ ,              H2    (LX1,LY1,LZ1,LELV)
C

      if (igeom.eq.1) then
c
c        Old geometry, old velocity
c
         call nlfsi_makefp
c
      else
c
c        New geometry, new velocity

         intype = -1
         call sethlm   (h1,h2,intype)
         call adjonbc_2 (h2)
!        prabal   
         call nlfsi_cresvipp (resv1,resv2,resv3,h1,h2)

         call ophinv   (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxh)

         call opadd2   (vxp(1,jp),vyp(1,jp),vzp(1,jp),dv1,dv2,dv3)

!         call incomprp (vxp(1,jp),vyp(1,jp),vzp(1,jp),prp(1,jp))
!         call incomprp_uzawa (vxp(1,jp),vyp(1,jp),vzp(1,jp),prp(1,jp))
         call nlfsi_incomprp_uzawa (vxp(1,jp),vyp(1,jp),vzp(1,jp)
     $                              ,prp(1,jp))

      endif

      return
      end subroutine nlfsi_fluid_pert 

c-----------------------------------------------------------------------
      subroutine nlfsi_makefp
c
c     Make rhs for velocity perturbation equation
c
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'
      include 'ADJOINT'
      include 'NLFSI'
                                                    call makeufp
      if (ifnav.and.(.not.ifchar).and.(.not.ifadj)) call advabp
      if (ifnav.and.(.not.ifchar).and.(     ifadj)) call advabp_adjoint
!     Routine changed so that we don't update lag arrays      
      if (iftran)                                   call nlfsi_makextp
                                                    call makebdfp



      return
      end subroutine nlfsi_makefp
c--------------------------------------------------------------------
      subroutine nlfsi_incomprp_uzawa (ux,uy,uz,up)
c
c     Project U onto the closest incompressible field
c
c     Input:  U     := (ux,uy,uz)
c
c     Output: updated values of U, iproj, proj; and
c             up    := pressure correction req'd to impose div U = 0
c
c
c     Dependencies: ifield ==> which "density" (vtrans) is used.
c
c     Comments (Mattias, 2014-12-23):
c     1.  This routine has been modified relative to incomprp to enable 
c         the Uzawa method, i.e. E=DH^(-1)D^T in addition to 
c         E=(dt/bd(1))DB^(-1)D^T. See e.g. Fischer (1997), JCP
c
c     2.  For Uzawa method no pressure projection is implemented!
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

c
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)
      COMMON /SCRCH/ PREXTR(LX2,LY2,LZ2,LELV)
      logical ifprjp

      real bdti,scaledt,scaledi,dtbd
c
      if (icalld.eq.0) tpres=0.0
      icalld=icalld+1
      npres=icalld
      etime1 = dnekclock()
c
      ntot1  = nx1*ny1*nz1*nelv
      ntot2  = nx2*ny2*nz2*nelv

      call opdiv   (dp,ux,uy,uz)
      if (IFUZAWA) then
         intype = -1

         intloc = -1
         call sethlm  (h1,h2,intloc)
         call rzero   (h2inv,ntot1)
      else
         intype = 1

         call rzero   (h1,ntot1)
         call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
         call invers2 (h2inv,h2,ntot1)

         bdti   = -bd(1)/dt
         call cmult  (dp,bdti,ntot2)
      endif
      call add2col2(dp,bm2,usrdiv,ntot2) ! User-defined divergence.
      call ortho   (dp)

      ifprjp=.false.    ! project out previous pressure solutions?

      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0.and.(.not.IFUZAWA)) 
     $     ifprjp=.true.

!     Most likely, the following can be commented out. (pff, 1/6/2010)
      if (npert.gt.1.or.ifbase) ifprjp=.false.

!     Projections      
      if (ifprjp)   call setrhs  (dp,h1,h2,h2inv)

!     Pressure solve      
      if (IFUZAWA) then
        call esolver (dp,h1,h2,h2inv,intype)
      else  
        scaledt = dt/bd(1)
        scaledi = 1./scaledt
        call cmult(dp,scaledt,ntot2)        ! scale for tol
        call esolver(dp,h1,h2,h2inv,intype)
        call cmult(dp,scaledi,ntot2)
      endif

!     Add back projections      
      if (ifprjp) call gensoln (dp,h1,h2,h2inv)

!     Correct pressure      
      call add2(up,dp,ntot2)

!     Correct velocity      
      call opgradt (w1 ,w2 ,w3 ,dp)
      if (IFUZAWA) then
        call ophinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h1 ,h2 ,tolhv ,nmxh)
        call opadd2  (ux ,uy ,uz ,dv1,dv2,dv3)
      else
        call opbinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
        dtb  = dt/bd(1)
        call opadd2cm (ux ,uy ,uz ,dv1,dv2,dv3, dtb )
      endif

!     timing      
      tpres=tpres+(dnekclock()-etime1)

c
!      call extrapprp(prextr)
!      call lagpresp
!      call add3(up,prextr,dp,ntot2)
c
      return
      end subroutine nlfsi_incomprp_uzawa
c------------------------------------------------------------------------

      subroutine nlfsi_makextp
c
c     Add extrapolation terms to perturbation source terms
c
c     (nek5 equivalent for velocity is "makeabf")
c
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'NLFSI'
C
      common /scrns/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)
c
      ntot1 = nx1*ny1*nz1*nelv
c

      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)
      call add3s2 (ta1,exx1p(1,jp),exx2p(1,jp),ab1,ab2,ntot1)
      call add3s2 (ta2,exy1p(1,jp),exy2p(1,jp),ab1,ab2,ntot1)
!      call copy   (exx2p(1,jp),exx1p(1,jp),ntot1)
!      call copy   (exy2p(1,jp),exy1p(1,jp),ntot1)
!      call copy   (exx1p(1,jp),bfxp (1,jp),ntot1)
!      call copy   (exy1p(1,jp),bfyp (1,jp),ntot1)
!     new array
      call copy(ABX0,BFXP(1,jp),NTOT1)
      call copy(ABY0,BFYP(1,jp),NTOT1)

      call add2s1 (bfxp(1,jp),ta1,ab0,ntot1)
      call add2s1 (bfyp(1,jp),ta2,ab0,ntot1)
      if (if3d) then
         call add3s2 (ta3,exz1p(1,jp),exz2p(1,jp),ab1,ab2,ntot1)
!         call copy   (exz2p(1,jp),exz1p(1,jp),ntot1)
!         call copy   (exz1p(1,jp),bfzp (1,jp),ntot1)
!        new array            
         call copy(ABZ0,BFZP(1,jp),NTOT1)
         call add2s1 (bfzp(1,jp),ta3,ab0,ntot1)
      endif
c
      return
      end subroutine nlfsi_makextp
c--------------------------------------------------------------------
      subroutine nlfsi_cresvipp (resv1,resv2,resv3,h1,h2)

c     Account for inhomogeneous Dirichlet boundary contributions 
c     in rhs of perturbation eqn.
c                                               n
c     Also, subtract off best estimate of grad p
c

      include 'SIZE'
      include 'TOTAL'

      real           resv1 (lx1,ly1,lz1,1)
      real           resv2 (lx1,ly1,lz1,1)
      real           resv3 (lx1,ly1,lz1,1)
      real           h1    (lx1,ly1,lz1,1)
      real           h2    (lx1,ly1,lz1,1)
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             prextr(lx2,ly2,lz2,lelv) 


      ntot1 = nx1*ny1*nz1*nelv
      ntot2 = nx2*ny2*nz2*nelv

      call bcdirvc (vxp(1,jp),vyp(1,jp),vzp(1,jp),
     $              v1mask,v2mask,v3mask)

!     prabal. Should check this. Does not seem to be included in
!     perturbation mode.
!      if (ifstrs)  call bcneutr

!     prlagp not updated here      
      call nlfsi_extrapp (prp(1,jp),prlagp(1,1,jp))

      call opgradt (resv1,resv2,resv3,prp(1,jp))
      CALL opadd2  (resv1,resv2,resv3,bfxp(1,jp),bfyp(1,jp),bfzp(1,jp))
      CALL ophx    (w1,w2,w3,vxp(1,jp),vyp(1,jp),vzp(1,jp),h1,h2)
      CALL opsub2  (resv1,resv2,resv3,w1,w2,w3)
C
      return
      end subroutine nlfsi_cresvipp
c-----------------------------------------------------------------------
      subroutine nlfsi_lagfieldp
c
c     Keep old Vp-field(s)
c
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'NLFSI'

      integer ilag

      do ilag=nbdinp-1,2,-1
         call opcopy
     $     (vxlagp(1,ilag  ,jp),vylagp(1,ilag  ,jp),vzlagp(1,ilag  ,jp)
     $     ,vxlagp(1,ilag-1,jp),vylagp(1,ilag-1,jp),vzlagp(1,ilag-1,jp))
      enddo
      call opcopy(vxlagp(1,1,jp),vylagp(1,1,jp),vzlagp(1,1,jp)
     $           ,vx0,vy0,vz0)
c
      return
      end subroutine nlfsi_lagfieldp
c-----------------------------------------------------------------------

      subroutine nlfsi_lagbfp

      implicit none

      INCLUDE 'SIZE_DEF'
      INCLUDE 'SIZE'
      INCLUDE 'SOLN_DEF'
      INCLUDE 'SOLN'
      INCLUDE 'NLFSI'

      INTEGER NTOT1

      NTOT1 = NX1*NY1*NZ1*NELV

      call copy   (exx2p(1,jp),exx1p(1,jp),ntot1)
      call copy   (exx1p(1,jp),ABX0,NTOT1)

      call copy   (exy2p(1,jp),exy1p(1,jp),ntot1)
      call copy   (exy1p(1,jp),ABY0,NTOT1)

      IF (NDIM.EQ.3) THEN
        call copy   (exz2p(1,jp),exz1p(1,jp),ntot1)
        call copy   (exz1p(1,jp),ABZ0,NTOT1)
      ENDIF

      return
      end subroutine nlfsi_lagbfp
c-----------------------------------------------------------------------

      subroutine nlfsi_lagprp

      implicit none

      INCLUDE 'SIZE_DEF'
      INCLUDE 'SIZE'
      INCLUDE 'SOLN_DEF'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'
      INCLUDE 'NLFSI'

      INTEGER NTOT2

      if (nbdinp.eq.3) then
        ntot2 = nx2*ny2*nz2*nelv
        call copy (prlagp(1,1,jp),pr0,ntot2)
      endif

      return
      end subroutine nlfsi_lagprp
c-----------------------------------------------------------------------

      subroutine nlfsi_torque_pert(scale,x0,ifdout,iftout,velx,
     $           vely,velz,press,velxp,velyp,velzp,pressp)
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
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'
      INCLUDE 'GEOM_DEF'
      INCLUDE 'GEOM'
      INCLUDE 'PARALLEL_DEF'
      INCLUDE 'PARALLEL'
      INCLUDE 'NLFSI' 

      real flow_rate,base_flow,domain_length,xsec,scale_vf
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)

c
      real x0(3),w1(0:maxobj)
      logical ifdout,iftout
c
      real pm1,sij,xm0,ym0,zm0
      common /scrns/         sij (lx1*ly1*lz1*(3*ldim-3)*lelv)
      common /scrcg/         pm1 (lx1,ly1,lz1,lelv)
      common /scrsf/         xm0(lx1,ly1,lz1,lelt)
     $,                      ym0(lx1,ly1,lz1,lelt)
     $,                      zm0(lx1,ly1,lz1,lelt)
c
      integer lr
      parameter (lr=lx1*ly1*lz1)
      real ur,us,ut,vr,vs,vt,wr,ws,wt
      common /scruz/         ur(lr),us(lr),ut(lr)
     $                     , vr(lr),vs(lr),vt(lr)
     $                     , wr(lr),ws(lr),wt(lr)
c

      logical ifpout                      ! if perturbation components
                                          ! output
      real dgtq_p(3,8)                    ! Saving all individual terms
      real dgtq_psum(3,8,maxobj)          ! Sum individual terms over
      real dgtq_wk(3,8,maxobj)            ! all elements

      real scale
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
      real vdiff(lt)

      real glmax,glmin        ! functions
      integer iglmax          ! functions

      real torq_timer
      real dnekclock          ! function

      real sa                 ! local area
      real sarea(0:maxobj)    ! total area

      integer icalld
      save icalld
      data icalld /0/

      integer i,i0,ie,ieg,ifc,ii,iobj,mem,memtot,n,nij

!     Perturbation variables
      integer lvp,lvp2
      parameter(lvp=lpx1*lpy1*lpz1*lpelv)
      parameter(lvp2=lpx2*lpy2*lpz2*lpelv)
      real velxp(lvp),velyp(lvp),velzp(lvp),pressp(lvp2)
      real pm0_wk(lt)
      real sij0(lx1*ly1*lz1*(3*ldim-3)*lelv)
      save sij0                           ! Save since this does not
                                          ! need to be recalculated 
                                          ! at each timestep

!     Check if all conditions are met for this subroutine      
!      if ((ifpert).and.(iffsi)) then
!!        continue
!      else
!        if (nio.eq.0) then
!          write(6,*) 'torque_calc_axis_pert called without proper 
!     $       conditions'
!          write(6,*) 'IFPERT/IFFSI/FSI_IFROT',ifpert,iffsi,fsi_ifrot
!        endif
!        call exitt
!      endif

      ifpout = .false. 

      torq_timer = dnekclock()      
c
      n = nx1*ny1*nz1*nelv
c
!      call mappr(pm1,pr,xm0,ym0) ! map pressure onto Mesh 1
      call mappr(pm0_wk,press,xm0,ym0)       ! xm0,ym0 used as work arrays
      call mappr(pm1,pressp,xm0,ym0)      ! xm0,ym0 used as work arrays

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

      call add2s2(pm0_wk,xm1,dpdx_mean,n)  ! Doesn't work if object is cut by 
      call add2s2(pm0_wk,ym1,dpdy_mean,n)  ! periodicboundary.  In this case,
      call add2s2(pm0_wk,zm1,dpdz_mean,n)  ! set ._mean=0 and compensate in
                                        ! Compute sij
c
      nij = 3
      if (if3d.or.ifaxis) nij=6

!     Only need to calculate mean stresses once.    
      if (icalld.eq.0) then
        call rzero(sij0,lt*(3*ldim-3))
        call comp_sij(sij0,nij,velx,vely,velz,
     $       ur,us,ut,vr,vs,vt,wr,ws,wt)
      endif
!     Perturbation stress
      call comp_sij(sij,nij,velxp,velyp,velzp,
     $     ur,us,ut,vr,vs,vt,wr,ws,wt)

c     Fill up viscous array w/ default
c

      call cfill(vdiff,param(2),n)
      if (icalld.eq.0) then
        xmx = glmax(xm1,n)
        xmn = glmin(xm1,n)
        ymx = glmax(ym1,n)
        ymn = glmin(ym1,n)
        zmx = glmax(zm1,n)
        zmn = glmin(zm1,n)
      endif

!     Routines are set up such that the rotation axis is along 'Z'      
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
      call rzero(dgtq_psum,24*maxobj)

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
                  call nlfsi_clcdcm_pert(dgtq,dgtq_p,sa,xm0,ym0,zm0,
     $                             sij,sij0,pm1,pm0_wk,vdiff,ifc,ie)

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
                  call cmult(dgtq_p,scale,24)
                  call add2(dgtq_psum(1,1,iobj),dgtq_p,24)

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

      call gop(dgtq_psum,dgtq_wk,'+  ',maxobj*24)
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
          write(6,*) 'Drag/Torque calculations'
          if (if3d.or.ifaxis) then
           if (ifdout) then
            write(6,16) istep,time,torq_timer,
     $              dragx(i),dragpx(i),dragvx(i),sarea(i),i,'dragx'
            write(6,16) istep,time,torq_timer,
     $              dragy(i),dragpy(i),dragvy(i),sarea(i),i,'dragy'
            write(6,16) istep,time,torq_timer,
     $              dragz(i),dragpz(i),dragvz(i),sarea(i),i,'dragz'
           endif
           if (iftout) then
            write(6,16) istep,time,torq_timer,
     $              torqx(i),torqpx(i),torqvx(i),sarea(i),i,'torqx'
            write(6,16) istep,time,torq_timer,
     $              torqy(i),torqpy(i),torqvy(i),sarea(i),i,'torqy'
            write(6,16) istep,time,torq_timer,
     $              torqz(i),torqpz(i),torqvz(i),sarea(i),i,'torqz'
           endif
          else
           if (ifdout) then
            write(6,16) istep,time,torq_timer,
     $              dragx(i),dragpx(i),dragvx(i),sarea(i),i,'dragx'
            write(6,16) istep,time,torq_timer,
     $              dragy(i),dragpy(i),dragvy(i),sarea(i),i,'dragy'
           endif
           if (iftout) then
            write(6,16) istep,time,torq_timer,
     $              torqz(i),torqpz(i),torqvz(i),sarea(i),i,'torqz'
           endif
          endif
        endif
   16   format(i8,1p6e16.8,1x,i3.1,a6)
      enddo

      if ((ifpout).and.(nio.eq.0)) then
        write(6,17) istep,time,(dgtq_psum(1,ii,1),ii=1,8),'dgtq_p1'
        write(6,17) istep,time,(dgtq_psum(2,ii,1),ii=1,8),'dgtq_p2'
        write(6,17) istep,time,(dgtq_psum(3,ii,1),ii=1,8),'dgtq_p3'
      endif
   17 format(i8,9(e17.8E3,1x),1x,a7)

      icalld = icalld+1

      return
      end subroutine nlfsi_torque_pert
!----------------------------------------------------------------------
      subroutine nlfsi_clcdcm_pert(dgtq,dgtq_p,a,xm0,ym0,zm0,
     $                        sij,sij0,pm1,pm0_wk,visc,f,e)

!     This routine must only be called when the following conditions are
!     met:
!     IFFSI = .TRUE.
!     IFPERT = .TRUE.
!
!     In this scenario, the routine calculates the linearized moments
!     along 'Z' on the defined object.            
            
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
      INCLUDE 'NLFSI'
c
      real dgtq(3,4)                      ! Total everything
      real dgtq_p(3,8)                    ! Saving all individual terms
      real xm0 (lx1,ly1,lz1,lelt)
      real ym0 (lx1,ly1,lz1,lelt)
      real zm0 (lx1,ly1,lz1,lelt)
      real sij (lx1,ly1,lz1,3*ldim-3,lelv)! Perturbation shear stress
      real sij0(lx1,ly1,lz1,3*ldim-3,lelv)! Mean Shear stress
      real pm1 (lx1,ly1,lz1,lelv)         ! perturbation pressure mapped
      real pm0_wk(lx1,ly1,lz1,lelv)         ! mean pressure mapped
      real visc(lx1,ly1,lz1,lelv)

      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1)
c

      real dg(3,6) !(p.n0),(Tauij.n0),(Tauij.[Sn].n0),
                   !(tauij.n0),(tauij.[Sn].n0),(P.[Sn].n0),
c
      integer f,e
      real    n1,n2,n3,n1p,n2p,n3p
      real n0(3),nr(3)        ! surface normals and rotated normals
      real r0(3),rr(3)        ! distance from axis. (Original and rotated)

      real a                  ! total (local) area

      integer i,l,k
      integer j1,j2,js1,js2,jf1,jf2,jskip1,jskip2
      integer pf
      real s11,s21,s31,s12,s22,s32,s13,s23,s33
      real v
      real r1,r2,r3,r1p,r2p,r3p

!     Perturbation rotational matrix      
      real rot_sn(3,3),rot_cs(3,3)


!     Check if all conditions are met for this subroutine      
      if (ifpert.and.ifnlfsi) then
!        continue
!          write(6,*) 'clcdcm_pert'
      else
        if (nio.eq.0) then
          write(6,*) 'clcdcm_pert called without proper conditions'
          write(6,*) 'IFPERT/IFNLFSI/FSI_IFROT',ifpert,ifnlfsi,
     $                   nlfsi_ifrot
        endif
        call exitt
      endif      

!     In this case cosider rotation about 'Z' axis
!     Rewrite rotation matrix as sum of two matrices.
!     One with just the cosine terms and one with just the sine terms.      
!     Rot = Cs + Sn
!     With linearization: 
!           psi << 1      
!           sin(psi) ~ psi
!           cos(psi) ~ 1
!     The Cs becomes an Identity matrix.
!     The Sn only has psi terms.
!     Sn can therefore be considered a 'perturbation matrix'

      call rzero(rot_sn,9)            ! sine terms
      call rzero(rot_cs,9)            ! cosine terms
!     Matrix of sine terms
!     For a clockwise rotation for a positive psi
      if (nlfsi_ifrot) then
        rot_sn(1,2) =  psi_s
        rot_sn(2,1) = -psi_s
      else
!       If there is only translation motion,
!       Then the direction of normals does not change
!       Also the moment arm of points on the surface does not change
!       As long as it is a rigid body motion. 
        continue
      endif  

!     Matrix of cosine terms 
      rot_cs(1,1) = 1
      rot_cs(2,2) = 1
      rot_cs(3,3) = 1

!     Since rot_cs is Identity. We skip all multiplications involving
!     rot_cs.      

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
      call rzero(dgtq_p,24)
      call rzero(dg,12)
      call rzero(nr,3)
      call rzero(rr,3)
c
      if (if3d.or.ifaxis) then
       i = 0
       a = 0
       do j2=js2,jf2,jskip2
       do j1=js1,jf1,jskip1
         i = i+1
         n0(1) = unx(i,1,f,e)*area(i,1,f,e)
         n0(2) = uny(i,1,f,e)*area(i,1,f,e)
         n0(3) = unz(i,1,f,e)*area(i,1,f,e)

!        Calculate rotated normals            
         call mxm(rot_sn,3,n0,3,nr,1)
         n1p = nr(1)
         n2p = nr(2)
         n3p = nr(3)

!        Regular normals   
         n1 = n0(1)
         n2 = n0(2)
         n3 = n0(3)

!        Calculate distance from axis   
         r0(1) = xm0(j1,j2,1,e)
         r0(2) = ym0(j1,j2,1,e)
         r0(3) = zm0(j1,j2,1,e)

!        Rotated vectors   
         call mxm(rot_sn,3,r0,3,rr,1)
         r1p = rr(1)
         r2p = rr(2)
         r3p = rr(3)
!        original distance
         r1 = r0(1)
         r2 = r0(2)
         r3 = r0(3)

         a  = a + area(i,1,f,e)
         v  = visc(j1,j2,1,e)
c
!        Mean Shear Stesses             
         s11 = sij0(j1,j2,1,1,e)
         s21 = sij0(j1,j2,1,4,e)
         s31 = sij0(j1,j2,1,6,e)
c
         s12 = sij0(j1,j2,1,4,e)
         s22 = sij0(j1,j2,1,2,e)
         s32 = sij0(j1,j2,1,5,e)
c
         s13 = sij0(j1,j2,1,6,e)
         s23 = sij0(j1,j2,1,5,e)
         s33 = sij0(j1,j2,1,3,e)
c
!        (Perturbation pressure)*(Original normals)
!        p.n0        
         dg(1,1) = pm1(j1,j2,1,e)*n1
         dg(2,1) = pm1(j1,j2,1,e)*n2
         dg(3,1) = pm1(j1,j2,1,e)*n3
c
!        (Mean Stress)*(Original normals)            
!        Tauij*n0 
         dg(1,2) = -v*(s11*n1 + s12*n2 + s13*n3) ! viscous drag
         dg(2,2) = -v*(s21*n1 + s22*n2 + s23*n3)
         dg(3,2) = -v*(s31*n1 + s32*n2 + s33*n3)

!        (Mean Stress)*(Rotated normals)            
!        Tauij*[Sn]*n0            
         dg(1,3) = -v*(s11*n1p + s12*n2p + s13*n3p) ! viscous drag
         dg(2,3) = -v*(s21*n1p + s22*n2p + s23*n3p)
         dg(3,3) = -v*(s31*n1p + s32*n2p + s33*n3p)

!        Perturbation Shear Stesses             
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

!        (Perturbation Stress)*(Original normals)            
!        tauij*n0            
         dg(1,4) = -v*(s11*n1 + s12*n2 + s13*n3) ! viscous drag
         dg(2,4) = -v*(s21*n1 + s22*n2 + s23*n3)
         dg(3,4) = -v*(s31*n1 + s32*n2 + s33*n3)

!        (Mean pressure)*(Rotated normals)
!        P.[Sn].n0        
         dg(1,5) = pm0_wk(j1,j2,1,e)*n1p 
         dg(2,5) = pm0_wk(j1,j2,1,e)*n2p
         dg(3,5) = pm0_wk(j1,j2,1,e)*n3p

!        (Perturbation Stress)*(Rotated normals) ! Non-linear term           
!        tauij*[Sn]*n0            
         dg(1,6) = -v*(s11*n1p + s12*n2p + s13*n3p)
         dg(2,6) = -v*(s21*n1p + s22*n2p + s23*n3p)
         dg(3,6) = -v*(s31*n1p + s32*n2p + s33*n3p)

!        Sum up pressure and viscous contributions to drag            
         do k=1,3
!           For linearized drag only two components play a role.
            dgtq(k,1) = dgtq(k,1) + dg(k,1)
            dgtq(k,2) = dgtq(k,2) + dg(k,3)+ dg(k,4)
!           Mean component (Sij.n0) is needed for calculating perturbation
!           moments. 

!           Individual components
            dgtq_p(k,1) = dgtq_p(k,1) + dg(k,1)
            dgtq_p(k,2) = dgtq_p(k,2) + dg(k,2)
            dgtq_p(k,3) = dgtq_p(k,3) + dg(k,3)
            dgtq_p(k,4) = dgtq_p(k,4) + dg(k,4)
         enddo

!        Sum up pressure and viscous contributions to torque            
!        Perturbation pressure torque.            
         dgtq(1,3) = dgtq(1,3) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq(2,3) = dgtq(2,3) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq(3,3) = dgtq(3,3) + (r1*dg(2,1)-r2*dg(1,1))

!        Saving individual components            
         dgtq_p(1,5) = dgtq_p(1,5) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq_p(2,5) = dgtq_p(2,5) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq_p(3,5) = dgtq_p(3,5) + (r1*dg(2,1)-r2*dg(1,1))

!!       Viscous torque
!        (Rotated normals) x [(Mean Stress)*(Original normals)]
!        [Sn*r0] x (Tauij*n0)            
         dgtq(1,4) = dgtq(1,4) + (r2p*dg(3,2)-r3p*dg(2,2))
         dgtq(2,4) = dgtq(2,4) + (r3p*dg(1,2)-r1p*dg(3,2))
         dgtq(3,4) = dgtq(3,4) + (r1p*dg(2,2)-r2p*dg(1,2))

!        Saving individual components            
         dgtq_p(1,6) = dgtq_p(1,6) + (r2p*dg(3,2)-r3p*dg(2,2))
         dgtq_p(2,6) = dgtq_p(2,6) + (r3p*dg(1,2)-r1p*dg(3,2))
         dgtq_p(3,6) = dgtq_p(3,6) + (r1p*dg(2,2)-r2p*dg(1,2))

!        (Original normals) x [(Mean Stress)*(Rotated normals)]
!        (r0) x (Tauij*[Sn]*n0)            
         dgtq(1,4) = dgtq(1,4) + (r2*dg(3,3)-r3*dg(2,3))
         dgtq(2,4) = dgtq(2,4) + (r3*dg(1,3)-r1*dg(3,3))
         dgtq(3,4) = dgtq(3,4) + (r1*dg(2,3)-r2*dg(1,3))

!        Saving individual components            
         dgtq_p(1,7) = dgtq_p(1,7) + (r2*dg(3,3)-r3*dg(2,3))
         dgtq_p(2,7) = dgtq_p(2,7) + (r3*dg(1,3)-r1*dg(3,3))
         dgtq_p(3,7) = dgtq_p(3,7) + (r1*dg(2,3)-r2*dg(1,3))
        
!        (Original normals) x [(Perturbation Stress)*(Original normals)]
!        (r0) x (tauij*n0)            
         dgtq(1,4) = dgtq(1,4) + (r2*dg(3,4)-r3*dg(2,4))
         dgtq(2,4) = dgtq(2,4) + (r3*dg(1,4)-r1*dg(3,4))
         dgtq(3,4) = dgtq(3,4) + (r1*dg(2,4)-r2*dg(1,4))

!        Saving individual components            
         dgtq_p(1,8) = dgtq_p(1,8) + (r2*dg(3,4)-r3*dg(2,4))
         dgtq_p(2,8) = dgtq_p(2,8) + (r3*dg(1,4)-r1*dg(3,4))
         dgtq_p(3,8) = dgtq_p(3,8) + (r1*dg(2,4)-r2*dg(1,4))

       enddo
       enddo

      else ! 2D

       i = 0
       a = 0
       do j2=js2,jf2,jskip2
       do j1=js1,jf1,jskip1
         i = i+1
         n0(1) = unx(i,1,f,e)*area(i,1,f,e)
         n0(2) = uny(i,1,f,e)*area(i,1,f,e)
         n0(3) = 0.

!        Calculate rotated normals            
         call mxm(rot_sn,3,n0,3,nr,1)
         n1p = nr(1)
         n2p = nr(2)
         n3p = 0.

!        Regular normals   
         n1 = n0(1)
         n2 = n0(2)
         n3 = n0(3)

!        Calculate distance from axis   
         r0(1) = xm0(j1,j2,1,e)
         r0(2) = ym0(j1,j2,1,e)
         r0(3) = 0.

!        Rotated vectors   
         call mxm(rot_sn,3,r0,3,rr,1)
         r1p = rr(1)
         r2p = rr(2)
         r3p = 0.
!        original distance
         r1 = r0(1)
         r2 = r0(2)
         r3 = 0.

         a  = a + area(i,1,f,e)
         v  = visc(j1,j2,1,e)

!        Mean Shear stresses            
         s11 = sij0(j1,j2,1,1,e)
         s12 = sij0(j1,j2,1,3,e)
         s21 = sij0(j1,j2,1,3,e)
         s22 = sij0(j1,j2,1,2,e)

!        (Perturbation pressure)*(Original normals)
!        p.n0        
         dg(1,1) = pm1(j1,j2,1,e)*n1
         dg(2,1) = pm1(j1,j2,1,e)*n2
         dg(3,1) = 0.
c
!        (Mean Stress)*(Original normals)            
!        Tauij*n0            
         dg(1,2) = -v*(s11*n1 + s12*n2)
         dg(2,2) = -v*(s21*n1 + s22*n2)
         dg(3,2) = 0.

!        (Mean Stress)*(Rotated normals)            
!        Tauij*[Sn]*n0            
         dg(1,3) = -v*(s11*n1p + s12*n2p)
         dg(2,3) = -v*(s21*n1p + s22*n2p)
         dg(3,3) = 0. 

!        Perturbation Shear stress            
         s11 = sij(j1,j2,1,1,e)
         s12 = sij(j1,j2,1,3,e)
         s21 = sij(j1,j2,1,3,e)
         s22 = sij(j1,j2,1,2,e)

!        (Perturbation Stress)*(Original normals)            
!        tauij*n0            
         dg(1,4) = -v*(s11*n1 + s12*n2)
         dg(2,4) = -v*(s21*n1 + s22*n2)
         dg(3,4) = 0. 

!        (Mean pressure)*(Rotated normals)
!        P.[Sn].n0        
         dg(1,5) = pm0_wk(j1,j2,1,e)*n1p
         dg(2,5) = pm0_wk(j1,j2,1,e)*n2p
         dg(3,5) = 0.

!        (Perturbation Stress)*(Rotated normals) ! Non-linear term           
!        tauij*[Sn]*n0            
         dg(1,6) = -v*(s11*n1p + s12*n2p + s13*n3p)
         dg(2,6) = -v*(s21*n1p + s22*n2p + s23*n3p)
         dg(3,6) = 0.

!        Sum up pressure and viscous contributions to drag            
         do k=1,3
            dgtq(k,1) = dgtq(k,1) + dg(k,1)
!           For linearized drag only two components play a role.
!           Mean component is needed for calculating perturbation
!           moments.            
            dgtq(k,2) = dgtq(k,2) + dg(k,3) + dg(k,4)
           
!           Individual components
            dgtq_p(k,1) = dgtq_p(k,1) + dg(k,1)
            dgtq_p(k,2) = dgtq_p(k,2) + dg(k,2)
            dgtq_p(k,3) = dgtq_p(k,3) + dg(k,3)
            dgtq_p(k,4) = dgtq_p(k,4) + dg(k,4)

         enddo

!        Sum up pressure and viscous contributions to torque            
!        Perturbation pressure torque.            
         dgtq(1,3) = 0. ! dgtq(1,3) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq(2,3) = 0. ! dgtq(2,3) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq(3,3) = dgtq(3,3) + (r1*dg(2,1)-r2*dg(1,1))

!        Saving individual components            
         dgtq_p(1,5) = 0. ! dgtq_p(1,5) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq_p(2,5) = 0. ! dgtq_p(2,5) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq_p(3,5) = dgtq_p(3,5) + (r1*dg(2,1)-r2*dg(1,1))


!!       Viscous torque
!        (Rotated normals) x [(Mean Stress)*(Original normals)]
!        [Sn*r0] x (Tauij*n0)            
         dgtq(1,4) = 0. ! dgtq(1,4) + (r2p*dg(3,2)-r3p*dg(2,2))
         dgtq(2,4) = 0. ! dgtq(2,4) + (r3p*dg(1,2)-r1p*dg(3,2))
         dgtq(3,4) = dgtq(3,4) + (r1p*dg(2,2)-r2p*dg(1,2))

!        Saving individual components            
         dgtq_p(1,6) = 0. ! dgtq_p(1,6) + (r2p*dg(3,2)-r3p*dg(2,2))
         dgtq_p(2,6) = 0. ! dgtq_p(2,6) + (r3p*dg(1,2)-r1p*dg(3,2))
         dgtq_p(3,6) = dgtq_p(3,6) + (r1p*dg(2,2)-r2p*dg(1,2))

!        (Original normals) x [(Mean Stress)*(Rotated normals)]
!        (r0) x (Tauij*[Sn]*n0)            
         dgtq(1,4) = 0. ! dgtq(1,4) + (r2*dg(3,3)-r3*dg(2,3))
         dgtq(2,4) = 0. ! dgtq(2,4) + (r3*dg(1,3)-r1*dg(3,3))
         dgtq(3,4) = dgtq(3,4) + (r1*dg(2,3)-r2*dg(1,3))

!        Saving individual components            
         dgtq_p(1,7) = 0. ! dgtq_p(1,7) + (r2*dg(3,3)-r3*dg(2,3))
         dgtq_p(2,7) = 0. ! dgtq_p(2,7) + (r3*dg(1,3)-r1*dg(3,3))
         dgtq_p(3,7) = dgtq_p(3,7) + (r1*dg(2,3)-r2*dg(1,3))

!        (Original normals) x [(Perturbation Stress)*(Original normals)]
!        (r0) x (tauij*n0)            
         dgtq(1,4) = 0. ! dgtq(1,4) + (r2*dg(3,4)-r3*dg(2,4))
         dgtq(2,4) = 0. ! dgtq(2,4) + (r3*dg(1,4)-r1*dg(3,4))
         dgtq(3,4) = dgtq(3,4) + (r1*dg(2,4)-r2*dg(1,4))

!        Saving individual components            
         dgtq_p(1,8) = 0. ! dgtq_p(1,8) + (r2*dg(3,4)-r3*dg(2,4))
         dgtq_p(2,8) = 0. ! dgtq_p(2,8) + (r3*dg(1,4)-r1*dg(3,4))
         dgtq_p(3,8) = dgtq_p(3,8) + (r1*dg(2,4)-r2*dg(1,4))

       enddo
       enddo
      endif       ! if 3d

      return
      end subroutine nlfsi_clcdcm_pert

!-----------------------------------------------------------------------
      




