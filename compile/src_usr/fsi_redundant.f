      subroutine fsi_nek_advancep

      implicit none

      INCLUDE 'SIZE_DEF'
      INCLUDE 'SIZE'
      INCLUDE 'SOLN_DEF'
      INCLUDE 'SOLN'
      INCLUDE 'INPUT_DEF'
      INCLUDE 'INPUT'
      INCLUDE 'CTIMER_DEF'
      INCLUDE 'CTIMER'
      INCLUDE 'GEOM_DEF'
      INCLUDE 'GEOM'
      INCLUDE 'MVGEOM_DEF'
      INCLUDE 'MVGEOM'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'
      INCLUDE 'FSI'

      integer igeom
      common /cgeom/ igeom    ! Apparently this common block is in NekNek

      integer icalld2
      save icalld2
      data icalld2 /0/

!     Variation in base-state due to transport      

      integer lt2
      parameter (lt2=lx2*ly2*lz2*lelv)

      real dttmp
      
!     Only for FSI perturbation solve
      if (.not.ifpert.and.iffsi) then
        if (nio.eq.0) 
     $    write(6,*) 'only for FSI Linear Stability'
        call exitt
      endif  
!     Not enabled for MHD
      if (ifmhd) then
        if (nio.eq.0) 
     $    write(6,*) 'MHD not enabled for FSI Linear Stability'
        call exitt
      endif  

!     Not enabled for NekNek
      if (ifneknekm) then
        if (nio.eq.0) 
     $    write(6,*) 'NekNek not enabled for FSI Linear Stability'
        call exitt
      endif  
!     Not enabled for ifheat
      if (ifheat) then
        if (nio.eq.0) 
     $    write(6,*) 'ifheat not enabled for FSI Linear Stability'
        call exitt
      endif  
!     Not enabled for CMT
      if (ifcmt) then
        if (nio.eq.0) 
     $    write(6,*) 'CMT not enabled for FSI Linear Stability'
        call exitt
      endif
!     Not enabled for PN/PN
      if (ifsplit) then
        if (nio.eq.0) 
     $    write(6,*) 'PnPn not enabled for FSI Linear Stability'
        call exitt
      endif

!     Npert.gt.1 not implemented      
      if (npert.gt.2) then
        if (nio.eq.0) 
     $    write(6,*) 'Only implemented for FSI npert<=2'
        call exitt
      endif

      call nekgsync
      if (iftran) call settime
        
      call setsolv
      call comment
            
!     PN-2/PN-2 formulation
!     Perturbation solve
      call setprop

      do igeom=1,ngeom

!       In perturbation mode we don't change the mesh 
        if (ifgeom) then
          call fsi_gengeom (igeom)
          call geneig  (igeom)
        endif

        do jp=1,npert

          fsi_jp = jp
          if (nio.eq.0.and.igeom.eq.2) write(6,1) istep,time,jp
   1      format(i9,1pe14.7,' Perturbation Solve:',i5)
     
          if (ifflow) call perturbv(igeom)

        enddo     ! jp=1,npert
        jp = 0
        fsi_jp=jp

      enddo     ! igeom

!     Correct for Added-mass effects
      call fsi_main_pert

      call fsi_meshv_pert2

      return
      end subroutine fsi_nek_advancep
c-----------------------------------------------------------------------


      subroutine fsi_main_pert_1

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
      include 'FSI'
      include 'CHKPOINT'

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

      integer i,lvp,lvp2
      parameter(lvp=lpx1*lpy1*lpz1*lpelv)
      parameter(lvp2=lpx2*lpy2*lpz2*lpelv)
      real wk_vxp(lvp),wk_vyp(lvp),wk_vzp(lvp),wk_prp(lvp2)

      real wk_vx0(lvp),wk_vy0(lvp),wk_vz0(lvp),wk_pr0(lvp2)
      real tmp_pr(lvp)

      integer jpold

      real*8 dnekclock          ! function. timing
      real*8 steptime           ! timing for step

      real rot_sn(3,3)
      real rot_cs(3,3)
      real dxyz(3),r0(3)
      real eta_sx,eta_sy,eta_sz

      real rxyz 
      common /scrsf/ rxyz(lx1*ly1*lz1*lelt,3)


!     Start timing for this step   
      steptime = dnekclock()

      if (icalld.eq.0) then
!       Initialize timer            
        fsi_timea = 0.
        fsi_timei = 0.

        call fsi_init

        icalld=icalld+1
      endif  

      if (.not.IFFSI) then
        etav_s = 0.
        etav   = 0.
        etaa   = 0.
        eta_s  = 0.
        if (ifusermv) call fsi_meshv

        return
      endif  

      if (ifchkptrst.and.istep.lt.chkptnrsf) then
!       Read saved values for restart
!       If the simulation has been restarted
        call fsi_rstread
        if (ifusermv.and.ifpert) call fsi_meshv_pert2

!       placeholder             
        etaa = 0.
        eta_s= eta

        call fsi_output

        fsi_timei=dnekclock()-steptime
        fsi_timea=fsi_timea+fsi_timei
        return

      elseif (istep.eq.0) then
!       If there are no restarts              
!       At first time step the projected velocity is
!       just the initialized velocity
        etav_s = etav_ini
        etav   = etav_ini
        eta    = eta_ini

!       prabal. placeholders 
        etaa = 0.
        eta_s  = etav_s*DT

        if (ifusermv.and.ifpert) call fsi_meshv_pert2

        call fsi_output

        fsi_timei=dnekclock()-steptime
        fsi_timea=fsi_timea+fsi_timei
        return
      endif

      x0(1) = fsi_x0
      x0(2) = fsi_y0
      if (ndim.eq.3) x0(3) = 1000.         ! set outside the domain

      ifdout=.false.
      iftout=.false.

      call opzero(wk_vxp,wk_vyp,wk_vzp)
      call rzero(wk_prp,lvp2)
      scale=fsi_rescale
      if (fsi_base_forc) then  
!       base flow forces 
        call torque_calc_axis_pert(scale,x0,ifdout,iftout,
     $      vx,vy,vz,pr,wk_vxp,wk_vyp,wk_vzp,wk_prp)

        call copy(fb_dragx (0),dragx (0),maxobj+1)
        call copy(fb_dragpx(0),dragpx(0),maxobj+1)
        call copy(fb_dragvx(0),dragpx(0),maxobj+1)
        call copy(fb_dragy (0),dragy (0),maxobj+1)
        call copy(fb_dragpy(0),dragpy(0),maxobj+1)
        call copy(fb_dragvy(0),dragvy(0),maxobj+1)
        call copy(fb_dragz (0),dragz (0),maxobj+1)
        call copy(fb_dragpz(0),dragpz(0),maxobj+1)
        call copy(fb_dragvz(0),dragvz(0),maxobj+1)
        call copy(fb_torqx (0),torqx (0),maxobj+1)
        call copy(fb_torqpx(0),torqpx(0),maxobj+1)
        call copy(fb_torqvx(0),torqpx(0),maxobj+1)
        call copy(fb_torqy (0),torqy (0),maxobj+1)
        call copy(fb_torqpy(0),torqpy(0),maxobj+1)
        call copy(fb_torqvy(0),torqvy(0),maxobj+1)
        call copy(fb_torqz (0),torqz (0),maxobj+1)
        call copy(fb_torqpz(0),torqpz(0),maxobj+1)
        call copy(fb_torqvz(0),torqvz(0),maxobj+1)
      else
        call rzero(fb_dragx (0),maxobj+1)
        call rzero(fb_dragpx(0),maxobj+1)
        call rzero(fb_dragvx(0),maxobj+1)
        call rzero(fb_dragy (0),maxobj+1)
        call rzero(fb_dragpy(0),maxobj+1)
        call rzero(fb_dragvy(0),maxobj+1)
        call rzero(fb_dragz (0),maxobj+1)
        call rzero(fb_dragpz(0),maxobj+1)
        call rzero(fb_dragvz(0),maxobj+1)
        call rzero(fb_torqx (0),maxobj+1)
        call rzero(fb_torqpx(0),maxobj+1)
        call rzero(fb_torqvx(0),maxobj+1)
        call rzero(fb_torqy (0),maxobj+1)
        call rzero(fb_torqpy(0),maxobj+1)
        call rzero(fb_torqvy(0),maxobj+1)
        call rzero(fb_torqz (0),maxobj+1)
        call rzero(fb_torqpz(0),maxobj+1)
        call rzero(fb_torqvz(0),maxobj+1)
      endif ! fsi_base_forc 

!     Assuming intrinsic perturbations are on jp=1
      call torque_calc_axis_pert(scale,x0,ifdout,iftout,
     $      wk_vxp,wk_vyp,wk_vzp,wk_prp,
     $      vxp(1,1),vyp(1,1),vzp(1,1),prp(1,1))

      call copy(fs_dragx (0),dragx (0),maxobj+1)
      call copy(fs_dragpx(0),dragpx(0),maxobj+1)
      call copy(fs_dragvx(0),dragpx(0),maxobj+1)
      call copy(fs_dragy (0),dragy (0),maxobj+1)
      call copy(fs_dragpy(0),dragpy(0),maxobj+1)
      call copy(fs_dragvy(0),dragvy(0),maxobj+1)
      call copy(fs_dragz (0),dragz (0),maxobj+1)
      call copy(fs_dragpz(0),dragpz(0),maxobj+1)
      call copy(fs_dragvz(0),dragvz(0),maxobj+1)
      call copy(fs_torqx (0),torqx (0),maxobj+1)
      call copy(fs_torqpx(0),torqpx(0),maxobj+1)
      call copy(fs_torqvx(0),torqpx(0),maxobj+1)
      call copy(fs_torqy (0),torqy (0),maxobj+1)
      call copy(fs_torqpy(0),torqpy(0),maxobj+1)
      call copy(fs_torqvy(0),torqvy(0),maxobj+1)
      call copy(fs_torqz (0),torqz (0),maxobj+1)
      call copy(fs_torqpz(0),torqpz(0),maxobj+1)
      call copy(fs_torqvz(0),torqvz(0),maxobj+1)

!     Calculate Impulsive response (Green's function)      
      call fsi_plan3

!     AMP forces on the object
      scale=fsi_rescale
      call opzero(wk_vxp,wk_vyp,wk_vzp)
      call rzero(wk_prp,lvp2) 
      call torque_calc_axis_pert(scale,x0,ifdout,iftout,
     $      wk_vxp,wk_vyp,wk_vzp,wk_prp,
     $      amp_vx,amp_vy,amp_vz,amp_pr)

      call copy(fg_dragx (0),dragx (0),maxobj+1)
      call copy(fg_dragpx(0),dragpx(0),maxobj+1)
      call copy(fg_dragvx(0),dragpx(0),maxobj+1)
      call copy(fg_dragy (0),dragy (0),maxobj+1)
      call copy(fg_dragpy(0),dragpy(0),maxobj+1)
      call copy(fg_dragvy(0),dragvy(0),maxobj+1)
      call copy(fg_dragz (0),dragz (0),maxobj+1)
      call copy(fg_dragpz(0),dragpz(0),maxobj+1)
      call copy(fg_dragvz(0),dragvz(0),maxobj+1)
      call copy(fg_torqx (0),torqx (0),maxobj+1)
      call copy(fg_torqpx(0),torqpx(0),maxobj+1)
      call copy(fg_torqvx(0),torqpx(0),maxobj+1)
      call copy(fg_torqy (0),torqy (0),maxobj+1)
      call copy(fg_torqpy(0),torqpy(0),maxobj+1)
      call copy(fg_torqvy(0),torqvy(0),maxobj+1)
      call copy(fg_torqz (0),torqz (0),maxobj+1)
      call copy(fg_torqpz(0),torqpz(0),maxobj+1)
      call copy(fg_torqvz(0),torqvz(0),maxobj+1)

!     Displacement forces on the object
      if (ifpert.and.fsi_ifdispl.and.npert.eq.1) then

!       Zero contribution for Mean field   
        call opzero(wk_vx0,wk_vy0,wk_vz0)
        call rzero(wk_pr0,lvp2) 

        if (fsi_ifrot) then

          call rzero(rot_sn,9)            ! sine terms
          call rzero(rot_cs,9)            ! cosine terms
!         Matrix of sine terms
!         For a clockwise rotation for a positive eta
          rot_sn(1,2) =  eta_s
          rot_sn(2,1) = -eta_s

!         Matrix of cosine terms
!         In a linear case this is just Identity      
          rot_cs(1,1) = 1
          rot_cs(2,2) = 1
          rot_cs(3,3) = 1

!         Distance from Rotational axis
          call opcopy(rxyz(1,1),rxyz(1,2),rxyz(1,3),xm1,ym1,zm1) 
          call cadd(rxyz(1,1),-fsi_x0,lvp)
          call cadd(rxyz(1,2),-fsi_y0,lvp)
          call rzero(rxyz(1,3),lvp)       ! rotation is about 'z'

          do i=1,lvp 
            r0(1) = rxyz(i,1) 
            r0(2) = rxyz(i,2)
            r0(3) = rxyz(i,3)

!           perturbation in position
            call mxm(rot_sn,3,r0,3,dxyz,1)
            eta_sx = dxyz(1)
            eta_sy = dxyz(2)
            eta_sz = dxyz(3) 

!           Calculate Grad(U0).eta            
            wk_vxp(i) =  fsi_grvx0(i,1,1,1,1)*eta_sx + 
     $           fsi_grvx0(i,1,1,1,2)*eta_sy +
     $           fsi_grvx0(i,1,1,1,3)*eta_sz
            wk_vyp(i) =  fsi_grvy0(i,1,1,1,1)*eta_sx + 
     $           fsi_grvy0(i,1,1,1,2)*eta_sy +
     $           fsi_grvy0(i,1,1,1,3)*eta_sz
            if (if3d) then
              wk_vzp(i) = fsi_grvz0(i,1,1,1,1)*eta_sx + 
     $               fsi_grvz0(i,1,1,1,2)*eta_sy +
     $               fsi_grvz0(i,1,1,1,3)*eta_sz
            else
              wk_vzp(i)=0.
            endif

!           Calculate Grad(P0).eta            
            tmp_pr(i) =  fsi_grpr0(i,1,1,1,1)*eta_sx + 
     $           fsi_grpr0(i,1,1,1,2)*eta_sy +
     $           fsi_grpr0(i,1,1,1,3)*eta_sz
          enddo 
!         Map this to pressure grid            
          call fsi_map12(wk_prp,tmp_pr) 
        else    

!         Calculate Grad(U0).eta            
          call opcopy(wk_vxp,wk_vyp,wk_vzp,fsi_grvx0(1,1,1,1,2),
     $                fsi_grvy0(1,1,1,1,2),fsi_grvz0(1,1,1,1,2))
          call opcmult(wk_vxp,wk_vyp,wk_vzp,eta_s) 

!         Calculate Grad(P0).eta            
          call copy(tmp_pr,fsi_grpr0(1,1,1,1,2),lvp)
          call cmult(tmp_pr,eta_s,lvp)
!         Map this to pressure grid            
          call fsi_map12(wk_prp,tmp_pr)
        endif   ! fsi_ifrot 

!       These quantities are treated as perturbation fields   
        call torque_calc_axis_pert(scale,x0,ifdout,iftout,
     $     wk_vx0,wk_vy0,wk_vz0,wk_pr0,wk_vxp,wk_vyp,wk_vzp,wk_prp)
        
        call copy(fdx_dragx (0),dragx (0),maxobj+1)
        call copy(fdx_dragpx(0),dragpx(0),maxobj+1)
        call copy(fdx_dragvx(0),dragpx(0),maxobj+1)
        call copy(fdx_dragy (0),dragy (0),maxobj+1)
        call copy(fdx_dragpy(0),dragpy(0),maxobj+1)
        call copy(fdx_dragvy(0),dragvy(0),maxobj+1)
        call copy(fdx_dragz (0),dragz (0),maxobj+1)
        call copy(fdx_dragpz(0),dragpz(0),maxobj+1)
        call copy(fdx_dragvz(0),dragvz(0),maxobj+1)
        call copy(fdx_torqx (0),torqx (0),maxobj+1)
        call copy(fdx_torqpx(0),torqpx(0),maxobj+1)
        call copy(fdx_torqvx(0),torqpx(0),maxobj+1)
        call copy(fdx_torqy (0),torqy (0),maxobj+1)
        call copy(fdx_torqpy(0),torqpy(0),maxobj+1)
        call copy(fdx_torqvy(0),torqvy(0),maxobj+1)
        call copy(fdx_torqz (0),torqz (0),maxobj+1)
        call copy(fdx_torqpz(0),torqpz(0),maxobj+1)
        call copy(fdx_torqvz(0),torqvz(0),maxobj+1)

      elseif (ifpert.and.fsi_ifdispl.and.npert.eq.2) then

        call opzero(wk_vxp,wk_vyp,wk_vzp)
        call rzero(wk_prp,lvp2)
        scale=fsi_rescale
!       Assuming displacement perturbations are on jp=2            
        call torque_calc_axis_pert(scale,x0,ifdout,iftout,
     $        wk_vxp,wk_vyp,wk_vzp,wk_prp,
     $        vxp(1,2),vyp(1,2),vzp(1,2),prp(1,2))

        call copy(fdx_dragx (0),dragx (0),maxobj+1)
        call copy(fdx_dragpx(0),dragpx(0),maxobj+1)
        call copy(fdx_dragvx(0),dragpx(0),maxobj+1)
        call copy(fdx_dragy (0),dragy (0),maxobj+1)
        call copy(fdx_dragpy(0),dragpy(0),maxobj+1)
        call copy(fdx_dragvy(0),dragvy(0),maxobj+1)
        call copy(fdx_dragz (0),dragz (0),maxobj+1)
        call copy(fdx_dragpz(0),dragpz(0),maxobj+1)
        call copy(fdx_dragvz(0),dragvz(0),maxobj+1)
        call copy(fdx_torqx (0),torqx (0),maxobj+1)
        call copy(fdx_torqpx(0),torqpx(0),maxobj+1)
        call copy(fdx_torqvx(0),torqpx(0),maxobj+1)
        call copy(fdx_torqy (0),torqy (0),maxobj+1)
        call copy(fdx_torqpy(0),torqpy(0),maxobj+1)
        call copy(fdx_torqvy(0),torqvy(0),maxobj+1)
        call copy(fdx_torqz (0),torqz (0),maxobj+1)
        call copy(fdx_torqpz(0),torqpz(0),maxobj+1)
        call copy(fdx_torqvz(0),torqvz(0),maxobj+1)
      else
        call rzero(fdx_dragx (0),maxobj+1)
        call rzero(fdx_dragpx(0),maxobj+1)
        call rzero(fdx_dragvx(0),maxobj+1)
        call rzero(fdx_dragy (0),maxobj+1)
        call rzero(fdx_dragpy(0),maxobj+1)
        call rzero(fdx_dragvy(0),maxobj+1)
        call rzero(fdx_dragz (0),maxobj+1)
        call rzero(fdx_dragpz(0),maxobj+1)
        call rzero(fdx_dragvz(0),maxobj+1)
        call rzero(fdx_torqx (0),maxobj+1)
        call rzero(fdx_torqpx(0),maxobj+1)
        call rzero(fdx_torqvx(0),maxobj+1)
        call rzero(fdx_torqy (0),maxobj+1)
        call rzero(fdx_torqpy(0),maxobj+1)
        call rzero(fdx_torqvy(0),maxobj+1)
        call rzero(fdx_torqz (0),maxobj+1)
        call rzero(fdx_torqpz(0),maxobj+1)
        call rzero(fdx_torqvz(0),maxobj+1)
      endif     ! fsi_ifdispl   

      call fsi_calc_implicit_pert

      fsi_timei=dnekclock()-steptime
      fsi_timea=fsi_timea+fsi_timei

      return
      end subroutine fsi_main_pert_1 
!---------------------------------------------------------------------- 
      subroutine admeshvp
C
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
C
      COMMON /SCRUZ/ FM1(LX1,LY1,LZ1,LELT)
     $             , FM2(LX1,LY1,LZ1,LELT)
     $             , FM3(LX1,LY1,LZ1,LELT)
     $             , PHI(LX1,LY1,LZ1,LELT)
C
      NTOT1=NX1*NY1*NZ1*NELV
C
!      CALL RZERO (FM1,NTOT1)
!      CALL RZERO (FM2,NTOT1)
!      CALL RZERO (FM3,NTOT1)
!C
!      CALL DIVWS (FM1,VXP(1,jp),PHI,NELV,1)
!      CALL DIVWS (FM2,VYP(1,JP),PHI,NELV,2)
!      CALL ADD2  (BFXP(1,jp),FM1,NTOT1)
!      CALL ADD2  (BFYP(1,jp),FM2,NTOT1)
!      IF (NDIM.EQ.3) THEN
!         CALL DIVWS (FM3,VZP(1,jp),PHI,NELV,3)
!         CALL ADD2  (BFZP(1,JP),FM3,NTOT1)
!      ENDIF
C

      CALL RZERO (FM1,NTOT1)
      CALL RZERO (FM2,NTOT1)
      CALL RZERO (FM3,NTOT1)
C
      CALL DIVWS (FM1,VX,PHI,NELV,1)
      CALL DIVWS (FM2,VY,PHI,NELV,2)
      CALL ADD2  (BFXP(1,jp),FM1,NTOT1)
      CALL ADD2  (BFYP(1,jp),FM2,NTOT1)
      IF (NDIM.EQ.3) THEN
         CALL DIVWS (FM3,VZ,PHI,NELV,3)
         CALL ADD2  (BFZP(1,JP),FM3,NTOT1)
      ENDIF


      return
      end
c-----------------------------------------------------------------------

      subroutine incomprp_uzawa_lpert (ux,uy,uz,up)
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

      implicit none
      
      include 'SIZE_DEF'
      include 'SIZE'
      INCLUDE 'INPUT_DEF'
      INCLUDE 'INPUT'               ! IFUZAWA,param(95)
      include 'SOLN_DEF'
      INCLUDE 'SOLN'                ! VTRANS
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'               ! IFIELD,nmxh,tolhv   
      INCLUDE 'CTIMER_DEF'
      include 'CTIMER'

c
      real w1,w2,w3,dv1,dv2,dv3,dp
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)
      real h1,h2,h2inv      
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      real           ux    (lx1,ly1,lz1,lelv)
     $ ,             uy    (lx1,ly1,lz1,lelv)
     $ ,             uz    (lx1,ly1,lz1,lelv)
     $ ,             up    (lx1,ly1,lz1,lelv)

      logical ifprjp
      integer ntot1,ntot2
      integer intloc,intype,istart

      real dtbd

c
!     icalld is declared and saved in CTIMER      
      if (icalld.eq.0) tpres=0.0
      icalld=icalld+1
      npres=icalld
c
      ntot1  = nx1*ny1*nz1*nelv
      ntot2  = nx2*ny2*nz2*nelv

      if (IFUZAWA) then
         intype = -1

         intloc = -1
         call sethlm  (h1,h2,intloc)
         call rzero   (h2inv,ntot1)
      else
         intype = 1
         dtbd   = bd(1)/dt

         call rzero   (h1,ntot1)
         call cmult2  (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)
         call invers2 (h2inv,h2,ntot1)
      endif

      call opdiv   (dp,ux,uy,uz)
      call chsign  (dp,ntot2)
      call ortho   (dp)

C******************************************************************

      ifprjp=.false.    ! project out previous pressure solutions?

      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0.and.(.not.IFUZAWA)) 
     $     ifprjp=.true.

      ! Most likely, the following can be commented out. (pff, 1/6/2010)
!      if (npert.gt.1.or.ifbase)            ifprjp=.false.

      if (ifprjp)   call setrhs_p  (dp,h1,h2,h2inv)
                    call esolver (dp,h1,h2,h2inv,intype)
      if (ifprjp)   call gensolnp_p (dp,h1,h2,h2inv)

C******************************************************************

      call opgradt (w1 ,w2 ,w3 ,dp)
      if (IFUZAWA) then
         call ophinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h1 ,h2 ,tolhv ,nmxh)
      else
         call opbinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
      endif
      call opadd2  (ux ,uy ,uz ,dv1,dv2,dv3)
c
!     prabal. 
!     I've changed the pressure extrapolation to be consistent
!     between non-linear and linear cases.
!     Hence the pressure correction here is also changed.      
      call add2(up,dp,ntot2)
c
      return
      end subroutine incomprp_uzawa_lpert
c------------------------------------------------------------------------
      subroutine setrhs_p(p,h1,h2,h2inv)
C
C     Project rhs onto best fit in the "E" norm.
C
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'MASS_DEF'
      include 'MASS'
      include 'TSTEP_DEF'
      include 'TSTEP'

!     Pointer to perturbation arrays.
!     Included in SOLN      
      integer jp
      common /ppointr/ jp

      REAL             P    (LX2,LY2,LZ2,LELV)
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
C
      logical ifdump
      save    ifdump
      data    ifdump /.false./
C

      REAL RHS,Pbar,Pnew,Pbrr
      REAL ALPHA,WORK,ALPHAN,DTLAST
      INTEGER Nprev,Mprev

      INTEGER LTOT2
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)

      INTEGER LPTOT2
      PARAMETER (LPTOT2=LPX2*LPY2*LPZ2*LPELV)
!     These two common blocks need to extended for multiple perturbation
!     mode case (NPERT>1). Or when both pert and base flow evolution are
!     being evaluated. Therefore new common blocks have been defined. 
      COMMON /ORTHOVP/ RHS(LPTOT2,MXPREV,lpert)
      COMMON /ORTHOIP/ Nprev(lpert),Mprev(lpert)

!     Following two common blocks are the same as the ones used in the
!     Std. Nek case (non-linear).      
      COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2),
     $                Pbrr(LTOT2)
      COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), 
     $                ALPHAN, DTLAST

      real vlsc2
      
      real alphas
      integer NTOT2
      integer I,ierr,INTETYPE

      integer icalld
      save    icalld
      data    icalld/0/

      integer icalld_pert(lpert)
      save icalld_pert

C
C     First call, we have no vectors to orthogonalize against.
      if (icalld.eq.0) then
        call izero(icalld_pert,lpert)
        icalld=icalld+1
      endif  

      if (icalld_pert(jp).eq.0) then
         Nprev(jp)=0
         Mprev(jp)=param(93)
         Mprev(jp)=min(Mprev(jp),Mxprev)
         icalld_pert(jp)=icalld_pert(jp)+1
      endif
C
      NTOT2  = LPTOT2
C
C     Update rhs's if E-matrix has changed
C
      CALL UPDRHSE_P(P,H1,H2,H2INV,ierr)
      if (ierr.eq.1) Nprev(jp)=0
C
C     Perform Gram-Schmidt for previous rhs's.
C
      DO 10 I=1,Nprev(jp)
         ALPHA(i) = VLSC2(P,RHS(1,i,jp),NTOT2)
   10 CONTINUE
C
      IF (Nprev(jp).GT.0) CALL gop(alpha,WORK,'+  ',Nprev)
C
      CALL RZERO(Pbar,NTOT2)
      DO 20 I=1,Nprev(jp)
         alphas = alpha(i)
         CALL ADD2S2(Pbar,RHS(1,i,jp),alphas,NTOT2)
   20 CONTINUE
C
      if (Nprev(jp).gt.0) then
         INTETYPE = 1
         CALL CDABDTP(Pnew,Pbar,H1,H2,H2INV,INTETYPE)
         CALL SUB2   (P,Pnew,NTOT2)
      endif
C
      return
      end subroutine setrhs_p
c-----------------------------------------------------------------------

      subroutine updrhse_p(p,h1,h2,h2inv,ierr)
C
C     Update rhs's if E-matrix has changed
C
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'      
      include 'INPUT'
      include 'MASS_DEF'
      include 'MASS'
      include 'TSTEP_DEF'
      include 'TSTEP'

!     Pointer to perturbation arrays.
!     Included in SOLN      
      integer jp
      common /ppointr/ jp

      integer LTOT2
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)

      LOGICAL IFNEWE
      COMMON /ORTHOLP/ IFNEWE

      REAL ALPHA,WORK,ALPHAN,DTLAST
      COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), 
     $                ALPHAN, DTLAST

      INTEGER Nprev,Mprev
      COMMON /ORTHOIP/ Nprev(lpert),Mprev(lpert)

      REAL             P    (LX2,LY2,LZ2,LELV)
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
C
      integer icalld
      save    icalld
      data    icalld/0/

      integer NTOT2
      integer Iprev,Nprevt,ierr

!     From FSI include file.
!     Can't include the whole file because of some variable name
!     conflicts.
!     Just including the logical block
      logical iffsi           ! if we are doing fsi simulations
      logical fsi_ifinit      ! if fsi has been initialized
      logical fsi_iftermso    ! if output terms for debugging
      logical fsi_pert_ifmvbd ! if move boundaries in fsi perturbation mode
      Common /FSI_SOLVEl/ iffsi,fsi_ifinit,fsi_iftermso,
     $                    fsi_pert_ifmvbd


      NTOT2=NX2*NY2*NZ2*NELV
C
C     First, we have to decide if the E matrix has changed.
C
      IF (icalld.eq.0) THEN
         icalld=1
         DTlast=DT
      ENDIF
C
      IFNEWE=.FALSE.
      IF (fsi_pert_ifmvbd) THEN
         IFNEWE=.TRUE.
         CALL INVERS2(bm2inv,bm2,Ntot2)
      ELSEIF (DTlast.ne.DT) THEN
         IFNEWE=.TRUE.
         DTlast=DT
      ENDIF
      IF (IFNEWE.and.nio.eq.0) write(6,*) 'reorthogo,jp:',nprev(jp),jp
C     
C     Next, we reconstruct a new rhs set.
C     
      IF (IFNEWE) THEN
c
c        new idea...
c        if (nprev.gt.0) nprev=1
c        call copy(rhs,pnew,ntot2)
c
         Nprevt = Nprev(jp)
         DO 100 Iprev=1,Nprevt
C           Orthogonalize this rhs w.r.t. previous rhs's
            CALL ECONJP_P (Iprev,H1,H2,H2INV,ierr)
            if (ierr.eq.1) then
               Nprev(jp) = 0
               return
            endif
  100    CONTINUE
C
      ENDIF
C
      RETURN
      end subroutine updrhse_p
c-----------------------------------------------------------------------
      subroutine econjp_p(kprev,h1,h2,h2inv,ierr)
C
C     Orthogonalize the rhs wrt previous rhs's for which we already
C     know the soln.
C
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'      
      include 'INPUT'
      include 'MASS_DEF'
      include 'MASS'
      include 'TSTEP_DEF'
      include 'TSTEP'

!     Pointer to perturbation arrays.
!     Included in SOLN      
      integer jp
      common /ppointr/ jp

C
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)

      INTEGER LTOT2
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
      REAL RHS,Pbar,Pnew,Pbrr
      REAL ALPHA,WORK,ALPHAN,DTLAST
      INTEGER Nprev,Mprev

      INTEGER LPTOT2
      PARAMETER (LPTOT2=LPX2*LPY2*LPZ2*LPELV)
!     These two common blocks need to extended for multiple perturbation
!     mode case (NPERT>1). Or when both pert and base flow evolution are
!     being evaluated. Therefore new common blocks have been defined. 
      COMMON /ORTHOVP/ RHS(LPTOT2,MXPREV,lpert)
      COMMON /ORTHOIP/ Nprev(lpert),Mprev(lpert)

!     Following two common blocks are the same as the ones used in the
!     Std. Nek case (non-linear).      
      COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2),
     $                Pbrr(LTOT2)
      COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), 
     $                ALPHAN, DTLAST

      integer ierr,NTOT2,INTETYPE
      integer ipass,npass
      real Alphad,alpham

      real glsc2        ! function
      real vlsc2        ! function

      integer i,Kprev,Kprev1

      ierr  = 0
      NTOT2 = NX2*NY2*NZ2*NELV
      INTETYPE=1
C
C     Gram Schmidt, w re-orthogonalization
C
      npass=1
      if (abs(param(105)).eq.2) npass=2
      do ipass=1,npass
c
         CALL CDABDTP(Pbrr,RHS(1,Kprev,jp),H1,H2,H2INV,INTETYPE)
C
C        Compute part of the norm
         Alphad = GLSC2(RHS(1,Kprev,jp),Pbrr,NTOT2)
C
C        Gram-Schmidt
         Kprev1=Kprev-1
         DO 10 I=1,Kprev1
            ALPHA(I) = VLSC2(Pbrr,RHS(1,i,jp),NTOT2)
   10    CONTINUE
         IF (Kprev1.GT.0) CALL gop(alpha,WORK,'+  ',Kprev1)
C
         DO 20 I=1,Kprev1
            alpham = -alpha(i)
            CALL ADD2S2(RHS(1,Kprev,jp),RHS(1,i,jp),alpham,NTOT2)
            Alphad = Alphad - alpha(i)**2
   20    CONTINUE
      enddo
C
C    .Normalize new element in P~
C
      if (ALPHAd.le.0.0) then
         if (nio.eq.0) write(6,*) 
     $   'ERROR:  alphad .le. 0 in ECONJ',alphad,Kprev
         ierr = 1
         return
      endif
      ALPHAd = 1.0/SQRT(ALPHAd)
      ALPHAN = Alphad
      CALL CMULT(RHS(1,Kprev,jp),alphan,NTOT2)
C
      return
      end subroutine econjp_p
c-----------------------------------------------------------------------
      subroutine gensolnp_p(p,h1,h2,h2inv)
C
C     Reconstruct the solution to the original problem by adding back
C     the previous solutions
C     know the soln.
C

      implicit none

      include 'SIZE_DEF'            
      include 'SIZE'

!     Pointer to perturbation arrays.
!     Included in SOLN      
      integer jp
      common /ppointr/ jp

      INTEGER LTOT2
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)

      REAL RHS,Pbar,Pnew,Pbrr
      REAL ALPHA,WORK,ALPHAN,DTLAST
      INTEGER Nprev,Mprev

      INTEGER LPTOT2
      PARAMETER (LPTOT2=LPX2*LPY2*LPZ2*LPELV)
!     These two common blocks need to extended for multiple perturbation
!     mode case (NPERT>1). Or when both pert and base flow evolution are
!     being evaluated. Therefore new common blocks have been defined. 
      COMMON /ORTHOVP/ RHS(LPTOT2,MXPREV,lpert)
      COMMON /ORTHOIP/ Nprev(lpert),Mprev(lpert)

!     Following two common blocks are the same as the ones used in the
!     Std. Nek case (non-linear).      
      COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2),
     $                Pbrr(LTOT2)
      COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), 
     $                ALPHAN, DTLAST

      REAL             P    (LX2,LY2,LZ2,LELV)
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)

      integer NTOT2,ierr

C
      NTOT2=LPTOT2
C
C     First, save current solution
C
      CALL COPY (Pnew,P,NTOT2)
C
C     Reconstruct solution
C
      CALL ADD2(P,Pbar,NTOT2)
C
C     Update the set of <p,rhs>
C
      CALL UPDTSETP_P(P,H1,H2,H2INV,ierr)
      if (ierr.eq.1) Nprev(jp) = 0
c
      return
      end subroutine gensolnp_p
c-----------------------------------------------------------------------
      subroutine updtsetp_p(p,h1,h2,h2inv,IERR)
C
C     Update the set of rhs's and the corresponding p-set:
C
C        . Standard case is to add P_new, and RHS_new = E*P_new
C
C        . However, when Nprev=Mprev (max. allowed), we throw out
C          the old set, and set P_1 = P, RHS_1=E*P_1
C
C        . Other schemes are possible, e.g., let's save a bunch of
C          old vectors, perhaps chosen wisely via P.O.D.
C
C
      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
!      include 'MASS'

!     Pointer to perturbation arrays.
!     Included in SOLN      
      integer jp
      common /ppointr/ jp

      INTEGER LTOT2
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)

      REAL RHS,Pbar,Pnew,Pbrr
      REAL ALPHA,WORK,ALPHAN,DTLAST
      INTEGER Nprev,Mprev

      INTEGER LPTOT2
      PARAMETER (LPTOT2=LPX2*LPY2*LPZ2*LPELV)
!     These two common blocks need to extended for multiple perturbation
!     mode case (NPERT>1). Or when both pert and base flow evolution are
!     being evaluated.      
      COMMON /ORTHOVP/ RHS(LPTOT2,MXPREV,lpert)
      COMMON /ORTHOIP/ Nprev(lpert),Mprev(lpert)

!     Following two common blocks are the same as the ones used in the
!     Std. Nek case (non-linear).      
      COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2),
     $                Pbrr(LTOT2)
      COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), 
     $                ALPHAN, DTLAST

      REAL             P    (LX2,LY2,LZ2,LELV)
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)

      INTEGER NTOT2,ierr

C
      NTOT2=LPTOT2
C
      IF (Nprev(jp).EQ.Mprev(jp)) THEN
         CALL COPY(Pnew,P,NTOT2)
         Nprev(jp)=0
      ENDIF
C
C     Increment solution set
      Nprev(jp) = Nprev(jp)+1
C
      CALL COPY   (RHS(1,Nprev(jp),jp),Pnew,NTOT2)
C
C     Orthogonalize rhs against previous rhs and normalize
C
      CALL ECONJP_P(Nprev(jp),H1,H2,H2INV,ierr)
C
c     Save last sol'n
      CALL COPY(Pnew,P,NTOT2)
C
      return
      end subroutine updtsetp_p
c-----------------------------------------------------------------------
      subroutine fsi_meshv_pert1

!     This routine is called at the end of perturbation computations.
!     To update the wall/mesh velocities for the Base flow solve.            

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'          ! xm1,ym1,zm1
      include 'MVGEOM_DEF'
      include 'MVGEOM'        ! wx,wy,wz
      include 'WING_MVMSH'
      include 'FSI'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'

      real ucx,ucy,ucz        ! mesh/boundary velocities
      real dx,dy,dz

      integer i,ix,e,n

      n = nx1*ny1*nz1

      do e=1,nelv

        do i=1,n

          ix = (e-1)*n + i

          if (.not.fsi_ifrot) then

!           Set Mesh velocities
            if (npert.gt.1) then                
              if (jp.eq.1) then
!               This sets mesh velocities for jp=2 solve                  
                ucx =  0.
                ucy =  etav
                ucz =  0.
              elseif (jp.eq.2) then
!               This sets up mesh velocities for jp=1 solve 
                ucx =  0.
                ucy =  0.
                ucz =  0.
              else
                if (nio.eq.0) write(6,*) 
     $               'Unknown Mesh velocity conditions for jp=',jp
                call exitt
              endif
            else
!             Assuming we don't want to decouple the two solutions
              ucx =  0.
              ucy =  etav
              ucz =  0.
            endif  
            wx(i,1,1,e) = basev(ix)*ucx
            wy(i,1,1,e) = basev(ix)*ucy
            if (ndim.eq.3) wz(i,1,1,e) = basev(ix)*ucz

!           Set Boundary conditions
            if (npert.gt.1) then 
              if (jp.eq.1) then
!               This sets Boundary conditions for jp=2 solve                  
                ucx =  0.
                ucy =  0.
                ucz =  0.
              elseif (jp.eq.2) then
!               This sets Boundary conditions for jp=1 solve 
                ucx =  0.
                ucy =  etav_s
                ucz =  0.
              else
                if (nio.eq.0) write(6,*) 
     $               'Unknown Boundary Conditions for jp=',jp
                call exitt
              endif
            else 
!             Assuming we don't want to decouple the two solutions
              ucx =  0.
              ucy =  etav_s
              ucz =  0.
            endif
            umeshx(i,1,1,e) = basev(ix)*ucx
            umeshy(i,1,1,e) = basev(ix)*ucy
            if (ndim.eq.3) umeshz(i,1,1,e) = basev(ix)*ucz

          else
!           Rotational FSI      
            dx = xm1(i,1,1,e) - fsi_x0
            dy = ym1(i,1,1,e) - fsi_y0

!           Set Mesh velocities
            if (npert.gt.1) then                
              if (jp.eq.1) then
!               This sets mesh velocities for jp=2 solve                  
                ucx =  etav*dy
                ucy = -etav*dx
                ucz =  0.0
              elseif (jp.eq.2) then
!               This sets up mesh velocities for jp=1 solve 
                ucx =  0.
                ucy =  0.
                ucz =  0.
              else
                if (nio.eq.0) write(6,*) 
     $               'Unknown Mesh velocity conditions for jp=',jp
                call exitt
              endif
            else
!             Assuming we don't want to decouple the two solutions
              ucx =  etav*dy
              ucy = -etav*dx
              ucz =  0.0
            endif 
            wx(i,1,1,e) = basev(ix)*ucx
            wy(i,1,1,e) = basev(ix)*ucy
            if (ndim.eq.3) wz(i,1,1,e) = basev(ix)*ucz

!           Set Boundary conditions
            if (npert.gt.1) then                
              if (jp.eq.1) then
!               This sets Boundary conditions for jp=2 solve                  
                ucx =  0.
                ucy =  0.
                ucz =  0.
              elseif (jp.eq.2) then
!               This sets Boundary Conditions for jp=1 solve 
                ucx =  etav_s*dy
                ucy = -etav_s*dx
                ucz =  0.0
              else
                if (nio.eq.0) write(6,*) 
     $               'Unknown Boundary conditions for jp=',jp
                call exitt
              endif
            else
!             Assuming we don't want to decouple the two solutions
              ucx =  etav_s*dy
              ucy = -etav_s*dx
              ucz =  0.0
            endif  
            umeshx(i,1,1,e) = basev(ix)*ucx
            umeshy(i,1,1,e) = basev(ix)*ucy
            if (ndim.eq.3) umeshz(i,1,1,e) = basev(ix)*ucz

          endif         ! fsi_ifrot

        enddo     ! i=1,n 
      enddo       ! e=1,nelv

      return
      end subroutine fsi_meshv_pert1
!----------------------------------------------------------------------
      subroutine fsi_meshv_pert2

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'          ! xm1,ym1,zm1
      include 'MVGEOM_DEF'
      include 'MVGEOM'        ! wx,wy,wz
      include 'WING_MVMSH'
      include 'FSI'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'

      real ucx,ucy,ucz        ! mesh velocities
      real dx,dy,dz
      real eta_sx,eta_sy,eta_sz

      integer i,n
      real rot_sn(3,3)
      real rot_cs(3,3)
      real dxyz(3),r0(3)


      call rzero(rot_sn,9)            ! sine terms
      call rzero(rot_cs,9)            ! cosine terms
!     Matrix of sine terms
!     For a clockwise rotation for a positive eta
      if (fsi_ifrot) then
        rot_sn(1,2) =  eta_s
        rot_sn(2,1) = -eta_s
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
        if (.not.fsi_ifrot) then

!         Mesh motion                
          if (fsi_ifale) then  
            ucx =  0.
            ucy =  etav
            ucz =  0.
          else
            ucx =  0.
            ucy =  0. 
            ucz =  0.
          endif  

          wx(i,1,1,1) = basev(i)*ucx
          wy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) wz(i,1,1,1) = basev(i)*ucz

!         Boundary condition            
          eta_sx = 0.
          eta_sy = eta_s
          eta_sz = 0.  

          if (fsi_bcdispl) then
            
            fsi_bcx(i,1,1,1) = (fsi_grvx0(i,1,1,1,1)*eta_sx + 
     $             fsi_grvx0(i,1,1,1,2)*eta_sy +
     $             fsi_grvx0(i,1,1,1,3)*eta_sz)
            fsi_bcy(i,1,1,1) =  (fsi_grvy0(i,1,1,1,1)*eta_sx + 
     $             fsi_grvy0(i,1,1,1,2)*eta_sy +
     $             fsi_grvy0(i,1,1,1,3)*eta_sz)
            
            ucx = 0. - fsi_bcx(i,1,1,1)
            ucy = etav_s - fsi_bcy(i,1,1,1)
            if (if3d) then
              fsi_bcz(i,1,1,1) =  (fsi_grvz0(i,1,1,1,1)*eta_sx + 
     $               fsi_grvz0(i,1,1,1,2)*eta_sy +
     $               fsi_grvz0(i,1,1,1,3)*eta_sz)
              ucz = 0. - fsi_bcz(i,1,1,1)
            else
              ucz = 0.
            endif 
          else
            ucx = 0.
            ucy = etav_s
            ucz = 0.
          endif          

          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz

        else
!         Rotational FSI

!         Distance from Rotational axis                
          dx = xm1(i,1,1,1) - fsi_x0
          dy = ym1(i,1,1,1) - fsi_y0
          dz = 0.

!         Mesh motion 
          if (fsi_ifale) then  
            ucx =  etav*dy
            ucy =  -etav*dx
            ucz =  0.
          else
            ucx =  0.
            ucy =  0.
            ucz =  0.
          endif  
  
          wx(i,1,1,1) = basev(i)*ucx
          wy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) wz(i,1,1,1) = basev(i)*ucz

          if (fsi_bcdispl) then

!           Distance from axis   
            r0(1) = dx 
            r0(2) = dy
            r0(3) = dz

!           perturbation in position
            call mxm(rot_sn,3,r0,3,dxyz,1)
            eta_sx = dxyz(1)
            eta_sy = dxyz(2)
            eta_sz = dxyz(3) 

            fsi_bcx(i,1,1,1) = (fsi_grvx0(i,1,1,1,1)*eta_sx + 
     $             fsi_grvx0(i,1,1,1,2)*eta_sy +
     $             fsi_grvx0(i,1,1,1,3)*eta_sz) 
            fsi_bcy(i,1,1,1) = (fsi_grvy0(i,1,1,1,1)*eta_sx + 
     $             fsi_grvy0(i,1,1,1,2)*eta_sy +
     $             fsi_grvy0(i,1,1,1,3)*eta_sz)

            ucx = etav_s*dy - fsi_bcx(i,1,1,1)
            ucy = -etav_s*dx - fsi_bcy(i,1,1,1)
            if (if3d) then
              fsi_bcz(i,1,1,1) = (fsi_grvz0(i,1,1,1,1)*eta_sx + 
     $               fsi_grvz0(i,1,1,1,2)*eta_sy +
     $               fsi_grvz0(i,1,1,1,3)*eta_sz)
              ucz = 0. - fsi_bcz(i,1,1,1)
            else
              ucz = 0.
            endif 
          else          ! not using displacement effects in BCs
            ucx =  etav_s*dy
            ucy =  -etav_s*dx
            ucz =  0.
          endif          
         
          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz
        endif
       
      enddo

      return
      end subroutine fsi_meshv_pert2
!---------------------------------------------------------------------- 
      subroutine torque_calc_axis_pert(scale,x0,ifdout,iftout,velx,
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
      INCLUDE 'FSI' 

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
      real pm0(lt)
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
      call mappr(pm0,press,xm0,ym0)       ! xm0,ym0 used as work arrays
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

      call add2s2(pm0,xm1,dpdx_mean,n)  ! Doesn't work if object is cut by 
      call add2s2(pm0,ym1,dpdy_mean,n)  ! periodicboundary.  In this case,
      call add2s2(pm0,zm1,dpdz_mean,n)  ! set ._mean=0 and compensate in
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
                  call clcdcm_pert(dgtq,dgtq_p,sa,xm0,ym0,zm0,sij,sij0,
     $                             pm1,pm0,vdiff,ifc,ie)

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
      end subroutine torque_calc_axis_pert
!----------------------------------------------------------------------
      subroutine clcdcm_pert(dgtq,dgtq_p,a,xm0,ym0,zm0,sij,sij0,
     $                       pm1,pm0,visc,f,e)

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
      INCLUDE 'FSI'
c
      real dgtq(3,4)                      ! Total everything
      real dgtq_p(3,8)                    ! Saving all individual terms
      real xm0 (lx1,ly1,lz1,lelt)
      real ym0 (lx1,ly1,lz1,lelt)
      real zm0 (lx1,ly1,lz1,lelt)
      real sij (lx1,ly1,lz1,3*ldim-3,lelv)! Perturbation shear stress
      real sij0(lx1,ly1,lz1,3*ldim-3,lelv)! Mean Shear stress
      real pm1 (lx1,ly1,lz1,lelv)         ! perturbation pressure mapped
      real pm0 (lx1,ly1,lz1,lelv)         ! mean pressure mapped
      real visc(lx1,ly1,lz1,lelv)

      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1)
c

      real dg(3,7) !(p.n0),(Tauij.n0),(Tauij.[Sn].n0),
                   !(tauij.n0),(P.[Sn].n0),(P.n0),(tauij.[Sn].n0)
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
      if (ifpert.and.iffsi) then
!        continue
!          write(6,*) 'clcdcm_pert'
      else
        if (nio.eq.0) then
          write(6,*) 'clcdcm_pert called without proper conditions'
          write(6,*) 'IFPERT/IFFSI',ifpert,iffsi
        endif
        call exitt
      endif      

!     In this case cosider rotation about 'Z' axis
!     Rewrite rotation matrix as sum of two matrices.
!     One with just the cosine terms and one with just the sine terms.      
!     Rot = Cs + Sn
!     With linearization: 
!           eta << 1      
!           sin(eta) ~ eta
!           cos(eta) ~ 1
!     The Cs becomes an Identity matrix.
!     The Sn only has eta terms.
!     Sn can therefore be considered a 'perturbation matrix'

      call rzero(rot_sn,9)            ! sine terms
      call rzero(rot_cs,9)            ! cosine terms
!     Matrix of sine terms
!     For a clockwise rotation for a positive eta
      if (.not.fsi_ifrot) then
!       If there is only vertical motion,
!       Then the direction of normals does not change
!       Also the moment arm of points on the surface does not change
!       As long as it is a rigid body motion. 
        continue
      else
        rot_sn(1,2) = eta_s
        rot_sn(2,1) = -eta_s
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

!        Original normals   
         n1 = n0(1)
         n2 = n0(2)
         n3 = n0(3)

!        Calculate rotated normals            
         call mxm(rot_sn,3,n0,3,nr,1)
         n1p = nr(1)
         n2p = nr(2)
         n3p = nr(3)

!        Calculate distance from axis   
         r0(1) = xm0(j1,j2,1,e)
         r0(2) = ym0(j1,j2,1,e)
         r0(3) = zm0(j1,j2,1,e)

!        Original distance
         r1 = r0(1)
         r2 = r0(2)
         r3 = r0(3)

!        Distance perturbation 
         call mxm(rot_sn,3,r0,3,rr,1)
         r1p = rr(1)
         r2p = rr(2)
         r3p = rr(3)

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
         dg(1,5) = pm0(j1,j2,1,e)*n1p 
         dg(2,5) = pm0(j1,j2,1,e)*n2p
         dg(3,5) = pm0(j1,j2,1,e)*n3p

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

!        Original normals   
         n1 = n0(1)
         n2 = n0(2)
         n3 = n0(3)

!        Calculate rotated normals            
         call mxm(rot_sn,3,n0,3,nr,1)
         n1p = nr(1)
         n2p = nr(2)
         n3p = 0.

!        prabal
!         write(6,12) 'Rotation:',istep,nid,e,n1,n2,eta_s,n1p,n2p
!   12 format(A9,1x,i3,1x,i4,1x,i4,1x,5(1pE14.7E2,1x))     

!        Calculate distance from axis   
         r0(1) = xm0(j1,j2,1,e)
         r0(2) = ym0(j1,j2,1,e)
         r0(3) = 0.

!        original distance
         r1 = r0(1)
         r2 = r0(2)
         r3 = 0.

!        Distance perturbation 
         call mxm(rot_sn,3,r0,3,rr,1)
         r1p = rr(1)
         r2p = rr(2)
         r3p = 0.

!        prabal
!         write(6,12) 'Distance:',istep,nid,e,r1,r2,eta_s,r1p,r2p

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

!        (Mean Stress)*(Normals perturbation)            
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

!        (Mean pressure)*(Perturbation normals)
!        P.[Sn].n0        
         dg(1,5) = pm0(j1,j2,1,e)*n1p
         dg(2,5) = pm0(j1,j2,1,e)*n2p
         dg(3,5) = 0.

!        (Mean pressure)*(Original normals)
!        P.n0        
         dg(1,6) = pm0(j1,j2,1,e)*n1
         dg(2,6) = pm0(j1,j2,1,e)*n2
         dg(3,6) = 0.

!        (Perturbation Stress)*(Perturbation normals) ! Non-linear term           
!        tauij*[Sn]*n0            
         dg(1,7) = -v*(s11*n1p + s12*n2p + s13*n3p)
         dg(2,7) = -v*(s21*n1p + s22*n2p + s23*n3p)
         dg(3,7) = 0.

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
!        (Original Distance) x [(Perturbation pressure)*(Original normals)]
!        [r0] x (p.n0)            
         dgtq(1,3) = 0. ! dgtq(1,3) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq(2,3) = 0. ! dgtq(2,3) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq(3,3) = dgtq(3,3) + (r1*dg(2,1)-r2*dg(1,1))

!        Saving individual components            
         dgtq_p(1,5) = 0. ! dgtq_p(1,5) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq_p(2,5) = 0. ! dgtq_p(2,5) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq_p(3,5) = dgtq_p(3,5) + (r1*dg(2,1)-r2*dg(1,1))

!        [[Sn]*r0] x (P.n0) and  [r0] x (P.[Sn]*n0) Should both be zero
!        for pure rotation. Adding them anyway

!        Mean pressure torque.            
!        (Perturbation Distance) x [(Mean pressure)*(Original normals)]
!        [[Sn]*r0] x (P.n0)            
         dgtq(1,3) = 0. ! dgtq(1,3) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq(2,3) = 0. ! dgtq(2,3) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq(3,3) = dgtq(3,3) + (r1p*dg(2,6)-r2p*dg(1,6))

!        Saving individual components            
         dgtq_p(1,5) = 0. ! dgtq_p(1,5) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq_p(2,5) = 0. ! dgtq_p(2,5) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq_p(3,5) = dgtq_p(3,5) + (r1p*dg(2,1)-r2p*dg(1,1))

!        Mean pressure torque.            
!        (Original Distance) x [(Mean pressure)*(Perturbation normals)]
!        [r0] x (P.[Sn]*n0)            
         dgtq(1,3) = 0. ! dgtq(1,3) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq(2,3) = 0. ! dgtq(2,3) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq(3,3) = dgtq(3,3) + (r1*dg(2,5)-r2*dg(1,5))

!        Saving individual components            
         dgtq_p(1,5) = 0. ! dgtq_p(1,5) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq_p(2,5) = 0. ! dgtq_p(2,5) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq_p(3,5) = dgtq_p(3,5) + (r1*dg(2,5)-r2*dg(1,5))

!!       Viscous torque
!        (Perturbation Distance) x [(Mean Stress)*(Original normals)]
!        [Sn*r0] x (Tauij*n0)            
         dgtq(1,4) = 0. ! dgtq(1,4) + (r2p*dg(3,2)-r3p*dg(2,2))
         dgtq(2,4) = 0. ! dgtq(2,4) + (r3p*dg(1,2)-r1p*dg(3,2))
         dgtq(3,4) = dgtq(3,4) + (r1p*dg(2,2)-r2p*dg(1,2))

!        Saving individual components            
         dgtq_p(1,6) = 0. ! dgtq_p(1,6) + (r2p*dg(3,2)-r3p*dg(2,2))
         dgtq_p(2,6) = 0. ! dgtq_p(2,6) + (r3p*dg(1,2)-r1p*dg(3,2))
         dgtq_p(3,6) = dgtq_p(3,6) + (r1p*dg(2,2)-r2p*dg(1,2))

!        (Original Distance) x [(Mean Stress)*(Perturbation normals)]
!        (r0) x (Tauij*[Sn]*n0)            
         dgtq(1,4) = 0. ! dgtq(1,4) + (r2*dg(3,3)-r3*dg(2,3))
         dgtq(2,4) = 0. ! dgtq(2,4) + (r3*dg(1,3)-r1*dg(3,3))
         dgtq(3,4) = dgtq(3,4) + (r1*dg(2,3)-r2*dg(1,3))

!        Saving individual components            
         dgtq_p(1,7) = 0. ! dgtq_p(1,7) + (r2*dg(3,3)-r3*dg(2,3))
         dgtq_p(2,7) = 0. ! dgtq_p(2,7) + (r3*dg(1,3)-r1*dg(3,3))
         dgtq_p(3,7) = dgtq_p(3,7) + (r1*dg(2,3)-r2*dg(1,3))

!        (Original Distance) x [(Perturbation Stress)*(Original normals)]
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
      end subroutine clcdcm_pert

!-----------------------------------------------------------------------
      subroutine clcdcm_pert2(dgtq,dgtq_p,a,xm0,ym0,zm0,sij,sij0,
     $                       pm1,pm0,visc,f,e)

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
      INCLUDE 'FSI'
c
      real dgtq(3,4)                      ! Total everything
      real dgtq_p(3,8)                    ! Saving all individual terms
      real xm0 (lx1,ly1,lz1,lelt)
      real ym0 (lx1,ly1,lz1,lelt)
      real zm0 (lx1,ly1,lz1,lelt)
      real sij (lx1,ly1,lz1,3*ldim-3,lelv)! Perturbation shear stress
      real sij0(lx1,ly1,lz1,3*ldim-3,lelv)! Mean Shear stress
      real pm1 (lx1,ly1,lz1,lelv)         ! perturbation pressure mapped
      real pm0 (lx1,ly1,lz1,lelv)         ! mean pressure mapped
      real visc(lx1,ly1,lz1,lelv)

      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1)
c

      real dg(3,7) !(p.n0),(Tauij.n0),(Tauij.[Sn].n0),
                   !(tauij.n0),(P.[Sn].n0),(P.n0),(tauij.[Sn].n0)
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
      if (ifpert.and.iffsi) then
!        continue
!          write(6,*) 'clcdcm_pert'
      else
        if (nio.eq.0) then
          write(6,*) 'clcdcm_pert called without proper conditions'
          write(6,*) 'IFPERT/IFFSI',ifpert,iffsi
        endif
        call exitt
      endif      

!     In this case cosider rotation about 'Z' axis
!     Rewrite rotation matrix as sum of two matrices.
!     One with just the cosine terms and one with just the sine terms.      
!     Rot = Cs + Sn
!     With linearization: 
!           eta << 1      
!           sin(eta) ~ eta
!           cos(eta) ~ 1
!     The Cs becomes an Identity matrix.
!     The Sn only has eta terms.
!     Sn can therefore be considered a 'perturbation matrix'

      call rzero(rot_sn,9)            ! sine terms
      call rzero(rot_cs,9)            ! cosine terms
!     Matrix of sine terms
!     For a clockwise rotation for a positive eta
      if (.not.fsi_ifrot) then
!       If there is only vertical motion,
!       Then the direction of normals does not change
!       Also the moment arm of points on the surface does not change
!       As long as it is a rigid body motion. 
        continue
      else
        rot_sn(1,2) = eta_s
        rot_sn(2,1) = -eta_s
!       prabal
        rot_sn(1,1) = 1
        rot_sn(2,2) = 1
        rot_sn(3,3) = 1
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

!        Original normals   
         n1 = n0(1)
         n2 = n0(2)
         n3 = n0(3)

!        Calculate rotated normals            
         call mxm(rot_sn,3,n0,3,nr,1)
         n1p = nr(1)
         n2p = nr(2)
         n3p = nr(3)

!        Calculate distance from axis   
         r0(1) = xm0(j1,j2,1,e)
         r0(2) = ym0(j1,j2,1,e)
         r0(3) = zm0(j1,j2,1,e)

!        Original distance
         r1 = r0(1)
         r2 = r0(2)
         r3 = r0(3)

!        Distance perturbation 
         call mxm(rot_sn,3,r0,3,rr,1)
         r1p = rr(1)
         r2p = rr(2)
         r3p = rr(3)

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
         dg(1,5) = pm0(j1,j2,1,e)*n1p 
         dg(2,5) = pm0(j1,j2,1,e)*n2p
         dg(3,5) = pm0(j1,j2,1,e)*n3p

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

!        Original normals   
         n1 = n0(1)
         n2 = n0(2)
         n3 = n0(3)

!        Calculate rotated normals            
         call mxm(rot_sn,3,n0,3,nr,1)
         n1p = nr(1)
         n2p = nr(2)
         n3p = 0.

!        prabal
         n1 = n1p
         n2 = n2p
         n3 = n3p

!        prabal
!         write(6,12) 'Rotation:',istep,nid,e,n1,n2,eta_s,n1p,n2p
!   12 format(A9,1x,i3,1x,i4,1x,i4,1x,5(1pE14.7E2,1x))     

!        Calculate distance from axis   
         r0(1) = xm0(j1,j2,1,e)
         r0(2) = ym0(j1,j2,1,e)
         r0(3) = 0.

!        original distance
         r1 = r0(1)
         r2 = r0(2)
         r3 = 0.

!        Distance perturbation 
         call mxm(rot_sn,3,r0,3,rr,1)
         r1p = rr(1)
         r2p = rr(2)
         r3p = 0.

!        prabal
         r1 = r1p
         r2 = r2p
         r3 = r3p

!        prabal
!         write(6,12) 'Distance:',istep,nid,e,r1,r2,eta_s,r1p,r2p

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

!        (Mean Stress)*(Normals perturbation)            
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

!        (Mean pressure)*(Perturbation normals)
!        P.[Sn].n0        
         dg(1,5) = pm0(j1,j2,1,e)*n1p
         dg(2,5) = pm0(j1,j2,1,e)*n2p
         dg(3,5) = 0.

!        (Mean pressure)*(Original normals)
!        P.n0        
         dg(1,6) = pm0(j1,j2,1,e)*n1
         dg(2,6) = pm0(j1,j2,1,e)*n2
         dg(3,6) = 0.

!        (Perturbation Stress)*(Perturbation normals) ! Non-linear term           
!        tauij*[Sn]*n0            
         dg(1,7) = -v*(s11*n1p + s12*n2p + s13*n3p)
         dg(2,7) = -v*(s21*n1p + s22*n2p + s23*n3p)
         dg(3,7) = 0.

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
!        (Original Distance) x [(Perturbation pressure)*(Original normals)]
!        [r0] x (p.n0)            
         dgtq(1,3) = 0. ! dgtq(1,3) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq(2,3) = 0. ! dgtq(2,3) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq(3,3) = dgtq(3,3) + (r1*dg(2,1)-r2*dg(1,1))

!        Saving individual components            
         dgtq_p(1,5) = 0. ! dgtq_p(1,5) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq_p(2,5) = 0. ! dgtq_p(2,5) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq_p(3,5) = dgtq_p(3,5) + (r1*dg(2,1)-r2*dg(1,1))

!        [[Sn]*r0] x (P.n0) and  [r0] x (P.[Sn]*n0) Should both be zero
!        for pure rotation. Adding them anyway

!        Mean pressure torque.            
!        (Perturbation Distance) x [(Mean pressure)*(Original normals)]
!        [[Sn]*r0] x (P.n0)            
         dgtq(1,3) = 0. ! dgtq(1,3) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq(2,3) = 0. ! dgtq(2,3) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq(3,3) = dgtq(3,3) + (r1p*dg(2,6)-r2p*dg(1,6))

!        Saving individual components            
         dgtq_p(1,5) = 0. ! dgtq_p(1,5) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq_p(2,5) = 0. ! dgtq_p(2,5) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq_p(3,5) = dgtq_p(3,5) + (r1p*dg(2,1)-r2p*dg(1,1))

!        Mean pressure torque.            
!        (Original Distance) x [(Mean pressure)*(Perturbation normals)]
!        [r0] x (P.[Sn]*n0)            
         dgtq(1,3) = 0. ! dgtq(1,3) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq(2,3) = 0. ! dgtq(2,3) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq(3,3) = dgtq(3,3) + (r1*dg(2,5)-r2*dg(1,5))

!        Saving individual components            
         dgtq_p(1,5) = 0. ! dgtq_p(1,5) + (r2*dg(3,1)-r3*dg(2,1))
         dgtq_p(2,5) = 0. ! dgtq_p(2,5) + (r3*dg(1,1)-r1*dg(3,1))
         dgtq_p(3,5) = dgtq_p(3,5) + (r1*dg(2,5)-r2*dg(1,5))

!!       Viscous torque
!        (Perturbation Distance) x [(Mean Stress)*(Original normals)]
!        [Sn*r0] x (Tauij*n0)            
         dgtq(1,4) = 0. ! dgtq(1,4) + (r2p*dg(3,2)-r3p*dg(2,2))
         dgtq(2,4) = 0. ! dgtq(2,4) + (r3p*dg(1,2)-r1p*dg(3,2))
         dgtq(3,4) = dgtq(3,4) + (r1p*dg(2,2)-r2p*dg(1,2))

!        Saving individual components            
         dgtq_p(1,6) = 0. ! dgtq_p(1,6) + (r2p*dg(3,2)-r3p*dg(2,2))
         dgtq_p(2,6) = 0. ! dgtq_p(2,6) + (r3p*dg(1,2)-r1p*dg(3,2))
         dgtq_p(3,6) = dgtq_p(3,6) + (r1p*dg(2,2)-r2p*dg(1,2))

!        (Original Distance) x [(Mean Stress)*(Perturbation normals)]
!        (r0) x (Tauij*[Sn]*n0)            
         dgtq(1,4) = 0. ! dgtq(1,4) + (r2*dg(3,3)-r3*dg(2,3))
         dgtq(2,4) = 0. ! dgtq(2,4) + (r3*dg(1,3)-r1*dg(3,3))
         dgtq(3,4) = dgtq(3,4) + (r1*dg(2,3)-r2*dg(1,3))

!        Saving individual components            
         dgtq_p(1,7) = 0. ! dgtq_p(1,7) + (r2*dg(3,3)-r3*dg(2,3))
         dgtq_p(2,7) = 0. ! dgtq_p(2,7) + (r3*dg(1,3)-r1*dg(3,3))
         dgtq_p(3,7) = dgtq_p(3,7) + (r1*dg(2,3)-r2*dg(1,3))

!        (Original Distance) x [(Perturbation Stress)*(Original normals)]
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
      end subroutine clcdcm_pert2

!-----------------------------------------------------------------------
      subroutine fsi_calc_implicit_pert1

!     Calculate superposition of Solutions to get implicit update

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
      include 'FSI'

      real uy                 ! output vertical velocity

      integer icalld
      save icalld
      data icalld /0/

      real ext_k,bd_eta,bd_etav
      real ab0,ab1,ab2
      integer iobj
      integer ilag
      real const
      real solid_area

      real lhs,rhs,rhs2
!      real alpha
      real Fs,Fb,Fk_ext,Fg,Fdx,Ft

      integer irst,NBDMSH,NABMSH
      integer istep_old

      irst=0
!     ABBD is not set if there is no moving boundary condition
      IF (.not.IFMVBD) THEN
         NBDMSH = 1
         NABMSH = PARAM(28)
         IF (NABMSH.GT.ISTEP .and. irst.le.0) NABMSH = ISTEP
         IF (IFSURT)          NABMSH = NBD
         CALL RZERO   (ABMSH,10)
         CALL SETABBD (ABMSH,DTLAG,NABMSH,NBDMSH)
      ENDIF

!     Only doing for one object right now      
      iobj=0
      if (.not.fsi_ifrot) then
        Fs = fs_dragy(iobj) 
        Fg = fg_dragy(iobj)     ! Green's function
        Fb = fb_dragy(iobj)     ! base flow forces
      else  
!       Rotational coordinate system is opposite to the one used in Nek            
        Fs = -fs_torqz(iobj) 
        Fg = -fg_torqz(iobj)    ! Green's function
        Fb = -fb_torqz(iobj)
      endif                            

!     For non-linear case this is zero      
      Fdx = 0.

      if (ifpert.and.fsi_ifdispl) then
        if (fsi_ifrot) then
!         Rotational coordinate system is opposite to the one used in Nek 
          Fdx = -fdx_torqz(iobj) 
        else  
          Fdx = fdx_dragy(iobj) 
        endif
      endif  

!     Position
!     Copy lag array for eta 
      do ilag=lorder-1,2,-1
        etalag(ilag)=etalag(ilag-1)
      enddo
      etalag(1)=eta
      eta  = eta + (ABMSH(1)*etav + ABMSH(2)*etavlag(1)
     $       + ABMSH(3)*etavlag(2))*DT

      Fk_ext = -fsi_stiff*eta

!     Backward differentiation
      bd_etav = etav*bd(2)
      do ilag=2,NBD
        bd_etav = bd_etav + bd(ilag+1)*etavlag(ilag-1)
      enddo

!     Copy lag array for etav 
      do ilag=lorder-1,2,-1
        etavlag(ilag)=etavlag(ilag-1)
      enddo
      etavlag(1)=etav

!     Backward differentiation coefficients already assume the terms
!     have been moved to the rhs. So a -ve sign is already included.      
      rhs = Fs + Fb + Fdx + Fk_ext - fsi_inertia/DT*bd(1)*etav_s
     $         + fsi_inertia/DT*bd_etav + fsi_damp*etav_s
      lhs = fsi_inertia/DT*bd(1)*etav_g - Fg -fsi_damp*etav_g
      fsi_alpha = rhs/lhs               ! Eqn (49) in paper.
      
      etav = etav_s + fsi_alpha*etav_g

      etaa = (bd(1)*etav - bd_etav)/DT
     
!     This assumes jp=1      
      call opadd2cm(vxp,vyp,vzp,amp_vx,amp_vy,amp_vz,fsi_alpha)
      call add2s2(prp,amp_pr,fsi_alpha,nx2*ny2*nz2*nelt)

!     Get approximate velocity for the next time-step.
!     Calculation is independent since we correct this in the Stokes' step anyway.
!     Good projection allows the correction to be small.
!     According to paper:
      bd_etav = etav*bd(2)
      do ilag=2,NBD
        bd_etav=bd_etav+bd(ilag+1)*etavlag(ilag-1)
      enddo
      etav_s=-fsi_inertia/DT*bd_etav/(fsi_damp-fsi_inertia*bd(1)/DT)

!     Find position for next step.
!     Here I'm assuming that DT is constant.
!     Otherwise we need to update dtlag array as well.
!     And then revert it back for the normal nek solve.      
      if (ifpert) then
         istep_old = istep
         istep = istep+1
         NBDMSH = 1
         NABMSH = PARAM(28)
         IF (NABMSH.GT.ISTEP .and. irst.le.0) NABMSH = ISTEP
         IF (IFSURT)          NABMSH = NBD
         CALL RZERO   (ABMSH,10)
         CALL SETABBD (ABMSH,DTLAG,NABMSH,NBDMSH)

!        This is the position for the next step.         
         eta_s  = eta + (ABMSH(1)*etav + ABMSH(2)*etavlag(1)
     $       + ABMSH(3)*etavlag(2))*DT
         istep = istep_old
      endif 

!     All forces      
      fsi_Fs = Fs

      fsi_Fb = Fb 

      fsi_Fg = Fg

      fsi_Ff = fsi_Fs + fsi_Fb + fsi_alpha*fsi_Fg

      fsi_Fd = fsi_damp*etav

      fsi_Fdx = Fdx

      fsi_Fk = Fk_ext

      fsi_Ft = fsi_Ff+fsi_Fd+fsi_Fdx+fsi_Fk

      call fsi_output

      call fsi_rstsave 

      return
      end subroutine fsi_calc_implicit_pert1

!---------------------------------------------------------------------- 
      subroutine fsi_gengeom (igeom)
C
C     Generate geometry data.
c     Without updating bm1lag arrays            
C
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'GEOM'
      include 'WZ'
C
      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
C

      if (nio.eq.0.and.istep.le.1) write(6,*) 'generate geometry data'

      IF (IGEOM.EQ.1) THEN
         RETURN
      ELSEIF (IGEOM.EQ.2) THEN
         CALL LAGMASS
         IF (ISTEP.EQ.0) CALL GENCOOR (XM3,YM3,ZM3)
!        prabal. Update coordinates only for non-linear case
         IF (ISTEP.GE.1) CALL FSI_UPDCOOR
         CALL GEOM1 (XM3,YM3,ZM3)         ! generate Geometrical data mesh 1
         CALL GEOM2                       ! generate Geometrical data mesh 2
         CALL UPDMSYS (1)
         CALL VOLUME
         CALL SETINVM
         CALL SETDEF
         CALL SFASTAX
         IF (ISTEP.GE.1) CALL EINIT
      ELSEIF (IGEOM.EQ.3) THEN
c
c        Take direct stiffness avg of mesh
c
         ifieldo = ifield
         if (.not.ifneknekm) CALL GENCOOR (XM3,YM3,ZM3)
         if (ifheat) then
            ifield = 2
            CALL dssum(xm3,nx3,ny3,nz3)
            call col2 (xm3,tmult,ntot3)
            CALL dssum(ym3,nx3,ny3,nz3)
            call col2 (ym3,tmult,ntot3)
            if (if3d) then
               CALL dssum(xm3,nx3,ny3,nz3)
               call col2 (xm3,tmult,ntot3)
            endif
         else
            ifield = 1
            CALL dssum(xm3,nx3,ny3,nz3)
            call col2 (xm3,vmult,ntot3)
            CALL dssum(ym3,nx3,ny3,nz3)
            call col2 (ym3,vmult,ntot3)
            if (if3d) then
               CALL dssum(xm3,nx3,ny3,nz3)
               call col2 (xm3,vmult,ntot3)
            endif
         endif
         CALL GEOM1 (XM3,YM3,ZM3)
         CALL GEOM2
         CALL UPDMSYS (1)
         CALL VOLUME
         CALL SETINVM
         CALL SETDEF
         CALL SFASTAX
         ifield = ifieldo
      ENDIF

      if (nio.eq.0.and.istep.le.1) then
        write(6,*) 'done :: generate geometry data' 
        write(6,*) ' '
      endif

      return
      end subroutine fsi_gengeom
c-----------------------------------------------------------------------
      subroutine fsi_updcoor
C
C     Subroutine to update geometry for moving boundary problems
c     Without updating lag arrays for mesh velocity            
C
C-----------------------------------------------------------------------

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
C
      integer nel
      integer ifld

      ifld = ifield
      ifield = 0
      nel    = nelfld(ifield)
      ifield=ifld

C     Update collocation points coordinates

      if (.not.ifpert) then
        CALL UPDXYZ (NEL)
      endif
C
C     Shift lagged mesh velocity
C
      CALL LAGMSHV (NEL)
C
      return
      end subroutine fsi_updcoor
c-----------------------------------------------------------------------
      subroutine fsi_map12(pm2,pm1)

      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'

      real pm1(lx1*ly1*lz1,lelv)
      real pm2(lx2*ly2*lz2,lelv)
      integer e

      do e=1,nelv
        call map12 (pm2(1,e),pm1(1,e),e)
      enddo
   
      return
      end subroutine fsi_map12
c-----------------------------------------------------------------------


