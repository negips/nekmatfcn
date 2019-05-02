!====================================================================== 
      subroutine nlfsi_main1

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


      if (icalld.eq.0) then
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
        if (ifusermv.and.ifpert) call nlfsi_meshv_pert
        if (ifusermv.and..not.ifpert) call nlfsi_meshv
        return
      endif  

!     Start timing for this step   
      steptime = dnekclock()

      if (ifchkptrst.and.istep.lt.chkptnrsf) then
!       Read saved values for restart
!       If the simulation has been restarted
        call nlfsi_rstread

        call nlfsi_output

        if (ifusermv) call nlfsi_meshv

        nlfsi_timei=dnekclock()-steptime
        nlfsi_timea=nlfsi_timea+nlfsi_timei
        return

      elseif (istep.eq.0) then
!       If there are no restarts              
!       At first time step the projected velocity is
!       just the initialized velocity
        psiv_s = psiv_ini
        psiv_m = psiv_ini
        psiv   = psiv_ini
        psia   = 0.
        psia_s = 0.    
        psi    = psi_ini
        psi_s  = psi + psiv_m*DT

        call nlfsi_output

        if (ifusermv.and.ifpert) call nlfsi_meshv_pert
        if (ifusermv.and..not.ifpert) call nlfsi_meshv

        nlfsi_timei=dnekclock()-steptime
        nlfsi_timea=nlfsi_timea+nlfsi_timei
        return
      endif

      x0(1) = nlfsi_x0
      x0(2) = nlfsi_y0
      if (ndim.eq.3) x0(3) = 1000.         ! set outside the domain

      ifdout=.true.
      iftout=.true.

      scale=nlfsi_rescale
      
      if (ifpert) then
        call nlfsi_torque_pert(scale,x0,ifdout,iftout,
     $      vx,vy,vz,pr,vxp(1,jp),vyp(1,jp),vzp(1,jp),prp(1,jp))
      else            
        call torque_calc_axis(scale,x0,ifdout,iftout,vx,vy,vz,pr)
      endif  

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

      if (ifpert.and.nlfsi_ifdispl) then

!       Zero contribution for Mean field   
        call opzero(wk_vx0,wk_vy0,wk_vz0)
        call rzero(wk_pr0,lvp2) 

        if (nlfsi_ifrot) then

          call rzero(rot_sn,9)            ! sine terms
          call rzero(rot_cs,9)            ! cosine terms
!         Matrix of sine terms
!         For a clockwise rotation for a positive psi
          rot_sn(1,2) =  psi_s
          rot_sn(2,1) = -psi_s

!         Matrix of cosine terms
!         In a linear case this is just Identity      
          rot_cs(1,1) = 1
          rot_cs(2,2) = 1
          rot_cs(3,3) = 1

!         Distance from Rotational axis
          call opcopy(rxyz(1,1),rxyz(1,2),rxyz(1,3),xm1,ym1,zm1) 
          call cadd(rxyz(1,1),-nlfsi_x0,lvp)
          call cadd(rxyz(1,2),-nlfsi_y0,lvp)
          call rzero(rxyz(1,3),lvp)       ! rotation is about 'z'

          do i=1,lvp 
            r0(1) = rxyz(i,1) 
            r0(2) = rxyz(i,2)
            r0(3) = rxyz(i,3)

!           perturbation in position
            call mxm(rot_sn,3,r0,3,dxyz,1)
            psi_sx = dxyz(1)
            psi_sy = dxyz(2)
            psi_sz = dxyz(3) 

!           Calculate Grad(U0).psi            
            wk_vxp(i) =  nlfsi_grvx0(i,1,1,1,1)*psi_sx + 
     $           nlfsi_grvx0(i,1,1,1,2)*psi_sy +
     $           nlfsi_grvx0(i,1,1,1,3)*psi_sz
            wk_vyp(i) =  nlfsi_grvy0(i,1,1,1,1)*psi_sx + 
     $           nlfsi_grvy0(i,1,1,1,2)*psi_sy +
     $           nlfsi_grvy0(i,1,1,1,3)*psi_sz
            if (if3d) then
              wk_vzp(i) = nlfsi_grvz0(i,1,1,1,1)*psi_sx + 
     $               nlfsi_grvz0(i,1,1,1,2)*psi_sy +
     $               nlfsi_grvz0(i,1,1,1,3)*psi_sz
            else
              wk_vzp(i)=0.
            endif

!           Calculate Grad(P0).psi            
            tmp_pr(i) =  nlfsi_grpr0(i,1,1,1,1)*psi_sx + 
     $           nlfsi_grpr0(i,1,1,1,2)*psi_sy +
     $           nlfsi_grpr0(i,1,1,1,3)*psi_sz
          enddo 
!         Map this to pressure grid            
          call nlfsi_map12(wk_prp,tmp_pr) 
        else    

!         Calculate Grad(U0).psi            
          call opcopy(wk_vxp,wk_vyp,wk_vzp,nlfsi_grvx0(1,1,1,1,2),
     $                nlfsi_grvy0(1,1,1,1,2),nlfsi_grvz0(1,1,1,1,2))
          call opcmult(wk_vxp,wk_vyp,wk_vzp,psi_s) 

!         Calculate Grad(P0).psi            
          call copy(tmp_pr,nlfsi_grpr0(1,1,1,1,2),lvp)
          call cmult(tmp_pr,psi_s,lvp)
!         Map this to pressure grid            
          call nlfsi_map12(wk_prp,tmp_pr)
        endif  

!       These quantities are treated as perturbation fields   
        call nlfsi_torque_pert(scale,x0,ifdout,iftout,
     $     wk_vx0,wk_vy0,wk_vz0,wk_pr0,wk_vxp,wk_vyp,wk_vzp,wk_prp)
        
        call copy(nlfdx_dragx (0),dragx (0),maxobj+1)
        call copy(nlfdx_dragpx(0),dragpx(0),maxobj+1)
        call copy(nlfdx_dragvx(0),dragpx(0),maxobj+1)
        call copy(nlfdx_dragy (0),dragy (0),maxobj+1)
        call copy(nlfdx_dragpy(0),dragpy(0),maxobj+1)
        call copy(nlfdx_dragvy(0),dragvy(0),maxobj+1)
        call copy(nlfdx_dragz (0),dragz (0),maxobj+1)
        call copy(nlfdx_dragpz(0),dragpz(0),maxobj+1)
        call copy(nlfdx_dragvz(0),dragvz(0),maxobj+1)
        call copy(nlfdx_torqx (0),torqx (0),maxobj+1)
        call copy(nlfdx_torqpx(0),torqpx(0),maxobj+1)
        call copy(nlfdx_torqvx(0),torqpx(0),maxobj+1)
        call copy(nlfdx_torqy (0),torqy (0),maxobj+1)
        call copy(nlfdx_torqpy(0),torqpy(0),maxobj+1)
        call copy(nlfdx_torqvy(0),torqvy(0),maxobj+1)
        call copy(nlfdx_torqz (0),torqz (0),maxobj+1)
        call copy(nlfdx_torqpz(0),torqpz(0),maxobj+1)
        call copy(nlfdx_torqvz(0),torqvz(0),maxobj+1)
      endif   

      call nlfsi_struc_solve_acc

      nlfsi_timei=dnekclock()-steptime
      nlfsi_timea=nlfsi_timea+nlfsi_timei

      return
      end subroutine nlfsi_main1 
!---------------------------------------------------------------------- 

      subroutine nlfsi_meshv_pert1

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
          if (nlfsi_ifale) then  
            ucx =  0.
            ucy =  psiv_m
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
          psi_sx = 0.
          psi_sy = psi_s
          psi_sz = 0.  

          if (nlfsi_bcdispl) then
            ucx =  0. - (nlfsi_grvx0(i,1,1,1,1)*psi_sx + 
     $             nlfsi_grvx0(i,1,1,1,2)*psi_sy +
     $             nlfsi_grvx0(i,1,1,1,3)*psi_sz) 
            ucy =  psiv_s - (nlfsi_grvy0(i,1,1,1,1)*psi_sx + 
     $             nlfsi_grvy0(i,1,1,1,2)*psi_sy +
     $             nlfsi_grvy0(i,1,1,1,3)*psi_sz)
            if (if3d) then
              ucz =  0. - (nlfsi_grvz0(i,1,1,1,1)*psi_sx + 
     $               nlfsi_grvz0(i,1,1,1,2)*psi_sy +
     $               nlfsi_grvz0(i,1,1,1,3)*psi_sz)
            else
              ucz = 0.
            endif 
          else
            ucx = 0.
            ucy = psiv_s
            ucz = 0.
          endif          

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
          if (nlfsi_ifale) then  
            ucx =  psiv_m*dy
            ucy =  -psiv_m*dx
            ucz =  0.
          else
            ucx =  0.
            ucy =  0.
            ucz =  0.
          endif  
  
          wx(i,1,1,1) = basev(i)*ucx
          wy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) wz(i,1,1,1) = basev(i)*ucz

          if (nlfsi_bcdispl) then

!           Calculate distance from axis   
            r0(1) = dx 
            r0(2) = dy
            r0(3) = dz

!           perturbation in position
            call mxm(rot_sn,3,r0,3,dxyz,1)
            psi_sx = dxyz(1)
            psi_sy = dxyz(2)
            psi_sz = dxyz(3) 

            ucx =  psiv_s*dy
            ucy =  -psiv_s*dx
            ucz =  0.

            ucx =  ucx - (nlfsi_grvx0(i,1,1,1,1)*psi_sx + 
     $             nlfsi_grvx0(i,1,1,1,2)*psi_sy +
     $             nlfsi_grvx0(i,1,1,1,3)*psi_sz) 
            ucy =  ucy - (nlfsi_grvy0(i,1,1,1,1)*psi_sx + 
     $             nlfsi_grvy0(i,1,1,1,2)*psi_sy +
     $             nlfsi_grvy0(i,1,1,1,3)*psi_sz)
            if (if3d) then
              ucz =  ucz - (nlfsi_grvz0(i,1,1,1,1)*psi_sx + 
     $               nlfsi_grvz0(i,1,1,1,2)*psi_sy +
     $               nlfsi_grvz0(i,1,1,1,3)*psi_sz)
            else
              ucz = 0.
            endif 
          else          ! not using displacement effects in BCs
            ucx =  psiv_s*dy
            ucy =  -psiv_s*dx
            ucz =  0.
          endif          
         
          umeshx(i,1,1,1) = basev(i)*ucx
          umeshy(i,1,1,1) = basev(i)*ucy
          if (ndim.eq.3) umeshz(i,1,1,1) = basev(i)*ucz
         
        endif
       
      enddo

      return
      end subroutine nlfsi_meshv_pert1
!---------------------------------------------------------------------- 

