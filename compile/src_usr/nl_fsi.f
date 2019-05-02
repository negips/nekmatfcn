!======================================================================
!     Routines for fluid structure interaction implicit solutions
!     using fixed-point method with dynamic relaxation      
!     Author: Prabal S. Negi
!======================================================================       
!     read parameters fluid-structure interaction 
      subroutine nlfsi_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'NLFSI'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /NLFSI/ ifnlfsi,psi_ini,psiv_ini,nlfsi_stiff,
     $                 nlfsi_damp,nlfsi_inertia,nlfsi_x0,nlfsi_y0,
     $                 nlfsi_rescale,nlfsi_tol,
     $                 nlfsi_rst_fli,nlfsi_iftermso

!     default values
      ifnlfsi             =  .FALSE.        ! if FSI using NL iterations  
      psi_ini             =  0.             ! initial position
      psiv_ini            =  0.             ! initial velocity
      nlfsi_stiff         =  1.             ! structural stiffness
      nlfsi_damp          = -1.             ! structural damping
      nlfsi_inertia       =  1.0            ! inertia
      nlfsi_x0            =  0.             ! elastic axis x0
      nlfsi_y0            =  0.             ! elastic axis y0
      nlfsi_rescale       =  1.0            ! scale the fluid forces
      nlfsi_tol           =  1.0E-08        ! Tolerance for structural solver
      nlfsi_rst_fli       =  0              ! restart file no.
      nlfsi_iftermso      = .FALSE.         ! if output terms for debugging

!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=NLFSI,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading NLFSI parameters.$')

!     broadcast data
      call bcast(ifnlfsi        , LSIZE)
      call bcast(psi_ini        ,WDSIZE)
      call bcast(psiv_ini       ,WDSIZE)
      call bcast(nlfsi_stiff    ,WDSIZE)
      call bcast(nlfsi_damp     ,WDSIZE)
      call bcast(nlfsi_inertia  ,WDSIZE)
      call bcast(nlfsi_x0       ,WDSIZE)
      call bcast(nlfsi_y0       ,WDSIZE)
      call bcast(nlfsi_rescale  ,WDSIZE)
      call bcast(nlfsi_tol      ,WDSIZE)
      call bcast(nlfsi_rst_fli  , ISIZE)
      call bcast(nlfsi_iftermso , LSIZE)

      nlfsi_ifinit = .false.

      return
      end
!----------------------------------------------------------------------
!     write parameters fluid-structure interaction
      subroutine nlfsi_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'NLFSI'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

      logical ifrot

!     namelists
      namelist /NLFSI/ ifnlfsi,psi_ini,psiv_ini,nlfsi_stiff,
     $                 nlfsi_damp,nlfsi_inertia,nlfsi_x0,nlfsi_y0,
     $                 nlfsi_rescale,nlfsi_tol,
     $                 nlfsi_rst_fli,nlfsi_iftermso,ifrot

      
!     This is set at compile time. Maybe I should change that?
!     Need an outpost in logfile if we need to check later      
      ifrot = nlfsi_ifrot

!     read the file
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=NLFSI,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing NLFSI parameters.$')

      return
      end
c-----------------------------------------------------------------------

      subroutine nlfsi_nek_advance

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'ADJOINT_DEF'
      include 'ADJOINT'       ! ifadj
      include 'NLFSI'


!     Only for FSI with non-linear iterations
      if (.not.ifnlfsi) then
        if (nio.eq.0) 
     $    write(6,*) 'Subroutine only for Non-Linear FSI'
        call exitt
      endif 

!     Not enabled for MHD
      if (ifmhd) then
        if (nio.eq.0) 
     $    write(6,*) 'MHD not enabled for Non-Linear FSI'
        call exitt
      endif  

!     Not enabled for NekNek
      if (ifneknekm) then
        if (nio.eq.0) 
     $    write(6,*) 'NekNek not enabled for Non-Linear FSI'
        call exitt
      endif  
!     Not enabled for ifheat
      if (ifheat) then
        if (nio.eq.0) 
     $    write(6,*) 'ifheat not enabled for Non-Linear FSI'
        call exitt
      endif  
!     Not enabled for CMT
      if (ifcmt) then
        if (nio.eq.0) 
     $    write(6,*) 'CMT not enabled for Non-Linear FSI'
        call exitt
      endif
!     Not enabled for PN/PN
      if (ifsplit) then
        if (nio.eq.0) 
     $    write(6,*) 'PnPn not enabled for Non-Linear FSI'
        call exitt
      endif  

!     Only for transient simulations
      if (ifsplit) then
        if (nio.eq.0) 
     $    write(6,*) 'Only for transient Non-Linear FSI(iftrans=.true.)'
        call exitt
      endif 

!     Rotational adjoint forces need to be derived (and implemented).      
      if (ifpert.and.ifadj.and.nlfsi_ifrot) then
        if (nio.eq.0) then
          write(6,*) 'Adjoint not implemented for rotational FSI'
          write(6,*) 'Adjoint forces not implemented.'
        endif  
        call exitt
      endif 

      if (ifpert) then
        call nlfsi_advancep
      else
        call nlfsi_advance
      endif  

      return
      end subroutine nlfsi_nek_advance

c-----------------------------------------------------------------------

      subroutine nlfsi_advance

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
      include 'NLFSI'

      include 'FSI_DEBUG'

      integer igeom
      common /cgeom/ igeom    ! Apparently this common block is in NekNek

      integer icalld2
      save icalld2
      data icalld2 /0/

      integer lt
      parameter (lt=lx1*ly1*lz1*lelv)

      integer lt2
      parameter (lt2=lx2*ly2*lz2*lelv)

      character str*14
      integer ilag

      real wk1(lt),wk2(lt),wk3(lt)

      call nekgsync
      if (iftran) call settime

!     PN-2/PN-2 formulation
      call setprop
        
      call setsolv
      call comment
            
      nlfsi_ifconv = .false.
      nlfsi_iter = 0
      str='NL FSI Solve: '

!     Save current time-step fields
      call opcopy(vx0,vy0,vz0,vx,vy,vz)
      call copy(pr0,pr,lt2)
      call copy(bm0,bm1,lt)
      call opcopy(nlfsi_wx,nlfsi_wy,nlfsi_wz,wx,wy,wz)

      do while (.not.nlfsi_ifconv)
    
        nlfsi_iter = nlfsi_iter+1

        do igeom=1,ngeom

          if (ifgeom) then
            call nlfsi_gengeom (igeom)
            call geneig  (igeom)
          endif

          ifield = 1
          imesh = 1
          call unorm
          call settolv

          if (ifflow) call nlfsi_fluid_plan3(igeom)

          if (igeom.ge.2) call vol_flow        ! check for fixed flow rate

        enddo
            
        call nlfsi_main

        if (nio.eq.0) write(6,11) str,istep,time,
     $             nlfsi_iter,nlfsi_timei,nlfsi_tol,
     $             nlfsi_rnorm(nlfsi_iter),psi_res(nlfsi_iter),
     $             psiv_res(nlfsi_iter),nlfsi_omegai(nlfsi_iter) 
   11   format(A14,I7,1x,1pE14.7E2,1x,I3,1x,1pE11.4E2,1x,1pE9.2E2,1x,
     $         4(1pE14.7E2,1x))

!       For debugging            
!        if (nlfsi_iter.eq.2) nlfsi_ifconv=.true.

        if (.not.nlfsi_ifconv) then       ! not converged
!         revert to old velocities                
          call opcopy(vx,vy,vz,vx0,vy0,vz0)
          call copy(pr,pr0,lt2)

!         revert to old mesh
          call opcmult(wx,wy,wz,-1.) 
          do ilag=1,nbd-1
            call opcmult(wxlag(1,1,1,1,ilag),wylag(1,1,1,1,ilag),
     $                   wzlag(1,1,1,1,ilag),-1.)
          enddo  
          do igeom=1,ngeom
            if (ifgeom) then
              call nlfsi_gengeom (igeom)
              call geneig  (igeom)
            endif
          enddo
!         Revert sign of mesh velocities            
          call opcmult(wx,wy,wz,-1.)  
          do ilag=1,nbd-1
            call opcmult(wxlag(1,1,1,1,ilag),wylag(1,1,1,1,ilag),
     $                   wzlag(1,1,1,1,ilag),-1.)
          enddo 

        endif 

        if (ifusermv) call nlfsi_meshv

      enddo       ! end of non-linear iterations

      call nlfsi_rstsave
      call nlfsi_output

      call setup_convect (igeom) ! Save convective velocity

      call nlfsi_lagvel
      call nlfsi_lagpr
      call nlfsi_lagbf
      call nlfsi_lagmass
      call nlfsi_lagmshv(nelv)

      return
      end subroutine nlfsi_advance
c-----------------------------------------------------------------------

      subroutine nlfsi_init

      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'NLFSI'

      integer lt2
      parameter (lt2=lx2*ly2*lz2*lelv)
      integer lt
      parameter (lt=lx1*ly1*lz1*lelv)

      real scale
      real x0(3)
      logical ifdout,iftout
      integer ierr

      real pm1(lx1,ly1,lz1,lelv)

      integer ifld,iface,nfaces,iel
      character cb*3

!     Initialize timer            
      nlfsi_timea = 0.
      nlfsi_timei = 0.

      nlfsi_iter = 0

!      Done in userchk
!      call set_obj
!      Put a check for number of objects. Prabal

!     This is ridiculous    
      call rzero(nlfs_dragx (0),maxobj+1)
      call rzero(nlfs_dragpx(0),maxobj+1)
      call rzero(nlfs_dragvx(0),maxobj+1)
      call rzero(nlfs_dragy (0),maxobj+1)
      call rzero(nlfs_dragpy(0),maxobj+1)
      call rzero(nlfs_dragvy(0),maxobj+1)
      call rzero(nlfs_dragz (0),maxobj+1)
      call rzero(nlfs_dragpz(0),maxobj+1)
      call rzero(nlfs_dragvz(0),maxobj+1)
      call rzero(nlfs_torqx (0),maxobj+1)
      call rzero(nlfs_torqpx(0),maxobj+1)
      call rzero(nlfs_torqvx(0),maxobj+1)
      call rzero(nlfs_torqy (0),maxobj+1)
      call rzero(nlfs_torqpy(0),maxobj+1)
      call rzero(nlfs_torqvy(0),maxobj+1)
      call rzero(nlfs_torqz (0),maxobj+1)
      call rzero(nlfs_torqpz(0),maxobj+1)
      call rzero(nlfs_torqvz(0),maxobj+1)

      call rzero(psilag,lorder-1)
      call rzero(psivlag,lorder-1)
      call rzero(psialag,lorder-1)
      call rzero(psivmlag,lorder-1)

      if (ifpert) then
        call gradm1(nlfsi_grvx0(1,1,1,1,1),nlfsi_grvx0(1,1,1,1,2),
     $              nlfsi_grvx0(1,1,1,1,3),vx)
        call gradm1(nlfsi_grvy0(1,1,1,1,1),nlfsi_grvy0(1,1,1,1,2),
     $              nlfsi_grvy0(1,1,1,1,3),vy)
        if (if3d) then
          call gradm1(nlfsi_grvz0(1,1,1,1,1),nlfsi_grvz0(1,1,1,1,2),
     $                nlfsi_grvz0(1,1,1,1,3),vz)
        else
          call rzero(nlfsi_grvz0(1,1,1,1,1),lt)        
          call rzero(nlfsi_grvz0(1,1,1,1,2),lt)        
          call rzero(nlfsi_grvz0(1,1,1,1,3),lt)        
          call rzero(nlfsi_grvx0(1,1,1,1,3),lt)        
          call rzero(nlfsi_grvy0(1,1,1,1,3),lt)        
        endif
          
      endif

      if (nid.eq.0) then

!       Open file to output psi,psiv,alpha               
        call IO_freeid(nlfsi_iunit1,ierr)
        open(nlfsi_iunit1,file=nlfsi_fname1,status='unknown',
     $       action='write',iostat=ierr)
        
        write(nlfsi_iunit1,'(A10,1x,5(A19,1x))') 
     $          'istep',"time",'psi','psiv','psia','psiv_m'
        flush(nlfsi_iunit1)


!       Open file to output terms of the structural solve
!       Mostly for debugging.
        if (nlfsi_iftermso) then          
          call IO_freeid(nlfsi_iunit2,ierr)
          open(nlfsi_iunit2,file=nlfsi_fname2,status='unknown',
     $         action='write',iostat=ierr)
          write(nlfsi_iunit2,'(A10,1x,8(A18,1x),A5)')
     $        'istep','Time','Fs','Fd','Fk','Ft','Res. Norm',
     $        'psiv Res.','psi Res.','Niter'
          flush(nlfsi_iunit2)
        endif  
      endif


!     Build a mask to remove only 'v  '
!     and preserving the 'mv ' points.
      if (ifpert) then
        call rone(nlfsi_mask,lx1*ly1*lz1*lelv)
        ifld = 1
        nfaces=2*ndim
        do iel=1,nelv
          do iface = 1,nfaces
            cb = cbc(iface,iel,ifld)
            if (cb.eq.'v  ') then
!             Put zeros on this face
              call facev(nlfsi_mask,iel,iface,0.,nx1,ny1,nz1)
            endif
          enddo ! iel
        enddo   ! iface

!       Weight for structural variables for inner-product in Arnoldi 
        nlfsi_psi_wt  = 1.0
        nlfsi_psiv_wt = nlfsi_inertia

      endif     ! ifpert  


!     initialization done      
      nlfsi_ifinit = .true.

      return 
      end subroutine nlfsi_init
!----------------------------------------------------------------------       

      subroutine nlfsi_fluid_plan3(igeom)
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
      IF (IGEOM.EQ.1) THEN
C
C        Old geometry
C
         CALL NLFSI_MAKEF
C
      ELSE
C
C        New geometry, new b.c.
C

         INTYPE = -1
         CALL SETHLM  (H1,H2,INTYPE)
         CALL NLFSI_CRESVIF (RESV1,RESV2,RESV3,H1,H2)

         mstep = abs(param(94))
         if (param(94).ne.0. .and. istep.ge.mstep) then
          call ophinvpr(dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxh)
c         CALL OPHINV  (DV1,DV2,DV3,RESV1,RESV2,RESV3,H1,H2,TOLHV,NMXH)
         else
           CALL OPHINV  (DV1,DV2,DV3,RESV1,RESV2,RESV3,H1,H2,TOLHV,NMXH)
         endif

         CALL OPADD2  (VX,VY,VZ,DV1,DV2,DV3)

         call incomprn_uzawa(vx,vy,vz,pr)

      ENDIF

      return
      end subroutine nlfsi_fluid_plan3    

c-----------------------------------------------------------------------

      subroutine nlfsi_makef
C
C     Compute and add: (1) user specified forcing function (FX,FY,FZ)
C                      (2) driving force due to natural convection
C                      (3) convection term
C
C     !! NOTE: Do not change the arrays BFX, BFY, BFZ until the
C              current time step is completed.
C
C----------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'

                                                call makeuf   ! user forcing 
      if (ifnatc)                               call natconv  ! natural convection
!      if (ifexplvis.and.ifsplit)                call explstrs
      if (ifnav.and.(.not.ifchar))              call advab    ! convection term
      if (ifmvbd)                               call admeshv  ! trilinear term due to moving mesh
      if (iftran)                               call nlfsi_makeabf  ! extrapolate RHS
                                                              ! without updating lag arrays for 
                                                              ! forcing (abx1..)
      if ((iftran.and..not.ifchar).or.
     $    (iftran.and..not.ifnav.and.ifchar))   call makebdf  ! backward time difference terms
!      if (ifnav.and.ifchar.and.(.not.ifmvbd))   call advchar  ! ifcharacteristics advection terms
!      if (ifmodel)                              call twallsh  ! RANS model terms

      return
      end subroutine nlfsi_makef
c-----------------------------------------------------------------------

      subroutine nlfsi_extrapp(p,plag)
C
C     Pressure extrapolation
C
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'

      real  p    (lx2,ly2,lz2,1)
     $     ,plag (lx2,ly2,lz2,1)

      common /cgeom/ igeom

      ntot2 = nx2*ny2*nz2*nelv

      if (nbd.eq.2.and.nbdinp.gt.2.and.igeom.le.2) then
         call copy(plag,p,ntot2)
      elseif (nbd.eq.3.and.igeom.le.2) then

         const = dtlag(1)/dtlag(2)

         do i=1,ntot2
            pnm1          = p   (i,1,1,1)
            pnm2          = plag(i,1,1,1)
            p   (i,1,1,1) = pnm1 + const*(pnm1-pnm2)
!           We save the lag pressure later 
!            plag(i,1,1,1) = pnm1
         enddo

      elseif (nbd.gt.3) then
         WRITE (6,*) 'Pressure extrapolation cannot be completed'
         WRITE (6,*) 'Try a lower-order temporal scheme'
         call exitt
      endif
      return
      end subroutine nlfsi_extrapp
c-----------------------------------------------------------------------

      subroutine nlfsi_cresvif (resv1,resv2,resv3,h1,h2)

!     The difference from the original is that the lag velocity
!     arrays are not updated yet            

      include 'SIZE'
      include 'TOTAL'
      REAL           RESV1 (LX1,LY1,LZ1,1)
      REAL           RESV2 (LX1,LY1,LZ1,1)
      REAL           RESV3 (LX1,LY1,LZ1,1)
      REAL           H1    (LX1,LY1,LZ1,1)
      REAL           H2    (LX1,LY1,LZ1,1)
      COMMON /SCRUZ/ W1    (LX1,LY1,LZ1,LELV)
     $ ,             W2    (LX1,LY1,LZ1,LELV)
     $ ,             W3    (LX1,LY1,LZ1,LELV)

      common /cgeom/ igeom

      NTOT1 = NX1*NY1*NZ1*NELV
      NTOT2 = NX2*NY2*NZ2*NELV
!      if (igeom.eq.2) CALL LAGVEL 
      CALL BCDIRVC (VX,VY,VZ,v1mask,v2mask,v3mask)
      IF (IFSTRS)  CALL BCNEUTR

!     prlag not updated here      
      call nlfsi_extrapp (pr,prlag)

      call opgradt (resv1,resv2,resv3,pr)
      CALL OPADD2  (RESV1,RESV2,RESV3,BFX,BFY,BFZ)
      CALL OPHX    (W1,W2,W3,VX,VY,VZ,H1,H2)
      CALL OPSUB2  (RESV1,RESV2,RESV3,W1,W2,W3)
C
      RETURN
      END
c-----------------------------------------------------------------------

      subroutine nlfsi_makeabf
C-----------------------------------------------------------------------
C
C     Sum up contributions to kth order extrapolation scheme.
c     NOTE: rho^{n+1} should multiply all the Sum_q{bpsi_q} term 
c           if rho is not constant!

C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'NLFSI'

C
      COMMON /SCRUZ/ TA1 (LX1,LY1,LZ1,LELV)
     $ ,             TA2 (LX1,LY1,LZ1,LELV)
     $ ,             TA3 (LX1,LY1,LZ1,LELV)
C
      NTOT1 = NX1*NY1*NZ1*NELV
C
      AB0 = AB(1)
      AB1 = AB(2)
      AB2 = AB(3)
      CALL ADD3S2 (TA1,ABX1,ABX2,AB1,AB2,NTOT1)
      CALL ADD3S2 (TA2,ABY1,ABY2,AB1,AB2,NTOT1)
!      CALL COPY   (ABX2,ABX1,NTOT1)
!      CALL COPY   (ABY2,ABY1,NTOT1)
!      CALL COPY   (ABX1,BFX,NTOT1)
!      CALL COPY   (ABY1,BFY,NTOT1)
!     prabal. New array ABX0,ABY0
      call copy(ABX0,BFX,NTOT1)
      call copy(ABY0,BFY,NTOT1)

      CALL ADD2S1 (BFX,TA1,AB0,NTOT1)
      CALL ADD2S1 (BFY,TA2,AB0,NTOT1)
      CALL COL2   (BFX,VTRANS,NTOT1)          ! multiply by density
      CALL COL2   (BFY,VTRANS,NTOT1)
      IF (NDIM.EQ.3) THEN
         CALL ADD3S2 (TA3,ABZ1,ABZ2,AB1,AB2,NTOT1)
!         CALL COPY   (ABZ2,ABZ1,NTOT1)
!         CALL COPY   (ABZ1,BFZ,NTOT1)
!        prabal. New array   
         call copy(ABZ0,BFZ,NTOT1)

         CALL ADD2S1 (BFZ,TA3,AB0,NTOT1)
         CALL COL2   (BFZ,VTRANS,NTOT1)
      ENDIF
C
      return
      end subroutine nlfsi_makeabf
C
c-----------------------------------------------------------------------

      subroutine nlfsi_lagvel
C
C     Keep old velocity field(s) 
C
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
C
      integer ilag,ntot1

      ntot1 = nx1*ny1*nz1*nelv
c
!      do 100 ilag=nbdinp-1,2,-1
      do 100 ilag=3-1,2,-1
         call copy (vxlag (1,1,1,1,ilag),vxlag (1,1,1,1,ilag-1),ntot1)
         call copy (vylag (1,1,1,1,ilag),vylag (1,1,1,1,ilag-1),ntot1)
         if (ndim.eq.3)
     $   call copy (vzlag (1,1,1,1,ilag),vzlag (1,1,1,1,ilag-1),ntot1)
 100  continue
c
      call opcopy (vxlag,vylag,vzlag,vx0,vy0,vz0)
c
      return
      end subroutine nlfsi_lagvel
c-----------------------------------------------------------------------

      subroutine nlfsi_lagpr

      implicit none

      INCLUDE 'SIZE_DEF'
      INCLUDE 'SIZE'
      INCLUDE 'SOLN_DEF'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'
      INCLUDE 'NLFSI'

      integer ntot2

      if (nbdinp.eq.3) then
         ntot2 = nx2*ny2*nz2*nelv
         call copy (prlag,pr0,ntot2)
      endif

      return
      end subroutine nlfsi_lagpr
c-----------------------------------------------------------------------

      subroutine nlfsi_lagmass
C
C     Lag the mass matrix (matrices)
C
      implicit none

      INCLUDE 'SIZE_DEF'      
      INCLUDE 'SIZE'
      INCLUDE 'MASS_DEF'
      INCLUDE 'MASS'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'
      INCLUDE 'NLFSI'

      integer ntot1,ilag
C
      ntot1 = nx1*ny1*nz1*nelt
      do 100 ilag=nbdinp-1,2,-1
         call copy (bm1lag(1,1,1,1,ilag),bm1lag(1,1,1,1,ilag-1),ntot1)
 100  continue
      call copy (bm1lag(1,1,1,1,1),bm0,ntot1)
c
      return
      end subroutine nlfsi_lagmass
c-----------------------------------------------------------------------

      subroutine nlfsi_lagbf

      implicit none

      INCLUDE 'SIZE_DEF'
      INCLUDE 'SIZE'
      INCLUDE 'SOLN_DEF'
      INCLUDE 'SOLN'
      INCLUDE 'NLFSI'

      INTEGER NTOT1

      NTOT1 = NX1*NY1*NZ1*NELV


      CALL COPY   (ABX2,ABX1,NTOT1)
      CALL COPY   (ABX1,ABX0,NTOT1)

      CALL COPY   (ABY2,ABY1,NTOT1)
      CALL COPY   (ABY1,ABY0,NTOT1)

      IF (NDIM.EQ.3) THEN
        CALL COPY   (ABZ2,ABZ1,NTOT1)
        CALL COPY   (ABZ1,ABZ0,NTOT1)
      ENDIF

      return
      end subroutine nlfsi_lagbf
c-----------------------------------------------------------------------

      subroutine nlfsi_lagmshv (nel)
C
C     Keep old mesh velocity
C
      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF' 
      include 'INPUT'
      include 'MVGEOM_DEF'
      include 'MVGEOM'
      include 'NLFSI'

      integer ilag,ntot1,nel

      ntot1 = nx1*ny1*nz1*nel

      do 100 ilag=nbdinp-1,2,-1
         call copy (wxlag(1,1,1,1,ilag),wxlag(1,1,1,1,ilag-1),ntot1)
         call copy (wylag(1,1,1,1,ilag),wylag(1,1,1,1,ilag-1),ntot1)
         if (ndim.eq.3)
     $   call copy (wzlag(1,1,1,1,ilag),wzlag(1,1,1,1,ilag-1),ntot1)
 100  continue

      call opcopy(wxlag,wylag,wzlag,nlfsi_wx,nlfsi_wy,nlfsi_wz)   

      return
      end subroutine nlfsi_lagmshv
c-----------------------------------------------------------------------

      subroutine nlfsi_gengeom (igeom)
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
!         CALL LAGMASS
         IF (ISTEP.EQ.0) CALL GENCOOR (XM3,YM3,ZM3)
!        prabal. Update coordinates without changing lag arrays
!        for mesh velocity         
         IF (ISTEP.GE.1) CALL NLFSI_UPDCOOR
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
      end subroutine nlfsi_gengeom
c-----------------------------------------------------------------------

      subroutine nlfsi_updcoor
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
!      CALL LAGMSHV (NEL)
C
      return
      end subroutine nlfsi_updcoor
c-----------------------------------------------------------------------

      subroutine nlfsi_output

      implicit none

      INCLUDE 'SIZE_DEF'
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'
      INCLUDE 'NLFSI'

      integer it

      it = nlfsi_iter
!     if things are not initialized      
      if (it.eq.0) it=1

      if (nid.eq.0) then
!       Output psi,psiv,alpha            
        write(nlfsi_iunit1,'(I10,1x,5(1pE19.11E3,1x))') istep, 
     $       time,psi,psiv,psia,psiv_m
        flush(nlfsi_iunit1) 

!       If debugging is on
        if (nlfsi_iftermso) then
          write(nlfsi_iunit2,'(I10,1x,8(1pE18.8E3,1x),I5)') istep, 
     $       time,nlfsi_Ff,nlfsi_Fd,nlfsi_Fk,nlfsi_Ft,
     $       nlfsi_rnorm(it),psiv_res(it),psi_res(it),nlfsi_iter
          flush(nlfsi_iunit2)
        endif           ! iftermso 
      endif             ! nid.eq.0

      return
      end subroutine nlfsi_output            
!----------------------------------------------------------------------       



