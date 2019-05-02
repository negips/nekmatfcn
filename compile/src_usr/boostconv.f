!=======================================================================
!     Clio

!     Tools for Immersed boundary fin
!     Parameters used by this set of subroutines:
!     LEVELSET:
!     DIFTIME, DIFSTEP, DIFK - diffusion of the mask boundary
!                              (smooth interface)
!     PHIMAX, PHIMIN - initial conditions for level set, max and min
!     AMPMSK1, AMPMSK2 - immersed boundary coefficients alpha(linear)
!                        and beta(integral)
!     IC_XMAX, IC_YMAX, IC_YMIN - initial geometry of the fin
!     IBM_PROP - coefficients for the governing equations
!=======================================================================
!***********************************************************************
!> @brief Read runtime parameters for arnoldi_arpack
!! @ingroup arnoldi_arpack
!! @note This interface is defined in @ref tstpr_param_get
!! @todo Check iostat value for missing namelist in the file
!! @see @ref readers_writers_page
!! @author Nicolas Offermans
      subroutine boostconv_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF'
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'BOOSTCONV'       ! BC_FREQ, BC_TOL, IFBOOSTCONV

!     local variables
      integer ierr, fid

!     namelists
      namelist /BOOSTCONV/ BC_FREQ, BC_TOL, IFBOOSTCONV,IFBOOSTVERBOSE
     $                    ,BOOSTVER       

!-----------------------------------------------------------------------
!     default values
      BC_FREQ       = 50
      BC_TOL        = 1.000E-10
      IFBOOSTCONV   = .true.
      IFBOOSTVERBOSE= .true.
      BOOSTVER      = 1             ! 1: PSN, 0:NOF

!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=BOOSTCONV,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading BOOSTCONV parameters.$')

!     broadcast data
      call bcast(BC_FREQ,ISIZE)
      call bcast(BC_TOL,WDSIZE)
      call bcast(IFBOOSTCONV,LSIZE)
      call bcast(IFBOOSTVERBOSE,LSIZE)
      call bcast(BOOSTVER,ISIZE)
     
      return
      end
!***********************************************************************
!     write parameters BOOSTCONV
      subroutine boostconv_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'BOOSTCONV'       ! BC_FREQ, BC_TOL, IFBOOSTCONV 

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /BOOSTCONV/ BC_FREQ, BC_TOL, IFBOOSTCONV, IFBOOSTVERBOSE
     $                    ,BOOSTVER  
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=BOOSTCONV,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing BOOSTCONV parameters.$')

      return
      end
!***********************************************************************
!     Boost convergence of an iterative algorithm:
!     estimates the next iteration 
!     as a linear combination of the previous iterations.
!     It is similar to a GMRES, 
!     but continuously updates the subspace with the base vectors.
!     
!     Parameters:
!     
!     ws: pointer to an internal workspace. ws=NULL is directly
!     initialized during the first call. At the end of the cycle
!     ws=FREE can be set or automatically expires at the end of
!     the program.
!     
!     r: ARRAY with the residuals of the iteration that we want
!     to accelerate. Can be substituted at the output in order
!     to accelerate the iteration.
!     
!     N: Size of the subspace.
      subroutine boostconv(rx_bc, ry_bc, rz_bc, rt_bc)
      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'
      include 'MASS_DEF'
      include 'MASS'            ! BM1
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP
      include 'INPUT_DEF'
      include 'INPUT'           ! ifheat    
      include 'BOOSTCONV'

!     The in_bc, out_bc and r arrays are distributed arrays

!     argument list
      real rx_bc(LX1,LY1,LZ1,LELV) 
      real ry_bc(LX1,LY1,LZ1,LELV)
      real rz_bc(LX1,LY1,LZ1,LELV)
      real rt_bc(LX1,LY1,LZ1,LELT)
      real ones(LX1,LY1,LZ1,LELT)
!     local variables
      integer rot_bc
      save rot_bc
      logical initialized
      save initialized
      data initialized /.false./

      integer iclean
      
      integer irbc(N_BC),icbc(N_BC)

      real dd_bc(N_BC,N_BC)
      save dd_bc    

      real dd_fact(N_BC,N_BC)
      real c_bc(N_BC)

      integer ipiv_bc(N_BC) 
      
      integer n_v,n_t,n_sub,j, info_bc, ntot1
!     functions
c     integer mod
      real op_glsc2_wt

!     Prabal
      integer i      
      
      n_v  = NX1*NY1*NZ1*NELV
      n_t  = NX1*NY1*NZ1*NELT
      n_sub= N_BC*N_BC

      call rone(ones,n_t)

      if (ISTEP.eq.BC_FREQ) initialized=.false.       ! Prabal. mod?
      
!     Executable statements
      if (.not.initialized) then

         do iclean = 1,N_BC
            call opzero(VINX_BC(1,1,1,1,iclean)
     $           ,VINY_BC(1,1,1,1,iclean),VINZ_BC(1,1,1,1,iclean))
            if (ifheat) call rzero(TIN_BC(1,1,1,1,iclean),n_t)
            call rzero(VOUTX_BC(1,1,1,1,iclean)
     $           ,VOUTY_BC(1,1,1,1,iclean),VOUTZ_BC(1,1,1,1,iclean))
            if (ifheat) call rzero(TOUT_BC(1,1,1,1,iclean),n_t)
         enddo
         call rzero(dd_bc,n_sub)
         
!        The dd_bc matrix
         do j = 1, N_BC
            dd_bc(j,j) = 1
         enddo
!        The index
         rot_bc = 1
         
         call opcopy(VOUTX_BC(1,1,1,1,rot_bc),
     $        VOUTY_BC(1,1,1,1,rot_bc),
     $        VOUTZ_BC(1,1,1,1,rot_bc),
     $        rx_bc,ry_bc,rz_bc)
         if (ifheat) call copy(TOUT_BC(1,1,1,1,rot_bc),rt_bc,n_t)

         call opcopy(VINX_BC(1,1,1,1,rot_bc),
     $        VINY_BC(1,1,1,1,rot_bc),
     $        VINZ_BC(1,1,1,1,rot_bc),
     $        rx_bc,ry_bc,rz_bc)
         if (ifheat) call copy(TIN_BC(1,1,1,1,rot_bc),rt_bc,n_t)

         initialized = .true.
         
      else
!        v_n = r_n-1 - r_n
         call opsub2 (VOUTX_BC(1,1,1,1,rot_bc),
     $        VOUTY_BC(1,1,1,1,rot_bc),
     $        VOUTZ_BC(1,1,1,1,rot_bc),
     $        rx_bc,ry_bc,rz_bc)
         if (ifheat) call sub2(TOUT_BC(1,1,1,1,rot_bc),rt_bc,n_t)

!        u_n - v_n = psi_n-1 - v_n
         call opsub2 (VINX_BC(1,1,1,1,rot_bc),
     $        VINY_BC(1,1,1,1,rot_bc),
     $        VINZ_BC(1,1,1,1,rot_bc),
     $        VOUTX_BC(1,1,1,1,rot_bc),
     $        VOUTY_BC(1,1,1,1,rot_bc),
     $        VOUTZ_BC(1,1,1,1,rot_bc))
         if (ifheat) call sub2(TIN_BC(1,1,1,1,rot_bc),
     $                   TOUT_BC(1,1,1,1,rot_bc),n_t)
!        D_ki = v_k * v_i

         if (ifheat) then
!            Prabal.
!            Needs the conjugate heat transfer routines.
!             do j = 1, N_BC
!                dd_bc(rot_bc,j)=cht_glsc2_wt(VOUTX_BC(1,1,1,1,rot_bc),
!     $               VOUTY_BC(1,1,1,1,rot_bc),
!     $               VOUTZ_BC(1,1,1,1,rot_bc),
!     $               TOUT_BC(1,1,1,1,rot_bc),
!     $               VOUTX_BC(1,1,1,1,j),
!     $               VOUTY_BC(1,1,1,1,j),
!     $               VOUTZ_BC(1,1,1,1,j),
!     $               TOUT_BC(1,1,1,1,j),BM1) ! (out(rot),out(j))
!                dd_bc(j,rot_bc) = dd_bc(rot_bc,j) ! make it symmetric
!             enddo
!!            t_k = v_k * y = v_k * r_n
!             do j = 1, N_BC 
!                c_bc(j)=cht_glsc2_wt(rx_bc,ry_bc,rz_bc,rt_bc,
!     $               VOUTX_BC(1,1,1,1,j),
!     $               VOUTY_BC(1,1,1,1,j),
!     $               VOUTZ_BC(1,1,1,1,j),
!     $               TOUT_BC(1,1,1,1,j),BM1) ! (r,out(j))
!             enddo
         else 
           do j = 1, N_BC
              dd_bc(rot_bc,j)=op_glsc2_wt(VOUTX_BC(1,1,1,1,rot_bc),
     $             VOUTY_BC(1,1,1,1,rot_bc),
     $             VOUTZ_BC(1,1,1,1,rot_bc),
     $             VOUTX_BC(1,1,1,1,j),
     $             VOUTY_BC(1,1,1,1,j),
     $             VOUTZ_BC(1,1,1,1,j),BM1) ! (out(rot),out(j))
              dd_bc(rot_bc,j)=dd_bc(rot_bc,j)/VOLVM1 ! Divide by volume
              dd_bc(j,rot_bc)=dd_bc(rot_bc,j)        ! make it symmetric
           enddo
!          t_k = v_k * y = v_k * r_n
           do j = 1, N_BC 
              c_bc(j)=op_glsc2_wt(rx_bc,ry_bc,rz_bc,
     $             VOUTX_BC(1,1,1,1,j),
     $             VOUTY_BC(1,1,1,1,j),
     $             VOUTZ_BC(1,1,1,1,j),BM1) ! (r,out(j))
              c_bc(j)=c_bc(j)/VOLVM1
           enddo
         endif          ! if (ifheat/ifflow)   

  102 format(A8,1x,5(E12.3E2,1x))

!        Prabal
!         if (nid.eq.0) then
!            do i=1,N_BC
!              write(6,102) 'dd_bc', (dd_bc(i,iclean), iclean=1,N_BC)
!            enddo
!            write(6,*) '-------------------------------'
!         endif   

         call copy(dd_fact,dd_bc,n_sub)
         call lu(dd_fact, N_BC, N_BC , irbc, icbc)
!        Prabal
!         if (nid.eq.0) then
!            do i=1,N_BC
!              write(6,102) 'dd_fact2',(dd_fact(i,iclean), iclean=1,N_BC)
!            enddo
!            write(6,*) '-------------------------------'
!         endif   
!!        Prabal
!         if (nid.eq.0) then
!            write(6,102) 'c_bc1', (c_bc(iclean), iclean=1,N_BC)
!         endif   

         call solve(c_bc,dd_fact,1,N_BC,N_BC,irbc,icbc)
!        Prabal
!         if (nid.eq.0) then
!            do i=1,N_BC
!              write(6,102) 'dd_fact3',(dd_fact(i,iclean), iclean=1,N_BC)
!            enddo
!            write(6,*) '-------------------------------'
!            write(6,102) 'c_bc2', (c_bc(iclean), iclean=1,N_BC)
!         endif   

!        increase the index
         rot_bc = mod(rot_bc, N_BC) + 1 ! 1 <= rot_bc <= n_bc

!        v_n+1 = r_n
         call opcopy(VOUTX_BC(1,1,1,1,rot_bc),
     $        VOUTY_BC(1,1,1,1,rot_bc),
     $        VOUTZ_BC(1,1,1,1,rot_bc),
     $        rx_bc,ry_bc,rz_bc)
         if (ifheat) call copy(TOUT_BC(1,1,1,1,rot_bc),rt_bc,n_t)

!        psi_n = r_n + sum_i(c_i (u_i - v_i))
         do j = 1, N_BC
            call opadd2cm(rx_bc,ry_bc,rz_bc,
     $           VINX_BC(1,1,1,1,j),
     $           VINY_BC(1,1,1,1,j),
     $           VINZ_BC(1,1,1,1,j),c_bc(j))
            if (ifheat) call add2s2(rt_bc,TIN_BC(1,1,1,1,j),c_bc(j),n_t)
         enddo
         call opcopy(VINX_BC(1,1,1,1,rot_bc),
     $        VINY_BC(1,1,1,1,rot_bc),
     $        VINZ_BC(1,1,1,1,rot_bc),
     $        rx_bc,ry_bc,rz_bc)
         if (ifheat) call copy(TIN_BC(1,1,1,1,rot_bc),rt_bc,n_t)
         
      endif
      
      return      
      end
!***********************************************************************
      subroutine boostconv_main

      implicit none
      include 'SIZE_DEF'
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP_DEF'
      INCLUDE 'TSTEP'
      INCLUDE 'BOOSTCONV'

      real coeff

c     Initial configuration
      if (istep.eq.0) then

         ! Get parameters from par file
!         call boostconv_param_get

         ! Fist order time integration required for boostconv
         if (IFBOOSTCONV) nbdinp = 1

        ! Set coefficient in cht routine
!         chcff_v = 1.0       ! hopefully not needed

         ! Set coeff variable (probably not necessary)
      endif

      coeff = 1.
c     Boost convergence to steady state
      if (istep .gt. 0 .and. IFBOOSTCONV) then
         if (BOOSTVER.eq.1) call boost_main(coeff)
!         if (BOOSTVER.eq.0) call nof_boost_main(coeff)
      endif         

c     Check for convergence
      if (coeff .lt. BC_TOL) then
         if (nio.eq.0) write(6,*) 'coeff: ', coeff
         lastep = 1
         nsteps = 0
      endif

      return
      end subroutine boostconv_main
!---------------------------------------------------------------------- 

      subroutine boost_main(norm2_bc)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! LX1, LY1, LZ1, LELT, NFIELD
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP, LASTEP
      include 'INPUT_DEF'       ! ifheat
      include 'INPUT'           ! IFUZAWA
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ], T
      include 'MASS_DEF'
      include 'MASS'            ! BM1
      include 'ADJOINT_DEF'
      include 'ADJOINT'         ! IFADJ
      include 'BOOSTCONV'       ! VXOLD,VYOLD,VZOLD,TOLD,VXOL1,VYOL1,
                                ! VZOL1,TOL1,DVX,DVY,DVZ,DTBC,BC_FREQ,
                                ! BC_TOL,IFBOOSTCONV

      real norm2_bc             ! residual norm (output)
      real dtboost
      save dtboost
      integer iswitch           ! Prabal. what is the point of this?
      integer nv,nt
      real op_glsc2_wt          ! function
      real opnorm2w             ! function

!     Prabal. Why didn't we use scratch arrays?      
      real vx_c(LX1,LY1,LZ1,LELV) 
      real vy_c(LX1,LY1,LZ1,LELV)
      real vz_c(LX1,LY1,LZ1,LELV)
      real t_c(LX1,LY1,LZ1,LELT)

!      logical IFUZAWA ! nof 18/8/2017 do not want to add it to INPUT
      integer n_v,n_t
!-----------------------------------------------------------------------
      
      n_v  = NX1*NY1*NZ1*NELV
      n_t  = NX1*NY1*NZ1*NELT

      if (ISTEP.EQ.1) then
         call opzero(VXOLD,VYOLD,VZOLD)
         if (ifheat) call rzero(TOLD,n_t)
         call opzero(VXOL1,VYOL1,VZOL1)
         if (ifheat) call rzero(TOL1,n_t)
         if (IFBOOSTVERBOSE) then
            if (nid.eq.0) then
!              first iteration: processor 0 (nid.eq.0) opens file
               open(UNIT=201,file='norm_dvdt.data',status='unknown',
     &              form='formatted')
            endif
         endif
      endif
      
      iswitch=mod(ISTEP,BC_FREQ)
      dtboost=dtboost+DT
      if (nid.eq.0) write(6,*) 'DTboost =',dtboost
      if (nid.eq.0) write(6,*) 'iswitch =',iswitch
      if (nid.eq.0) write(6,*) 'TOLboost=',BC_TOL
         
!     Calling Boostconv
      if(IFBOOSTCONV.and.(ISTEP.gt.1).and.(iswitch.eq.0)) then
         
         if (nid.eq.0) then
            write(6,*) '---------> calling boostconv '
         endif
         
!        For linear(adj) and not
         if(.not.IFADJ) then    !non-linear
            call opsub3(DVX,DVY,DVZ,VX,VY,VZ,
     $                  VXOL1,VYOL1,VZOL1) ! dv = v - vold1
            if (ifheat) call sub3(DTBC,T,TOL1,n_t)
         else                   !linear
            call opsub3(DVX,DVY,DVZ,VXP,VYP,VZP,
     $                  VXOL1,VYOL1,VZOL1) ! dv = v - vold1
            if (ifheat) call sub3(DTBC,TP,TOL1,n_t)
         endif
         
         call Boostconv(DVX,DVY,DVZ,DTBC)
         
!        For linear(adj) and not        
         if(.not.IFADJ) then    !non-linear
!           v = vold1 + dv 
            if (IFFLOW) then
               call add3(VX,VXOL1,DVX,n_v)
               call add3(VY,VYOL1,DVY,n_v)
               if(IF3D) call add3(VZ,VZOL1,DVZ,n_v)
            endif
            if (IFHEAT) call add3(T,TOL1,DTBC,n_t)
         else                   !linear 
            if (IFFLOW) then
               call add3(VXP,VXOL1,DVX,n_v)
               call add3(VYP,VYOL1,DVY,n_v)
               if(IF3D) call add3(VZP,VZOL1,DVZ,n_v)
            endif
            if (IFHEAT) call add3(TP,TOL1,DTBC,n_t)
         endif  
         
         IFUZAWA = .false.       ! why do we need this flag?
         
      endif                     ! end call to Boostconv
      
      
!     Check convergence
      if(IFBOOSTCONV.and.(iswitch.eq.0)) then !every Boostconv Dt
         
!        save the new ol1
         if(.not.IFADJ) then    !non-linear
            call opcopy(VXOL1,VYOL1,VZOL1,VX,VY,VZ)
            if (ifheat) call copy(TOL1,T,n_t)
         else                   !linear
            call opcopy(VXOL1,VYOL1,VZOL1,VXP,VYP,VZP)
            if (ifheat) call copy(TOL1,TP,n_t)
         endif
         dtboost=0
      else                      !every Dt
!        subtracting the new velocity field from the old one:
!        result stored in v[xyz]old
         if(.not.IFADJ) then    !non-linear
            call opsub2(VXOLD,VYOLD,VZOLD,VX,VY,VZ)
            if (ifheat) call sub2(TOLD,T,n_t)
         else                   !linear
            call opsub2(VXOLD,VYOLD,VZOLD,VXP,VYP,VZP)
            if (ifheat) call sub2(TOLD,TP,n_t)
         endif
         
!        norm of the difference
         if (ifheat) then
!           Prabal.
!            Needs the conjugate heat transfer routine
!            Not just an inner product. It is a norm.
!            Therefore divided by the volume.
!            norm2_bc=cht_glsc2_wt(VXOLD,VYOLD,VZOLD,TOLD
!     $           ,VXOLD,VYOLD,VZOLD,TOLD,BM1)
         else   
            norm2_bc=opnorm2w(VXOLD,VYOLD,VZOLD,BM1)
         endif
!         norm2_bc=sqrt(norm2_bc)   ! opnorm2bc already has sqrt.
         
!        normalization by DT (estimate of the norm of time derivative)
         norm2_bc=norm2_bc/DT
!        updating v[xyz]old with v[xyz] at present iteration
         if(.not.IFADJ) then    !non-linear
            call opcopy(VXOLD,VYOLD,VZOLD,VX,VY,VZ)
            if (ifheat) call copy(TOLD,T,n_t)
         else                   !linear
            call opcopy(VXOLD,VYOLD,VZOLD,VXP,VYP,VZP)
            if (ifheat) call copy(TOLD,TP,n_t)
         endif
         
         if (ISTEP.gt.1.and.iswitch.gt.2.and.IFBOOSTVERBOSE) then
!           write to file except first iteration
            if(nid.eq.0) write(201,101) ISTEP,TIME,norm2_bc,BC_TOL
  101       format(I5,1x,3(E13.5E2,1x))    
            call FLUSH(201)
         endif                  !istep.eq.0
         
      endif                     !if boost step or normal step
      
!     Exitt if l2 < tol in the FD_main loop
      return      
      end
!***********************************************************************
