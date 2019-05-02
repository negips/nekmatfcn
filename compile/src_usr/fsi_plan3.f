!======================================================================
!     Routines for fluid structure interaction implicit solutions
!     Author: Prabal S. Negi
!     Description: Solver routines for the AMP Stokes' solve.
!
!     Based on Fischer P., Schmitt M. and Tomboulides A. (2017) Recent
!     Developments in Spectral Element Simulations of Moving-Domain
!     Problems.
!======================================================================       

c-----------------------------------------------------------------------
      subroutine fsi_plan3
!     Compute pressure and velocity using consistent approximation spaces.     
!     Operator splitting technique.
!
!----------------------------------------------------------------------

      implicit none

      include 'SIZE_DEF'            
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'FSI'
      include 'FSI_DEBUG'
C

      real resv1,resv2,resv3,dv1,dv2,dv3

      COMMON /SCRNS/  RESV1 (LX1,LY1,LZ1,LELV)
     $ ,              RESV2 (LX1,LY1,LZ1,LELV)
     $ ,              RESV3 (LX1,LY1,LZ1,LELV)
     $ ,              DV1   (LX1,LY1,LZ1,LELV)
     $ ,              DV2   (LX1,LY1,LZ1,LELV)
     $ ,              DV3   (LX1,LY1,LZ1,LELV)

      real H1,H2
      COMMON /SCRVH/  H1    (LX1,LY1,LZ1,LELV)
     $ ,              H2    (LX1,LY1,LZ1,LELV)

      logical fsi_hlm

      integer icalld
      save icalld
      data icalld /0/

      real dttol
      parameter (dttol=1.0E-12)  ! Tolerance to check if DT has changed

      integer intype
      real mstep
      real tolv_old,tolp_old
      real tolv_new,tolp_new
      real tolv_max,tolp_max

!     In perturbation mode we don't change the geometry.
!     So we can calculate AMP arrays just once and save them      
!     Ideally don't need to go till 4th step.
      if (ifpert.and.istep.gt.4) then
!       If DT hasn't changed either                
        if (abs(DT-DTLAG(1)).lt.dttol) then
          return
        endif  
      endif  

      tolv_old = param(22)
      tolp_old = param(21)

      tolv_max = 1.0e-3
      tolp_max = 1.0e-3

      if (istep.lt.5) then
        continue
      else
        if (abs(fsi_alpha).gt.1e-12) then    
          tolv_new = (1.0e-5)*param(22)/abs(fsi_alpha)    ! Prabal. Arbitrarily set for now
          tolp_new = (1.0e-5)*param(21)/abs(fsi_alpha)
          param(21) = max(tolp_new,param(21))
          param(22) = max(tolv_new,param(22))
          param(21) = min(tolp_max,param(21))
          param(22) = min(tolv_max,param(22))
        endif          
!        if (nio.eq.0) write(6,*) 'Tol v/p',param(21),param(22),
!     $ fsi_alpha

      endif  

      if (nio.eq.0) write(6,11) 'Solving Stokes'' step', fsi_timea,
     $       fsi_timei
   11 format(A21,1x,2(E12.5E2,1x))
C
!      IF (IGEOM.EQ.1) THEN
C
C        Old geometry
C        No forcing for AMP.
!         CALL MAKEF
         call opzero(amp_vx,amp_vy,amp_vz)
         call rzero(amp_pr,nx2*ny2*nz2*nelt)
         call opzero(RESV1,RESV2,RESV3)
             
!      ELSE
C
C        New geometry, new b.c.
C

         INTYPE = -1
         CALL SETHLM  (H1,H2,INTYPE)
         CALL FSI_CRESVIF (RESV1,RESV2,RESV3,H1,H2)

         mstep = 0 !abs(param(94))
!         if (param(94).ne.0. .and. istep.ge.mstep) then
!           Projections. No projection of AMP right now. 
!           call ophinvpr(dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxh)
!         else
!          Just solve Helmholz equation right now
           fsi_hlm=.true. 
           CALL FSI_OPHINV  (DV1,DV2,DV3,RESV1,RESV2,RESV3,
     $                       H1,H2,TOLHV,NMXH,FSI_HLM)

!         endif

         CALL OPADD2  (AMP_VX,AMP_VY,AMP_VZ,DV1,DV2,DV3)

!        prabal   
!         call opcopy(dr1,dr2,dr3,amp_vx,amp_vy,amp_vz)
!         call copy(drp,amp_pr,lx2*ly2*lz2*lelv)

         call fsi_incomprn(amp_vx,amp_vy,amp_vz,amp_pr)

!         call chkdiv(amp_vx,amp_vy,amp_vz)   

!      ENDIF      ! IGEOM.eq.1
C
      param(22) = tolv_old            ! Prabal. Arbitrarily set for now
      param(21) = tolp_old

         
      return
      end subroutine fsi_plan3
C---------------------------------------------------------------------


      subroutine FSI_CRESVIF (resv1,resv2,resv3,h1,h2)
C---------------------------------------------------------------------
C
C     Compute startresidual/right-hand-side in the velocity solver
C     For the added mass step            
C
C---------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'FSI'

      REAL           RESV1 (LX1,LY1,LZ1,LELV)
      REAL           RESV2 (LX1,LY1,LZ1,LELV)
      REAL           RESV3 (LX1,LY1,LZ1,LELV)
      REAL           H1    (LX1,LY1,LZ1,LELV)
      REAL           H2    (LX1,LY1,LZ1,LELV)
      COMMON /SCRUZ/ W1    (LX1,LY1,LZ1,LELV)
     $ ,             W2    (LX1,LY1,LZ1,LELV)
     $ ,             W3    (LX1,LY1,LZ1,LELV)

      common /cgeom/ igeom

      NTOT1 = NX1*NY1*NZ1*NELV
      NTOT2 = NX2*NY2*NZ2*NELV
!      if (igeom.eq.2) CALL LAGVEL

!     Set Bdry conditions 
      CALL FSI_BCDIRVC (AMP_VX,AMP_VY,AMP_VZ,
     $       amp_v1mask,amp_v2mask,amp_v3mask)

!     We don't care about stresses at the outflow boundary
!      IF (IFSTRS)  CALL BCNEUTR

!     Zero previous pressures? Not entirely sure about this.
!      call extrapp (pr,prlag)
!      call opgradt (resv1,resv2,resv3,pr)

!     There's no RHS forcing
!      CALL OPADD2  (RESV1,RESV2,RESV3,BFX,BFY,BFZ)

      CALL OPHX (W1,W2,W3,AMP_VX,AMP_VY,AMP_VZ,H1,H2)
 
      call opzero  (RESV1,RESV2,RESV3) 
      CALL OPSUB2  (RESV1,RESV2,RESV3,W1,W2,W3)

C
      RETURN
      END
!---------------------------------------------------------------------- 

      subroutine fsi_ophinv (out1,out2,out3,inp1,inp2,inp3,
     $             h1,h2,tolh,nmxi,fsi_hlm)
!----------------------------------------------------------------------
!
!     OUT = (H1*A+H2*B)-1 * INP  (implicit)
!
!----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      include 'FSI'

      REAL OUT1 (LX1,LY1,LZ1,LELV)
      REAL OUT2 (LX1,LY1,LZ1,LELV)
      REAL OUT3 (LX1,LY1,LZ1,LELV)
      REAL INP1 (LX1,LY1,LZ1,LELV)
      REAL INP2 (LX1,LY1,LZ1,LELV)
      REAL INP3 (LX1,LY1,LZ1,LELV)
      REAL H1   (LX1,LY1,LZ1,LELV)
      REAL H2   (LX1,LY1,LZ1,LELV)
C
      logical fsi_hlm


      IMESH = 1
C
      if (ifstrs) then
        MATMOD = 0
        if (fsi_hlm) then
           CALL FSI_HMHZSF  ('NOMG',OUT1,OUT2,OUT3,INP1,INP2,INP3,
     $                    H1,H2,
     $                    AMP_V1MASK,AMP_V2MASK,AMP_V3MASK,VMULT,
     $                    TOLH,NMXI,MATMOD)
        endif
      elseif (ifcyclic) then
        continue
      elseif (ifpert) then
        CALL HMHOLTZ ('VELX',OUT1,INP1,H1,H2,AMP_V1MASK,VMULT,
     $                                  IMESH,TOLH,NMXI,1)
        CALL HMHOLTZ ('VELY',OUT2,INP2,H1,H2,AMP_V2MASK,VMULT,
     $                                  IMESH,TOLH,NMXI,2)
        IF (NDIM.EQ.3) 
     $  CALL HMHOLTZ ('VELZ',OUT3,INP3,H1,H2,AMP_V3MASK,VMULT,
     $                                  IMESH,TOLH,NMXI,3)

      else
        if (nio.eq.0) write(6,*) 'OPHINV not defined for AMP module'   
        continue   
      endif
C
      return
      end subroutine fsi_ophinv

!----------------------------------------------------------------------
      subroutine fsi_incomprn (ux,uy,uz,up)
c
c     Project U onto the closest incompressible field
c
c     Input:  U     := (ux,uy,uz)
c
c     Output: updated values of U, iproj, proj; and
c             up    := pressure currection req'd to impose div U = 0
c
c
c     Dependencies: ifield ==> which "density" (vtrans) is used.
c
c     Notes  1.  up is _not_ scaled by bd(1)/dt.  This should be done
c                external to incompr().
c
c            2.  up accounts _only_ for the perturbation pressure,
c                not the current pressure derived from extrapolation.
c
c

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'FSI_DEBUG'
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

!      parameter(nset = 1 + lbelv/lelv)
!      common /orthov/ pset(lx2*ly2*lz2*lelv*mxprev,nset)
!      common /orthbi/ nprv(2)

      logical ifprjp

!     If Project out previous pressure solutions?
      if (istep.ge.10) then
        ifprjp=.true.
      else
        Nprev=0
        ifprjp=.false.
      endif

!      if (icalld.eq.0) tpres=0.0
!      icalld = icalld+1
!      npres  = icalld
      etime1 = dnekclock()

      ntot1  = nx1*ny1*nz1*nelv
      ntot2  = nx2*ny2*nz2*nelv
      intype = 1

      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call invers2 (h2inv,h2,ntot1)

      call opdiv   (dp,ux,uy,uz)

      bdti = -bd(1)/dt
      call cmult   (dp,bdti,ntot2)

!      call add2col2(dp,bm2,usrdiv,ntot2) ! User-defined divergence.

      call ortho   (dp)

      if (ifprjp)  call fsi_setrhsp(dp,h1,h2,h2inv)
                   scaledt = dt/bd(1)
                   scaledi = 1./scaledt
                   call cmult(dp,scaledt,ntot2)        ! scale for tol
!                   call opcopy(dr1,dr2,dr3,ux,uy,uz)   ! prabal. debugging
!                   call copy(drp,dp,ntot2)             ! prabal. debugging
                   call esolver  (dp,h1,h2,h2inv,intype)
                   call cmult(dp,scaledi,ntot2)
      if (ifprjp)  call fsi_gensolnp(dp,h1,h2,h2inv)

      call add2(up,dp,ntot2)

      call opgradt  (w1 ,w2 ,w3 ,dp)
      call opbinv   (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
      dtb  = dt/bd(1)
      call opadd2cm (ux ,uy ,uz ,dv1,dv2,dv3, dtb )

      tpres=tpres+(dnekclock()-etime1)

      return
      end subroutine fsi_incomprn
c-----------------------------------------------------------------------

      subroutine fsi_bcdirvc(v1,v2,v3,mask1,mask2,mask3)
C
C     Apply Dirichlet boundary conditions to surface of vector (V1,V2,V3).
C     Use IFIELD as a guide to which boundary conditions are to be applied.
C
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'TOPOL'
      INCLUDE 'CTIMER'
      COMMON /SCRUZ/ TMP1(LX1,LY1,LZ1,LELV)
     $             , TMP2(LX1,LY1,LZ1,LELV)
     $             , TMP3(LX1,LY1,LZ1,LELV)
      COMMON /SCRMG/ TMQ1(LX1,LY1,LZ1,LELV)
     $             , TMQ2(LX1,LY1,LZ1,LELV)
     $             , TMQ3(LX1,LY1,LZ1,LELV)
C
      REAL V1(NX1,NY1,NZ1,LELV),V2(NX1,NY1,NZ1,LELV)
     $    ,V3(NX1,NY1,NZ1,LELV)
      real mask1(nx1,ny1,nz1,lelv),mask2(nx1,ny1,nz1,lelv)
     $    ,mask3(nx1,ny1,nz1,lelv)
c
      common  /nekcb/ cb
      character cb*3
      character*1 cb1(3)
      equivalence (cb1,cb)
c
      logical ifonbc

      real val                ! value at the faces
c
      etime1=dnekclock()
C
C
      NFACES=2*NDIM
      NXYZ  =NX1*NY1*NZ1
      NEL   =NELFLD(IFIELD)
      NTOT  =NXYZ*NEL
C
      CALL RZERO(TMP1,NTOT)
      CALL RZERO(TMP2,NTOT)
      IF (IF3D) CALL RZERO(TMP3,NTOT)
C
C     Velocity boundary conditions
C
c     write(6,*) 'BCDIRV: ifield',ifield
      DO 2100 ISWEEP=1,2
         DO 2000 IE=1,NEL
         DO 2000 IFACE=1,NFACES
            CB  = CBC(IFACE,IE,IFIELD)
            BC1 = BC(1,IFACE,IE,IFIELD)
            BC2 = BC(2,IFACE,IE,IFIELD)
            BC3 = BC(3,IFACE,IE,IFIELD)

            IF (CB.EQ.'ON ' .OR. CB.EQ.'on ') then   ! 5/21/01 pff
                ifonbc =.true.
                CALL FSI_FACEIV ('v  ',TMP1(1,1,1,IE),TMP2(1,1,1,IE),
     $                       TMP3(1,1,1,IE),IE,IFACE,NX1,NY1,NZ1)
            else 
                call fsi_faceiv (cb,tmp1(1,1,1,ie),tmp2(1,1,1,ie),
     $                       tmp3(1,1,1,ie),ie,iface,nx1,ny1,nz1)

                IF ( IFQINP(IFACE,IE) )
     $          CALL GLOBROT (TMP1(1,1,1,IE),TMP2(1,1,1,IE),
     $                        TMP3(1,1,1,IE),IE,IFACE)
            ENDIF

 2000    CONTINUE

C
C        Take care of Neumann-Dirichlet shared edges...
C
         if (isweep.eq.1) then
            call opdsop(tmp1,tmp2,tmp3,'MXA')
         else
            call opdsop(tmp1,tmp2,tmp3,'MNA')
         endif
 2100 CONTINUE
C
C     Copy temporary array to velocity arrays.
C
      IF ( .NOT.IFSTRS ) THEN
         CALL COL2(V1,mask1,NTOT)
         CALL COL2(V2,mask2,NTOT)
         IF (IF3D) CALL COL2(V3,mask3,NTOT)
         if (ifonbc) then
            call antimsk1(tmp1,mask1,ntot)
            call antimsk1(tmp2,mask2,ntot)
            if (if3d) call antimsk1(tmp3,mask3,ntot)
         endif
      ELSE
         IF (IFMODEL) THEN
             CALL COPY (TMQ1,TMP1,NTOT)
             CALL COPY (TMQ2,TMP2,NTOT)
             IF (NDIM.EQ.3) CALL COPY (TMQ3,TMP3,NTOT)
             CALL AMASK (TMP1,TMP2,TMP3,TMQ1,TMQ2,TMQ3,NELV)
         ENDIF
         CALL FSI_RMASK (V1,V2,V3,NELV)
      ENDIF

      CALL ADD2(V1,TMP1,NTOT)
      CALL ADD2(V2,TMP2,NTOT)
      IF (IF3D) CALL ADD2(V3,TMP3,NTOT)

      tusbc=tusbc+(dnekclock()-etime1)

      return
      end
c-----------------------------------------------------------------------

      subroutine fsi_faceiv (cb,v1,v2,v3,iel,iface,nx,ny,nz)
C
C     Assign fortran function boundary conditions to 
C     face IFACE of element IEL for vector (V1,V2,V3).
C
      INCLUDE 'SIZE'
      INCLUDE 'NEKUSE'
      INCLUDE 'PARALLEL'
C
      dimension v1(nx,ny,nz),v2(nx,ny,nz),v3(nx,ny,nz)
      character cb*3
c
      character*1 cb1(3)
c
      common  /nekcb/ cb3
      character*3 cb3
      cb3 = cb
c
      call chcopy(cb1,cb,3)
c
      ieg = lglel(iel)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
C
      IF (CB.EQ.'mv '.OR.CB.EQ.'mvn'.or.CB.EQ.'v  ') THEN
C
         DO 100 IZ=KZ1,KZ2
         DO 100 IY=KY1,KY2
         DO 100 IX=KX1,KX2
            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL FSI_USERBC  (IX,IY,IZ,IFACE,IEG)
            V1(IX,IY,IZ) = UX
            V2(IX,IY,IZ) = UY
            V3(IX,IY,IZ) = UZ
  100    CONTINUE
         RETURN
C
      else 
         do iz=kz1,kz2
         do iy=ky1,ky2
         do ix=kx1,kx2
            V1(IX,IY,IZ) = 0.
            V2(IX,IY,IZ) = 0.
            V3(IX,IY,IZ) = 0.
         enddo
         enddo
         enddo
         return
      ENDIF

C
      return
      end subroutine fsi_faceiv
c-----------------------------------------------------------------------

      subroutine fsi_userbc (ix,iy,iz,iside,ieg)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKUSE_DEF'
      include 'NEKUSE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'FSI'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'

      integer ix,iy,iz,ieg,iside,iel
      real dx,dy

      iel=gllel(ieg)

      if (ifpert) then      
!        if ((abs(x).lt.0.6).and.(abs(y).lt.0.6)) then  ! for a cylinder
        if (cbu.eq.'mv ') then 
          if (.not.fsi_ifrot) then
            ux = 0. 
            uy = etav_g 
            if (ndim.eq.3) uz = 0.
          else
!           We solve for rotational equations.
!           Meaning velocity components depend on current position (eta)
!           as well as xm1,ym1 distance from the rotational axis.
!           The coordinate system assumes positive rotational velocity creates a
!           clockwise motion.
            dx = (x-fsi_x0)
            dy = (y-fsi_y0)

            ux =  etav_g*dy
            uy = -etav_g*dx
            if (ndim.eq.3) uz = 0.
          endif               ! fsi_ifrot

        else            ! zero at external boundaries
          ux = 0.
          uy = 0.
          uz = 0.
        endif
      else              ! non-linear simulation
        if (cbu.eq.'mv ') then 
          if (.not.fsi_ifrot) then
            ux = 0. 
            uy = etav_g 
            if (ndim.eq.3) uz = 0.
          else
!           We solve for rotational equations.
!           Meaning velocity components depend on distance of points 
!           from the rotational axis.
!           The coordinate system assumes positive rotational velocity creates a
!           clockwise motion.
            dx = (x-fsi_x0)
            dy = (y-fsi_y0)

            ux =  etav_g*dy
            uy = -etav_g*dx
            if (ndim.eq.3) uz = 0.
          endif 

        else
          ux = 0.
          uy = 0.
          uz = 0.
        endif
      endif             ! ifpert

      return
      end
c-----------------------------------------------------------------------

      subroutine set_amp_mask
C
C     Zero out masks for all boundary points.
C     In the AMP solve all boundary conditions are set to 0.
C     Except Periodic boundary conditions
C

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'
      include 'SOLN_DEF'
      include 'SOLN'

      integer e,f
      integer nfaces
      integer nxyz

      nfaces=2*ndim
      nxyz  =nx1*ny1*nz1

!      call oprone(amp_v1mask,amp_v2mask,amp_v3mask)
      call rone(amp_v1mask,nxyz*nelv)
      call rone(amp_v2mask,nxyz*nelv)
      if (ndim.eq.3) call rone(amp_v3mask,nxyz*nelv)

      do e=1,nelv
        do f=1,nfaces
          if (cbc(f,e,1).ne.'E  '.and.cbc(f,e,1).ne.'P  '.and.
     $      cbc(f,e,1).ne.'o  '.and.cbc(f,e,1).ne.'O  ') then
            call facev(amp_v1mask,e,f,0.,nx1,ny1,nz1)     
            call facev(amp_v2mask,e,f,0.,nx1,ny1,nz1)     
            if (ndim.eq.3) call facev(amp_v3mask,e,f,0.,nx1,ny1,nz1)
          endif
        enddo
      enddo

      call opdsop(amp_v1mask,amp_v2mask,amp_v3mask,'MUL') ! no rotation for mul

      return
      end subroutine set_amp_mask

!---------------------------------------------------------------------- 

      subroutine fsi_rmask (r1,r2,r3,nel)
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'FSI'
C
      dimension r1  (lx1,ly1,lz1,nelt)
     $        , r2  (lx1,ly1,lz1,nelt)
     $        , r3  (lx1,ly1,lz1,nelt)
C

      call qmask (r1,r2,r3,amp_v1mask,amp_v2mask,amp_v3mask,nel)


      return
      end
!---------------------------------------------------------------------- 
      subroutine fsi_hmhzsf (name,u1,u2,u3,r1,r2,r3,h1,h2,
     $                   rmask1,rmask2,rmask3,rmult,
     $                   tol,maxit,matmod)
C-----------------------------------------------------------------------
C
C     Compute solution to coupled Helmholtz equations 
C     (stress formulation)
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
!      include 'SOLN'   ! For outpost diagnostic call
      include 'TSTEP'
      include 'ORTHOSTRS'

      DIMENSION U1(LX1,LY1,LZ1,LELV)
     $        , U2(LX1,LY1,LZ1,LELV)
     $        , U3(LX1,LY1,LZ1,LELV)
     $        , R1(LX1,LY1,LZ1,LELV)
     $        , R2(LX1,LY1,LZ1,LELV)
     $        , R3(LX1,LY1,LZ1,LELV)
     $        , H1(LX1,LY1,LZ1,LELV)
     $        , H2(LX1,LY1,LZ1,LELV)
     $        , RMASK1(LX1,LY1,LZ1,LELV)
     $        , RMASK2(LX1,LY1,LZ1,LELV)
     $        , RMASK3(LX1,LY1,LZ1,LELV)
     $        , RMULT (LX1,LY1,LZ1,LELV)
      CHARACTER NAME*4

      common /cpfjunk/ y(lx1*ly1*lz1*lelt,3)
      common /cpfjun2/ v(lx1*ly1*lz1*lelt,3)

      nel = nelfld(ifield)
      vol = volfld(ifield)
      n   = nx1*ny1*nz1*nel

      call fsi_rmask   (r1,r2,r3,nel)
      call opdssum (r1,r2,r3)
      call rzero3  (u1,u2,u3,n)

c     call set_up_h1_crs_strs(h1,h2,ifield,matmod)
      
      if (imesh.eq.1) then
         call chktcgs (r1,r2,r3,rmask1,rmask2,rmask3,rmult,binvm1
     $                ,vol,tol,nel)

!         if (matmod.lt.0) then
!          napprox(1) = 0
!          iproj      = 0 !param(94)
!          if (iproj.gt.0.and.istep.gt.iproj) napprox(1)=param(93)
!          napprox(1)=min(napprox(1),istep/3)
!          call strs_project_a(r1,r2,r3,h1,h2,rmult,ifield,ierr,matmod)

c         call opcopy(y(1,1),y(1,2),y(1,3),x(1),x(1+n),x(1+2*n))

!         endif

         call cggosf  (u1,u2,u3,r1,r2,r3,h1,h2,rmult,binvm1
     $                ,vol,tol,maxit,matmod)


!         if (matmod.lt.0)
!     $    call strs_project_b(u1,u2,u3,h1,h2,rmult,ifield,ierr)

      else
         call chktcgs (r1,r2,r3,rmask1,rmask2,rmask3,rmult,bintm1
     $                ,vol,tol,nel)
         call cggosf  (u1,u2,u3,r1,r2,r3,h1,h2,rmult,bintm1
     $                ,vol,tol,maxit,matmod)
      endif

      return
      end
!---------------------------------------------------------------------- 
      subroutine fsi_setrhsp(p,h1,h2,h2inv)
C
C     Project soln onto best fit in the "E" norm.
C
      implicit none
      
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'MASS_DEF'
      include 'MASS'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'FSI'

      real p    (lx2,ly2,lz2,lelv)
      real h1   (lx1,ly1,lz1,lelv)
      real h2   (lx1,ly1,lz1,lelv)
      real h2inv(lx1,ly1,lz1,lelv)
!      real pset (lx2,ly2,lz2,lelv,amp_mproj)

      integer ltot2
      parameter (ltot2=lx2*ly2*lz2*lelv)

      integer i,intetype,ntot2,ierr

      real vlsc2                ! function

!      real pbar,pnew,alpha,work
!      common /fsi_orthox/ pbar(ltot2),pnew(ltot2)
!      common /fsi_orthos/ alpha(amp_mproj),work(amp_mproj)

      if (Nprev.eq.0) return

      ntot2  = nx2*ny2*nz2*nelv

      ierr = 0
      call fsi_updrhse(p,h1,h2,h2inv,ierr) ! update rhs's if E-matrix has changed
      if (ierr.ne.0) then
        if (nio.eq.0) then
          write(6,*) 'Error in fsi_updrhse. Setting Nprev=0' 
        endif  
        Nprev=0           ! Doesn't happen w/ new formulation
      endif  

      do i=1,Nprev  ! Perform Gram-Schmidt for previous soln's.
         alpha(i) = vlsc2(p,amp_prproj(1,1,1,1,i),ntot2)
      enddo
      call gop(alpha,work,'+  ',Nprev)

      call rzero(pbar,ntot2)
      do i=1,Nprev
         call add2s2(pbar,amp_prproj(1,1,1,1,i),alpha(i),ntot2)
      enddo
C
      intetype = 1
      call cdabdtp(pnew,pbar,h1,h2,h2inv,intetype)
      call sub2   (p,pnew,ntot2)

      return
      end subroutine fsi_setrhsp
c-----------------------------------------------------------------------
      subroutine fsi_gensolnp(p,h1,h2,h2inv)
C
C     Reconstruct the solution to the original problem by adding back
C     the previous solutions
C
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'

      real p    (lx2,ly2,lz2,lelv)
      real h1   (lx1,ly1,lz1,lelv)
      real h2   (lx1,ly1,lz1,lelv)
      real h2inv(lx1,ly1,lz1,lelv)
!      real pset (lx2*ly2*lz2*lelv,amp_mproj)

      integer ltot2
      parameter (ltot2=lx2*ly2*lz2*lelv)

!      real pbar,pnew,alpha,work
!      common /fsi_orthox/ pbar(ltot2),pnew(ltot2)
!      common /fsi_orthos/ alpha(amp_mproj),work(amp_mproj)

      integer ierr,ntot2,mprv

      integer isave,isave_freq
      save isave
      data isave /0/

      ierr = 0

      mprv=amp_mproj

      isave_freq=1
      isave=mod(isave,isave_freq)

      ntot2=nx2*ny2*nz2*nelv

      if (Nprev.lt.mprv) then
        if (isave.eq.0) then
!          if (nid.eq.0) write(6,*) 'Saving new AMP projection'
          Nprev = Nprev + 1
          call copy  (amp_prproj(1,1,1,1,Nprev),p,ntot2)          ! Save current solution
          call econjp(amp_prproj,Nprev,h1,h2,h2inv,ierr)          ! Orthonormalize set
        endif  
        call add2  (p,pbar,ntot2)                               ! Reconstruct solution.

        if (ierr.eq.1) then
          if (nid.eq.0) write(6,*) 'Error in AMP projection', ierr
          isave = 0
          Nprev = 1
          call copy  (amp_prproj(1,1,1,1,Nprev),p,ntot2)        ! Save current solution
          call econjp(amp_prproj,Nprev,h1,h2,h2inv,ierr)        ! and orthonormalize.
        endif
      else                                                      ! (uses pnew).
!        if (nid.eq.0) write(6,*) 'Resetting AMP projections'
        Nprev = 1
        call add2  (p,pbar,ntot2)                               ! Reconstruct solution.
        call copy  (amp_prproj(1,1,1,1,Nprev),p,ntot2)          ! Save current solution
        call econjp(amp_prproj,Nprev,h1,h2,h2inv,ierr)          ! and orthonormalize.
      endif
      isave = isave + 1

      return
      end subroutine fsi_gensolnp
c-----------------------------------------------------------------------
      subroutine fsi_updrhse(p,h1,h2,h2inv,ierr)
C
C     Update rhs's if E-matrix has changed
C
C
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'FSI'
C
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
!      COMMON /FSI_ORTHOX/ Pbar(LTOT2),Pnew(LTOT2)
!      COMMON /FSI_ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
!      COMMON /FSI_ORTHOI/ Nprev,Mprev
      COMMON /ORTHOL/ IFNEWE
!      REAL ALPHA,WORK
      LOGICAL IFNEWE
C
C
      REAL             P    (LX2,LY2,LZ2,LELV)
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
C
      integer icalld
      save    icalld
      data    icalld/0/

      NTOT2=NX2*NY2*NZ2*NELV
C
C
C     First, we have to decide if the E matrix has changed.
C
      IF (icalld.eq.0) THEN
         icalld=1
         DTlast=DT
      ENDIF
C
      IFNEWE=.FALSE.
      IF (IFMVBD) THEN
         IFNEWE=.TRUE.
         CALL INVERS2(bm2inv,bm2,Ntot2)
      ELSEIF (DTlast.ne.DT) THEN
         IFNEWE=.TRUE.
         DTlast=DT
      ENDIF
      IF (IFNEWE.and.nio.eq.0) write(6,*) 'reorthogo AMP:',Nprev
C     
C     
C     Next, we reconstruct a new rhs set.
C     
      IF (IFNEWE) THEN
c
c        new idea...
c        if (nprev.gt.0) nprev=1
c        call copy(rhs,pnew,ntot2)
c
         Nprevt = Nprev
         DO 100 Iprev=1,Nprevt
C           Orthogonalize this rhs w.r.t. previous rhs's
            CALL fsi_econj (Iprev,H1,H2,H2INV,ierr)
            if (ierr.eq.1) then
               Nprev = 0
               return
            endif
  100    CONTINUE
C
      ENDIF
C
      return
      end subroutine fsi_updrhse
!---------------------------------------------------------------------- 

      subroutine fsi_econj(kprev,h1,h2,h2inv,ierr)
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
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'FSI'
C
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
C
      INTEGER LTOT2
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
!      COMMON /FSI_ORTHOX/ Pbar(ltot2),Pnew(ltot2),Pbrr(ltot2)
!      COMMON /FSI_ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
!      COMMON /FSI_ORTHOI/ Nprev,Mprev
!      REAL ALPHA,WORK
      real ALPHAd
      real alpham
      real GLSC2,VLSC2

      integer Kprev,Kprev1,ipass,npass,ntot2
      integer I,ierr
      integer intetype
C
C
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
         CALL CDABDTP(Pbrr,AMP_PRPROJ(1,1,1,1,Kprev),H1,H2,H2INV,
     $                  INTETYPE)
C
C        Compute part of the norm
         Alphad = GLSC2(AMP_PRPROJ(1,1,1,1,Kprev),Pbrr,NTOT2)
C
C        Gram-Schmidt
         Kprev1=Kprev-1
         DO 10 I=1,Kprev1
            ALPHA(I) = VLSC2(Pbrr,AMP_PRPROJ(1,1,1,1,i),NTOT2)
   10    CONTINUE
         IF (Kprev1.GT.0) CALL gop(alpha,WORK,'+  ',Kprev1)
C
         DO 20 I=1,Kprev1
            alpham = -alpha(i)
            CALL ADD2S2(AMP_PRPROJ(1,1,1,1,Kprev),AMP_PRPROJ(1,1,1,1,i),
     $             alpham,NTOT2)
            Alphad = Alphad - alpha(i)**2
   20    CONTINUE
      enddo
C
C    .Normalize new element in P~
C
      if (ALPHAd.le.0.0) then
         write(6,*) 'ERROR:  alphad .le. 0 in ECONJ',alphad,Kprev
         ierr = 1
         return
      endif
      ALPHAd = 1.0/SQRT(ALPHAd)
      ALPHAN = Alphad
      CALL CMULT(AMP_PRPROJ (1,1,1,1,Kprev),alphan,NTOT2)
C
      return
      end subroutine fsi_econj
!----------------------------------------------------------------------

