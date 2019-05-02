!======================================================================
!     Description: Routines for coupling with ARNOLDI
!     Author: Prabal S. Negi
!     Last Modified: 11-01-2019      
!======================================================================
!----------------------------------------------------------------------       
      real function fsi_glsc2_wt()

!     global weighted scalar product

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
      include 'SOLN_DEF'
      include 'SOLN'
      include 'MASS_DEF'
      include 'MASS'            ! VOLVM1, VOLTM1
      include 'FSI'             ! IFFSI

!     local variables
      integer ntotv, ntott
      real scal, f1, f2, f3, f4
      integer i
      real op_glsc2_wt
      real glsc3  

      integer nfsi

      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT
      nfsi  = NVECAFSI           ! nonzero only for nid=0

      scal = 0.
      scal = op_glsc2_wt(vxp,vyp,vzp,vxp,vyp,vzp,bm1)

!     Zero weight for pressure      

!     for conjugate heat transfer
      if (IFHEAT) then
        scal = scal + glsc3(tp,tp,bm1,ntott)
      endif  

!     Fluid-Structure-Interaction
      scal = scal + glsc3(eta,eta,fsi_eta_wt,nfsi)
      scal = scal + glsc3(etav,etav,fsi_etav_wt,nfsi)

      fsi_glsc2_wt = scal 
      
      return
      end function fsi_glsc2_wt
!---------------------------------------------------------------------- 

      subroutine fsi_opcmult (v1,v2,v3,t1,s1,s2,const)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
      include 'FSI'             ! IFFSI

!     argument list
      real v1(1),v2(1),v3(1),t1(1)
      real s1(1),s2(1)
      real const

!     local variables
      integer ntotv, ntott
      integer nfsi

      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT
      nfsi  = NVECAFSI

      if (IFFLOW) then
        call cmult(v1,const,ntotv)
        call cmult(v2,const,ntotv)
        if (IF3D) call cmult(v3,const,ntotv)
      endif

      if (IFHEAT) call cmult(t1,const,ntott)

      if (IFFSI)  then
         call cmult(s1,const,nfsi)
         call cmult(s2,const,nfsi)
      endif  

      return
      end subroutine fsi_opcmult

!---------------------------------------------------------------------- 

      subroutine fsi_weight_fcn

!     Called when IDOARP.eq.2
!     if y=Ax is the Arnoldi vector
!     This subroutine calculates w = My
!     and stores it in WORDA(IPNTARP(2))
!     Needed by Arnoldi to perform M-weighted
!     orthogonalization

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'FSI'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'
      include 'MASS_DEF'
      include 'MASS'                ! BM1
      include 'SOLN_DEF'
      include 'SOLN'                ! MASKs
      include 'INPUT_DEF'
      include 'INPUT'

      real bmv_wt(lx1,ly1,lz1,lelv)
      real bmt_wt(lx1,ly1,lz1,lelt)

      integer ntotv,ntott
      integer vecptr
      integer aptr

      ntotv = lx1*ly1*lz1*nelv
      ntott = lx1*ly1*lz1*nelt

      call copy(bmv_wt,BM1,ntotv)
      call copy(bmt_wt,BM1,ntott)

      vecptr = 0
      aptr   = IPNTARP(2)

!     velocity
      call copy(WORKDA(aptr),bmv_wt,NVECAV)
      vecptr=vecptr+NVECAV
      call copy(WORKDA(aptr+vecptr),bmv_wt,NVECAV)
      vecptr=vecptr+NVECAV
      if (IF3D) then
        call copy(WORKDA(aptr+vecptr),
     $     bmv_wt,NVECAV)
        vecptr=vecptr+NVECAV
      endif

!     Pressure. Zero weight. (if we are using a Semi-norm) 
      if (TSTIFPR) then
        call rzero(WORKDA(aptr+vecptr),NVECAP)
        vecptr=vecptr+NVECAP
      endif  

!     temperature
      if(IFHEAT) then
        call copy(WORKDA(aptr+vecptr),
     $       bmt_wt,NVECAT)
             vecptr=vecptr+NVECAT
      endif

      call copy(WORKDA(aptr+vecptr),fsi_eta_wt,NVECAFSI)
      vecptr=vecptr+NVECAFSI
      call copy(WORKDA(aptr+vecptr),fsi_etav_wt,NVECAFSI)
      vecptr=vecptr+NVECAFSI

      call col2(WORKDA(IPNTARP(2)),WORKDA(IPNTARP(1)),NVECAS)

      return
      end subroutine fsi_weight_fcn
!----------------------------------------------------------------------             

      subroutine fsi_arnoldi_restart

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'FSI'

!     Zero lag arrays      
      call rzero(etalag,lorder-1)
      call rzero(etavlag,lorder-1)

!     velocity/position extrapolation      
      etav_s = etav
      eta_s  = eta + etav_s*DT
      etaa   = 0.
      fsi_alpha = 0.

      call fsi_meshv_pert    ! set boundary conditions

      call fsi_output

      return
      end subroutine fsi_arnoldi_restart
!----------------------------------------------------------------------       

      subroutine fsi_arnoldi_rstsave

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'
      include 'CHKPOINT'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'

      integer i,is
      character*132 fname
      character*1   fname1(132)
      equivalence (fname1,fname)      

      character*1 exten0,exten1
      save exten0,exten1
      data exten0 /"0"/
      data exten1 /"1"/

      integer exten
      save exten
      data exten /0/

      character*1 dot
      save dot
      data dot /"."/
      character*16 outfmt

      integer iunit,ierr
      integer vecptr

      ierr=0
      fsi_rst_eta(1)  = eta
      fsi_rst_etav(1) = etav

      if (nid.eq.0) then
!       Create filename
        call blank(fname1,132)
        write(fname,'(A11,A1,i1)') fsi_rst_arnoldi,dot,ARNISTOP

!       Write output
        call IO_freeid(iunit,ierr)     
        open(iunit,file=fname,status='unknown',action='write',
     $    iostat=ierr)
        if (ierr.eq.0) then
          write(6,'(A24,1x,A13)') 'FSI ARN_RST: Saving file',
     $             fname
!!         ETA,ETAV  
!          write(iunit,*) (fsi_rst_eta(i),i=1,NVECAFSI)       ! eta
!          write(iunit,*) (fsi_rst_etav(i),i=1,NVECAFSI)      ! etav

!         RESIDA 
          write(iunit,*) (eta_resida(i),i=1,NVECAFSI)
          write(iunit,*) (etav_resida(i),i=1,NVECAFSI)

!         WORKDA
          write(iunit,*) (eta_workda(i),i=1,3*NVECAFSI)
          write(iunit,*) (etav_workda(i),i=1,3*NVECAFSI)

!         IPNTARP
          write(iunit,*) (IPNTARP(i),i=1,14)

          close(iunit)
          write(6,'(A4)')  'Done'
        endif   ! ierr.eq.0 
      endif       ! nid.eq.0

      call err_chk(ierr,'Error opening fsi_arn_rst file.')
      
      return 
      end subroutine fsi_arnoldi_rstsave

!---------------------------------------------------------------------- 

      subroutine fsi_arnoldi_rstread

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'FSI'
      include 'CHKPOINT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'
      
      character*1 exten0,exten1
      save exten0,exten1
      data exten0 /'0'/
      data exten1 /'1'/

      integer exten
      save exten
      data exten /0/

      character*132 fname
      character*1   fname1(132)
      equivalence (fname1,fname)      


      character*1 dot
      save dot
      data dot /'.'/
      character*16 outfmt

      integer ltrunc
      integer len

      real bd_etav
      integer ilag
      integer iunit,ierr

      integer workdapos(14)
      integer i,j,k

      integer vecptr0

      character*132 tmp


      ierr=0 
      if (istep.gt.0) return

!     Create filename
      call blank(fname1,132)
      write(fname,'(A11,A1,i1)') fsi_rst_arnoldi,dot,ARNISTART

!     Read restart file 
      if (nid.eq.0) then
        call IO_freeid(iunit,ierr)
        open(iunit, file=fname, status='old',action='read',
     $    iostat=ierr)
        if (ierr.eq.0) then
          write(6,'(A25,1x,A13)') 'FSI Restart: Opening file',
     $             fname
!!         ETA,ETAV
!          read(iunit,*) (fsi_rst_eta(i),i=1,NVECAFSI)
!          read(iunit,*) (fsi_rst_etav(i),i=1,NVECAFSI)

!         RESIDA 
          read(iunit,*) (eta_resida(i),i=1,NVECAFSI)
          read(iunit,*) (etav_resida(i),i=1,NVECAFSI)

!         WORKDA
          read(iunit,*) (eta_workda(i),i=1,3*NVECAFSI)
          read(iunit,*) (etav_workda(i),i=1,3*NVECAFSI)

!         IPNTARP
          read(iunit,*) (workdapos(i),i=1,14)

          close(iunit)
          write(6,'(A4)')  'Done'
        endif   ! ierr.eq.0 
      endif     ! nid.eq.0
      call err_chk(ierr,'Error reading fsi_arn_rst file.')

!     WORKDA
!     Probably need to ensure that workda is initialized to 0
      if (nid.eq.0) then
        do i=1,3
          do j=1,NVECAFSI
            k=workdapos(i) + j - 1    
            workda(k)=eta_workda(i+j-1)
!           prabal
!            write(6,*) 'workdapos,workda',k,workdapos(i),workda(k)
            k=workdapos(i) + NVECAFSI + j-1
            workda(k)=etav_workda(i+j-1)
!           prabal
!            write(6,*) 'workdapos,workda',k,workdapos(i),workda(k)
          enddo
        enddo
      endif

!     RESIDA
      vecptr0=ndim*NVECAV
      if (IFHEAT) vecptr0=vecptr0+NVECAT

      if (nid.eq.0) then
        j=vecptr0+1 
        call copy(RESIDA(j),eta_resida,NVECAFSI)     ! residual eta
        j=j+NVECAFSI
        call copy(RESIDA(j),etav_resida,NVECAFSI)    ! residual etav
        j=j+NVECAFSI
      else
        j=vecptr0+1 
        call rzero(RESIDA(j),NVECAFSI)     ! residual eta
        j=j+NVECAFSI
        call rzero(RESIDA(j),NVECAFSI)    ! residual etav
        j=j+NVECAFSI
      endif

!     ETA,ETAV        
!      call bcast(fsi_rst_eta, NVECAFSI*WDSIZE)
!      call bcast(fsi_rst_etav,NVECAFSI*WDSIZE)

      eta       = fsi_rst_eta(1) 
      etav      = fsi_rst_etav(1)
      etaa      = 0.
      fsi_alpha = 0.
      etav_s    = etav
      eta_s     = eta + etav_s*DT

      if (ifusermv) call fsi_meshv_pert


      return 
      end subroutine fsi_arnoldi_rstread

!---------------------------------------------------------------------- 
      subroutine fsi_workda2flds

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'SOLN_DEF'
      include 'SOLN'

      integer vecptr
      integer aptr
      integer i,j
      integer nfsi

      nfsi = NVECAFSI            ! nonzero only for nid=0

!     velocity
      vecptr= 0
      aptr  = IPNTARP(1)

      call copy(VXP,WORKDA(aptr),NVECAV)
      vecptr=vecptr+NVECAV
      call copy(VYP,WORKDA(aptr+vecptr),NVECAV)
      vecptr=vecptr+NVECAV
      if (IF3D) then
        call copy(VZP,WORKDA(aptr+vecptr),NVECAV)
        vecptr=vecptr+NVECAV
      endif

!     Pressure 
      if (TSTIFPR) then
        call copy(PRP,WORKDA(aptr+vecptr),NVECAP)
        vecptr=vecptr+NVECAP
      endif  

!     temperature
      if (IFHEAT) then
        call copy(TP,WORKDA(aptr+vecptr),NVECAT)
        vecptr=vecptr+NVECAT
      endif  

!     RESIDA
      j=vecptr+1 
      call copy(eta_resida,RESIDA(j),nfsi)     ! residual eta
      j=j+nfsi
      call copy(etav_resida,RESIDA(j),nfsi)    ! residual etav
      j=j+nfsi

!     WORKDA      
      do i=1,3
        j=vecptr
        call copy(eta_workda(i),WORKDA(IPNTARP(i)+j),nfsi)
        j=j+nfsi
        call copy(etav_workda(i),WORKDA(IPNTARP(i)+j),nfsi)
        j=j+nfsi
      enddo

!     Get how many degrees of freedom      
      call bcast(NVECAFSI, 1*ISIZE)
      nfsi = NVECAFSI

!     nid=0 is the relevant data point.
      call bcast(eta_resida,  nfsi*WDSIZE) 
      call bcast(etav_resida, nfsi*WDSIZE)

!     nid=0 is the relevant data point.      
      call bcast(eta_workda,  3*nfsi*WDSIZE) 
      call bcast(etav_workda, 3*nfsi*WDSIZE)

!     Now all processors have the structural values      
      eta  = eta_workda(1)
      etav = etav_workda(1)

!     Reset other nid values to zero 
      if (nid.ne.0) then
        call rzero(eta_resida,nfsi)
        call rzero(etav_resida,nfsi)
        call rzero(eta_workda,3*nfsi)
        call rzero(etav_workda,3*nfsi)
        NVECAFSI = 0
      endif  


      return
      end subroutine fsi_workda2flds
!---------------------------------------------------------------------- 
      subroutine fsi_flds2workda

!     Copy the vector resulting from
!     A*x to the arnoldi work array            

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'INPUT_DEF'
      include 'INPUT'           ! IFHEAT, IF3D
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, TP, V?MASK
      include 'MASS_DEF'
      include 'MASS'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !
      include 'FSI'                     ! eta,etav

      integer vecptr
      integer aptr

      vecptr=0
      aptr  = IPNTARP(2)      ! pointer for Arpack work array

      call col3(WORKDA(aptr+vecptr),VXP,FSI_MASK,NVECAV)
      vecptr=vecptr+NVECAV
      call col3(WORKDA(aptr+vecptr),VYP,FSI_MASK,NVECAV)
      vecptr=vecptr+NVECAV
      if (IF3D) then
        call col3(WORKDA(aptr+vecptr),VZP,FSI_MASK,NVECAV)
        vecptr=vecptr+NVECAV
      endif  

!     if we include pressure in the arnoldi vector
!     No masks here
      if (TSTIFPR) then
        call copy(WORKDA(aptr+vecptr),PRP,NVECAP)
        vecptr=vecptr+NVECAP
      endif  

!     temperature
      if (IFHEAT) then
!       Do we need a mask? Is it the same mask?            
        call col3(WORKDA(aptr+vecptr),TP,FSI_MASK,NVECAT)
        vecptr=vecptr+NVECAT
      endif

!     FSI
      if (IFFSI) then
        call copy(WORKDA(aptr+vecptr),eta,NVECAFSI)
        vecptr=vecptr+NVECAFSI
        call copy(WORKDA(aptr+vecptr),etav,NVECAFSI)
        vecptr=vecptr+NVECAFSI
      endif 

      return
      end subroutine fsi_flds2workda
!----------------------------------------------------------------------       
      subroutine fsi_flds2resida

!     Copy the vector to the arnoldi resid array
!     Only done once during initialization
!     (Not entirely sure why its needed) 

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'INPUT_DEF'
      include 'INPUT'           ! IFHEAT, IF3D
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, TP, V?MASK
      include 'MASS_DEF'
      include 'MASS'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !
      include 'FSI'                     ! eta,etav

      integer vecptr
      integer aptr

      integer ntott,ntotp

      real fsi_glsc2_wt           ! function
      real beta


!     Apply masks      
      call opcol2(vxp,vyp,vzp,fsi_mask,fsi_mask,fsi_mask)
      if (ifheat) then
!       Not sure if we need a mask here. Needs to be tested.
        ntott = nx1*ny1*nz1*nelt
        call col2(tp,fsi_mask,ntott)
      endif  

!     Normalize starting vector      
      beta  = fsi_glsc2_wt()
      beta  = sqrt(beta)
      call opcmult(vxp,vyp,vzp,1./beta)

      if (TSTIFPR) then
        ntotp=nx2*ny2*nz2*nelv
        call cmult(prp,1./beta,ntotp)
      endif  

      vecptr= 1

      call copy(RESIDA(vecptr),VXP,NVECAV)
      vecptr=vecptr+NVECAV
      call copy(RESIDA(vecptr),VYP,NVECAV)
      vecptr=vecptr+NVECAV
      if (IF3D) then
        call copy(RESIDA(vecptr),VZP,NVECAV)
        vecptr=vecptr+NVECAV
      endif  

!     if we include pressure in the arnoldi vector
      if (TSTIFPR) then
        call copy(RESIDA(vecptr),PRP,NVECAP)
        vecptr=vecptr+NVECAP
      endif

!     temperature
      if (IFHEAT) then
        call copy(RESIDA(vecptr),TP,NVECAT)
        vecptr=vecptr+NVECAT
      endif

!     FSI
      if (IFFSI) then
!       NVECAFSI=0 for nid.ne.0
        eta = eta/beta 
        call copy(RESIDA(vecptr),eta,NVECAFSI)
        vecptr=vecptr+NVECAFSI
        etav = etav/beta
        call copy(RESIDA(vecptr),etav,NVECAFSI)
        vecptr=vecptr+NVECAFSI
      endif

      return
      end subroutine fsi_flds2resida
!----------------------------------------------------------------------       
      subroutine fsi_outpostegv(iegv)
  
      implicit none
  
      include 'SIZE_DEF'
      include 'SIZE'
      include 'ARNOLDI_ARPACKD'
      include 'FSI'

      integer iegv              ! eigenvector no
      integer ierr

      character*132 fname
      character*1   fname1(132)
      equivalence (fname1,fname) 

      integer iunit
      save iunit 
        
      if (iegv.eq.1) then
!       Create filename
        call blank(fname1,132)
        write(fname,'(A15)') fsi_egv_arnoldi

!       only nid.eq.0 opens the file
        if (nid.eq.0) then
          call IO_freeid(iunit,ierr)     
          open(iunit,file=fname,status='unknown',action='write',
     $                  iostat=ierr)
          if (ierr.eq.0) then
            write(6,'(A24,1x,A17)') 'FSI ARN_EGV: Saving file',
     $               fname
          endif

!         write header      
          write(iunit,'(A5,2(1x,A25))') 'i','ETA', 'ETAV'
          write(iunit,'(I5,2(1x,E25.16E3))') iegv,eta,etav

          flush(iunit)
        endif   ! nid.eq.0

      else  

        if (nid.eq.0) then
!         Write structural component of eigenvector
          write(iunit,'(I5,2(1x,E25.16E3))') iegv,eta,etav
          flush(iunit)
        endif

      endif

      if (iegv.eq.ARNEGV) close(iunit)  ! done writing eigenvectors


      return
      end subroutine fsi_outpostegv
!----------------------------------------------------------------------       


!----------------------------------------------------------------------

