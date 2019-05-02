!======================================================================
!     Description: Routines for coupling with ARNOLDI
!     Author: Prabal S. Negi
!     Last Modified: 11-01-2019      
!======================================================================       
      real function nlfsi_glsc2_wt()

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
!      include 'CONHT'           ! CHCST_SC,CHCFF_V,CHCFF_T
      include 'NLFSI'             ! IFNLFSI

!     local variables
      integer ntotv, ntott
      real scal, tmp
      integer i
      real op_glsc2_wt
      real glsc3  

      integer nfsi

!     functions
      real glsum

      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT
      nfsi  = NVECA_NLFSI           ! nonzero only for nid=0

      scal = 0.
      scal = op_glsc2_wt(vxp,vyp,vzp,vxp,vyp,vzp,bm1)

!     Zero weight for pressure

!     for conjugate heat transfer
      if (IFHEAT) then
        scal = scal + glsc3(tp,tp,bm1,ntott)
      endif  

!     Fluid-Structure-Interaction
      tmp = glsc3(psi,psi,nlfsi_psi_wt,nfsi)
      tmp = tmp + glsc3(psiv,psiv,nlfsi_psiv_wt,nfsi)
      scal = scal + tmp

      nlfsi_glsc2_wt = scal 
      
      return
      end function nlfsi_glsc2_wt

!---------------------------------------------------------------------- 
      subroutine nlfsi_opcmult (v1,v2,v3,t1,s1,s2,const)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
      include 'NLFSI'             ! IFNLFSI

!     argument list
      real v1(1),v2(1),v3(1),p1(1),t1(1)
      real s1(1),s2(1)
      real const

!     local variables
      integer ntotv, ntotp, ntott
      integer nfsi

      ntotv = NX1*NY1*NZ1*NELV
      ntotp = NX2*NY2*NZ2*NELV
      ntott = NX1*NY1*NZ1*NELT
      nfsi  = NVECA_NLFSI           ! nonzero only for nid=0

      if (IFFLOW) then
        call cmult(v1,const,ntotv)
        call cmult(v2,const,ntotv)
        if (IF3D) call cmult(v3,const,ntotv)
      endif

      if (IFHEAT) call cmult(t1,const,ntott)

      if (IFNLFSI)  then
         call cmult(s1,const,nfsi)
         call cmult(s2,const,nfsi)
      endif  

      return
      end subroutine nlfsi_opcmult

!---------------------------------------------------------------------- 
      subroutine nlfsi_arnoldi_restart

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'NLFSI'
      include 'ADJOINT_DEF'
      include 'ADJOINT'

!     Zero lag arrays      
      call rzero(psilag,lorder-1)
      call rzero(psivlag,lorder-1)

!     set boundary conditions
      if (ifadj) then
!       velocity/position extrapolation
!       Non-linear iteration so it doesn't matter            
        psiv_s = psiv
        psi_s  = psi
        psia   = 0.
           
        call nlfsi_meshv_adj
      else
!       velocity/position extrapolation      
        psiv_s = psiv
        psi_s  = psi + DT*psiv_s
        psia   = 0.
           
        call nlfsi_meshv_pert
      endif

      call nlfsi_output

      return
      end subroutine nlfsi_arnoldi_restart
!----------------------------------------------------------------------       

      subroutine nlfsi_arnoldi_rstsave

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'NLFSI'
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
      integer nfsi

      ierr=0
      nfsi = NVECA_NLFSI

      nlfsi_rst_psi(1)  = psi
      nlfsi_rst_psiv(1) = psiv

      if (nid.eq.0) then
!       Create filename
        call blank(fname1,132)
        write(fname,'(A13,A1,i1)') nlfsi_rst_arnoldi,dot,ARNISTOP

!       Write output
        call IO_freeid(iunit,ierr)     
        open(iunit,file=fname,status='unknown',action='write',
     $    iostat=ierr)
        if (ierr.eq.0) then
          write(6,'(A26,1x,A13)') 'NLFSI ARN_RST: Saving file',
     $             fname

!         RESIDA 
          write(iunit,*) (psi_resida(i),i=1,nfsi)
          write(iunit,*) (psiv_resida(i),i=1,nfsi)

!         WORKDA
          write(iunit,*) (psi_workda(i),i=1,3*nfsi)
          write(iunit,*) (psiv_workda(i),i=1,3*nfsi)

!         IPNTARP
          write(iunit,*) (IPNTARP(i),i=1,14)

          close(iunit)
          write(6,'(A4)')  'Done'
        endif   ! ierr.eq.0 
      endif       ! nid.eq.0

      call err_chk(ierr,'Error opening nlfsi_arn_rst file.')
      
      return 
      end subroutine nlfsi_arnoldi_rstsave

!---------------------------------------------------------------------- 

      subroutine nlfsi_arnoldi_rstread

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'NLFSI'
      include 'CHKPOINT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'
      include 'ADJOINT_DEF'
      include 'ADJOINT'
      
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

      real bd_psiv
      integer ilag
      integer iunit,ierr

      integer workdapos(14)
      integer i,j,k

      integer vecptr0
      integer nfsi

      character*132 tmp


      ierr=0 
      nfsi = NVECA_NLFSI

      if (istep.gt.0) return

!     Create filename
      call blank(fname1,132)
      write(fname,'(A13,A1,i1)') nlfsi_rst_arnoldi,dot,ARNISTART

!     Read restart file 
      if (nid.eq.0) then
        call IO_freeid(iunit,ierr)
        open(iunit, file=fname, status='old',action='read',
     $    iostat=ierr)
        if (ierr.eq.0) then
          write(6,'(A27,1x,A13)') 'NLFSI Restart: Opening file',
     $             fname

!         RESIDA 
          read(iunit,*) (psi_resida(i),i=1,nfsi)
          read(iunit,*) (psiv_resida(i),i=1,nfsi)

!         WORKDA
          read(iunit,*) (psi_workda(i),i=1,3*nfsi)
          read(iunit,*) (psiv_workda(i),i=1,3*nfsi)

!         IPNTARP
          read(iunit,*) (workdapos(i),i=1,14)

          close(iunit)
          write(6,'(A4)')  'Done'
        endif   ! ierr.eq.0 
      endif     ! nid.eq.0
      call err_chk(ierr,'Error reading nlfsi_arn_rst file.')

!     WORKDA
!     Probably need to ensure that workda is initialized to 0
      if (nid.eq.0) then
        do i=1,3
          do j=1,nfsi
            k=workdapos(i) + j - 1    
            workda(k)=psi_workda(i+j-1)
!           prabal
!            write(6,*) 'workdapos,workda',k,workdapos(i),workda(k)
            k=workdapos(i) + nfsi + j-1
            workda(k)=psiv_workda(i+j-1)
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
        call copy(RESIDA(j),psi_resida,nfsi)     ! residual psi
        j=j+nfsi
        call copy(RESIDA(j),psiv_resida,nfsi)    ! residual psiv
        j=j+nfsi
      else
        j=vecptr0+1 
        call rzero(RESIDA(j),nfsi)     ! residual psi
        j=j+nfsi
        call rzero(RESIDA(j),nfsi)    ! residual psiv
        j=j+nfsi
      endif

      psi       = nlfsi_rst_psi(1) 
      psiv      = nlfsi_rst_psiv(1)
      psia      = 0.
      psiv_s    = psiv
      psi_s     = psi

!     Set boundary conditions
      if (ifadj) then          
        call nlfsi_meshv_adj
      else
        call nlfsi_meshv_pert
      endif

      return 
      end subroutine nlfsi_arnoldi_rstread

!----------------------------------------------------------------------

      subroutine nlfsi_weight_fcn

!     Called when IDOARP.eq.2
!     if y=Ax is the Arnoldi vector
!     This subroutine calculates w = My
!     and stores it in WORDA(IPNTARP(2))
!     Needed by Arnoldi to perform M-weighted
!     orthogonalization

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'NLFSI'
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
      integer nfsi
      integer vecptr,aptr


      ntotv = lx1*ly1*lz1*nelv
      ntott = lx1*ly1*lz1*nelt
      nfsi = NVECA_NLFSI            ! nonzero only for nid=0

      call copy(bmv_wt,BM1,ntotv)
      call copy(bmt_wt,BM1,ntott)

!      call rone(bmv_wt,ntotv)
!      call rone(bmt_wt,ntott)    

!     velocity
      vecptr=0
      aptr  = IPNTARP(2)

      call copy(WORKDA(aptr+vecptr),bmv_wt,NVECAV)
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

!     Temperature. 
!     Test before use
      if(IFHEAT) then
        call copy(WORKDA(aptr+vecptr),
     $       bmt_wt,NVECAT)
             vecptr=vecptr+NVECAT
      endif

!     nfsi=0 for nid.ne.0      
      call copy(WORKDA(aptr+vecptr),nlfsi_psi_wt,nfsi)
      vecptr=vecptr+nfsi
      call copy(WORKDA(aptr+vecptr),nlfsi_psiv_wt,nfsi)
      vecptr=vecptr+nfsi

      call col2(WORKDA(IPNTARP(2)),WORKDA(IPNTARP(1)),NVECAS)

      return
      end subroutine nlfsi_weight_fcn

!----------------------------------------------------------------------             

      subroutine nlfsi_workda2flds

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'NLFSI'
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

      real glsc2,tmp


      nfsi = NVECA_NLFSI            ! nonzero only for nid=0

!     velocity
      vecptr=0
      aptr=IPNTARP(1)

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
      call copy(psi_resida,RESIDA(j),nfsi)     ! residual psi
      j=j+nfsi
      call copy(psiv_resida,RESIDA(j),nfsi)    ! residual psiv
      j=j+nfsi

!     WORKDA 
      do i=1,3
        j=vecptr
        call copy(psi_workda(i),WORKDA(IPNTARP(i)+j),nfsi)
        j=j+nfsi
        call copy(psiv_workda(i),WORKDA(IPNTARP(i)+j),nfsi)
        j=j+nfsi
      enddo

!     Get how many degrees of freedom      
      call bcast(NVECA_NLFSI, 1*ISIZE)
      nfsi = NVECA_NLFSI

!     nid=0 is the relevant data point.
      call bcast(psi_resida,  nfsi*WDSIZE) 
      call bcast(psiv_resida, nfsi*WDSIZE)

!     nid=0 is the relevant data point.      
      call bcast(psi_workda,  3*nfsi*WDSIZE) 
      call bcast(psiv_workda, 3*nfsi*WDSIZE)

!     Now all processors have the structural values      
      psi  = psi_workda(1)
      psiv = psiv_workda(1)

!     Reset other nid values to zero 
      if (nid.ne.0) then
        call rzero(psi_resida,nfsi)
        call rzero(psiv_resida,nfsi)
        call rzero(psi_workda,3*nfsi)
        call rzero(psiv_workda,3*nfsi)
        NVECA_NLFSI = 0
      endif  

!     Boundary conditions are set later in the call to
!     nlfsi_arnoldi_restart 

      return
      end subroutine nlfsi_workda2flds
!---------------------------------------------------------------------- 
      subroutine nlfsi_flds2workda

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
      include 'NLFSI'                   ! psi,psiv

      integer vecptr
      integer aptr,aptr2

      real bmv_wt(lx1,ly1,lz1,lelv)
      real bmt_wt(lx1,ly1,lz1,lelt)

      integer ntotv,ntott

      vecptr= 0
      aptr  = IPNTARP(2)      ! pointer for Arpack work array

      call col3(WORKDA(aptr),VXP,NLFSI_MASK,NVECAV)
      vecptr=vecptr+NVECAV
      call col3(WORKDA(aptr+vecptr),VYP,NLFSI_MASK,NVECAV)
      vecptr=vecptr+NVECAV
      if (IF3D) then
        call col3(WORKDA(aptr+vecptr),VZP,NLFSI_MASK,NVECAV)
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
!       Not sure if we need a mask here. Needs to be tested.
        call col3(WORKDA(aptr+vecptr),TP,NLFSI_MASK,NVECAT)
        vecptr=vecptr+NVECAT
      endif

!     FSI
      if (IFNLFSI) then
!       NVECA_NLFSI=0 for nid.ne.0            
        call copy(WORKDA(aptr+vecptr),psi,NVECA_NLFSI)
        vecptr=vecptr+NVECA_NLFSI
        call copy(WORKDA(aptr+vecptr),psiv,NVECA_NLFSI)
        vecptr=vecptr+NVECA_NLFSI
      endif

      return
      end subroutine nlfsi_flds2workda
!----------------------------------------------------------------------       
      subroutine nlfsi_flds2resida

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
      include 'NLFSI'                   ! psi,psiv

      integer vecptr
      integer aptr,aptr2

      integer ntotv,ntott,ntotp

      real nlfsi_glsc2_wt           ! function
      real beta


!     Apply masks      
      call opcol2(vxp,vyp,vzp,nlfsi_mask,nlfsi_mask,nlfsi_mask)
      if (ifheat) then
!       Not sure if we need a mask here. Needs to be tested.
        ntott = nx1*ny1*nz1*nelt
        call col2(tp,nlfsi_mask,ntott)
      endif  
     
!     Normalize starting vector
      beta = nlfsi_glsc2_wt()
      beta = sqrt(beta)
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
      if (IFNLFSI) then
!       NVECA_NLFSI=0 for nid.ne.0
        psi = psi/beta 
        call copy(RESIDA(vecptr),psi,NVECA_NLFSI)
        vecptr=vecptr+NVECA_NLFSI
        psiv = psiv/beta
        call copy(RESIDA(vecptr),psiv,NVECA_NLFSI)
        vecptr=vecptr+NVECA_NLFSI
      endif

      return
      end subroutine nlfsi_flds2resida
!----------------------------------------------------------------------       
      subroutine nlfsi_outpostegv(iegv)
  
      implicit none
  
      include 'SIZE_DEF'
      include 'SIZE'
      include 'ARNOLDI_ARPACKD'
      include 'NLFSI'

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
        write(fname,'(A17)') nlfsi_egv_arnoldi

!       only nid.eq.0 opens the file
        if (nid.eq.0) then
          call IO_freeid(iunit,ierr)     
          open(iunit,file=fname,status='unknown',action='write',
     $                  iostat=ierr)
          if (ierr.eq.0) then
            write(6,'(A26,1x,A17)') 'NLFSI ARN_EGV: Saving file',
     $               fname
          endif

!         write header     
          write(iunit,'(A5,2(1x,A25))') 'i','PSI', 'PSIV'
          write(iunit,'(I5,2(1x,E25.16E3))') iegv,psi,psiv
          flush(iunit)
        endif   ! nid.eq.0

      else  

        if (nid.eq.0) then
!         Write structural component of eigenvector
          write(iunit,'(I5,2(1x,E25.16E3))') iegv,psi,psiv
          flush(iunit)
        endif

      endif

      if (iegv.eq.ARNEGV) close(iunit)  ! done writing eigenvectors
  

      return
      end subroutine nlfsi_outpostegv
!----------------------------------------------------------------------       

      subroutine nlfsi_apply_mask

!     Copy the vector resulting from
!     A*x to the arnoldi work array            

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'INPUT_DEF'
      include 'INPUT'           ! IFHEAT, IF3D
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, TP, V?MASK
      include 'NLFSI'                   ! psi,psiv

      integer ntotv,ntott

!     Apply masks      
      call opcol2(vxp,vyp,vzp,nlfsi_mask,nlfsi_mask,nlfsi_mask)
      if (ifheat) then
!       Not sure if we need a mask here. Needs to be tested.
        ntott = nx1*ny1*nz1*nelt
        call col2(tp,nlfsi_mask,ntott)
      endif  


      return
      end subroutine nlfsi_apply_mask
!----------------------------------------------------------------------       




