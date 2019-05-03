!       Author: Prabal
!       Description: Testing Orthogonalization of Arnoldi method
!       
!====================================================================== 
!---------------------------------------------------------------------- 
!     read parameters fluid-structure interaction 
      subroutine nek_arnoldi_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'NEK_ARNOLDI'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /NEKARN/ ifnekarn,northo,nek_uzawa,nekarn_ifpr,
     $                  ngs,pstep,pinistep 

!     default values
      ifnekarn          = .FALSE.        ! if perform my Arnoldi
      nek_uzawa         = .TRUE.         ! uzawa at the first time step?
      nekarn_ifpr       = .TRUE.         ! if include pressure in
                                         ! arnoldi vector
      northo            = arnkryl         ! no of vectors to save
      ngs               = 1              ! no of Gram-Schmidt passes
      pstep             = 100            ! No of steps between
                                         ! successive orthogonalization
      pinistep          = pstep          ! No of steps before taking the
                                         ! first krylov vector
                                         
!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=NEKARN,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading NEKARN parameters.$')

!     broadcast data
      call bcast(ifnekarn    , LSIZE)
      call bcast(nek_uzawa     , LSIZE)
      call bcast(nekarn_ifpr , LSIZE)
      call bcast(northo        , ISIZE)
      call bcast(ngs           , ISIZE)
      call bcast(pstep         , ISIZE)
      call bcast(pinistep      , ISIZE)

      return
      end subroutine nek_arnoldi_param_in
!----------------------------------------------------------------------
!     write parameters fluid-structure interaction
      subroutine nek_arnoldi_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'NEK_ARNOLDI'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /NEKARN/ ifnekarn,northo,nek_uzawa,nekarn_ifpr,
     $                  ngs,pstep,pinistep 

!     read the file
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=NEKARN,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing NEKARN parameters.$')

      return
      end subroutine nek_arnoldi_param_out
c-----------------------------------------------------------------------

      subroutine nek_arnoldi_main
      
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'             ! ifuzawa
      include 'NLFSI'
      include 'NEK_ARNOLDI'

      include 'MASS_DEF'
      include 'MASS'

      integer icalld
      data icalld /0/
      save icalld

      integer i,j,igs
      real hj,beta

      real nek_innerprod           ! function
      real nek_wt_innerprod

!     Do nothing if its not on
      if (.not.ifnekarn) return

      IFUZAWA = .false.
      if (northo.gt.arnkryl) then
        call exitti('northo > arnkryl, $', northo)
      endif

      if ((nkryl.eq.0).and.(istep.lt.pinistep)) return

      if (mod(istep,pstep).ne.0) then
        return
      endif

      if (icalld.eq.0) then
        call nek_arnoldi_init
        icalld=icalld+1
        return
      endif

      nsteps=nsteps+1
      lastep = 0

!     fill up Ax,Axtmp vector
      call getAx

!     Remove orthogonal projections
      i=nkryl
      do igs=1,ngs            ! No of GS passes
        do j=1,nkryl
          hj = nek_wt_innerprod(Ax,Qortho(1,j),vlen) ! find projection
          hessen(j,i)=hessen(j,i) + hj
          call add2s2(Axtmp,Qortho(1,j),-hj,vlen)    ! remove the projection
        enddo
        call copy(Ax,Axtmp,vlen)
      enddo
!     Add residual to hessenberg Matrix       
      beta = nek_wt_innerprod(Ax,Ax,vlen)
      beta = sqrt(beta)
      hessen(i+1,i)=beta
      
      call outhessen

!     Add unit vector to Qortho
      call cmult(Ax,1./beta,vlen)
      call copy(Qortho(1,i+1),Ax,vlen)
      nkryl=nkryl+1
      if (nio.eq.0) write(6,'(I5,1x,A9,E25.16E3)') nkryl,
     $                            'Residual=',beta
 
!      call qorthocheck

!     Reset time-stepping
      call nek_stepper_restart       

      if (nkryl.eq.northo+1) then
        call nek_arnoldi_finalize
      endif

      return
      end subroutine nek_arnoldi_main
!----------------------------------------------------------------------

      subroutine nek_arnoldi_init

!     Initialize my Arnoldi

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEK_ARNOLDI'

      integer i
      real vnorm
      real nek_wt_innerprod
      real glsc3

      integer ntot2

      ntot2=nx2*ny2*nz2*nelv

!     Set weights. 
!     Also sets Arnoldi vector length (vlen)
      call set_arnoldi_weight
!      call set_arnoldi_wtone

!     Create Masks
      call set_arnoldi_msk

!     Zero Hessenberg matrix
      i=arnkryl
      call rzero(hessen,(i+1)*i)

!     Zero Qortho matrix
      i=qlen2*arnkryl
      call rzero(qortho,i)

!     get starting vector        
      call getAx

!     Normalize starting vector
      vnorm = nek_wt_innerprod(Ax,Ax,vlen)
      vnorm = sqrt(vnorm)
      call cmult(Ax,1./vnorm,vlen)

      if (nio.eq.0) write(6,*) 'Initial Residual=',vnorm

!     Add starting vector to Krylov space
      call copy(Qortho(1,1),Ax,vlen)
      nkryl = 0       ! we don't want to pick the random initial condition
      pinistep = pstep  ! we don't want the long initial

!      call qorthocheck

!     Restart stepper
      call nek_stepper_restart


      return
      end subroutine nek_arnoldi_init

!----------------------------------------------------------------------

      subroutine set_arnoldi_msk

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'
      include 'NLFSI'
      include 'NEK_ARNOLDI'  
        
      integer ifld,nfaces,iel,iface
      character cb*3
      integer i,ntot,ntotp

      real mskfld(lx1,ly1,lz1,lelv)

!     Build a mask to remove only 'v  '
!     and preserving the 'mv ' points.
      call rone(mskfld,lx1*ly1*lz1*lelv)
      ifld = 1
      nfaces=2*ndim
      do iel=1,nelv
        do iface = 1,nfaces
          cb = cbc(iface,iel,ifld)
          if (cb.eq.'v  ') then
!           Put zeros on this face
            call facev(mskfld,iel,iface,0.,nx1,ny1,nz1)
          endif
        enddo ! iel
      enddo   ! iface

      ntot=nx1*ny1*nz1*nelv
      i=1
      call copy(Axmsk(i),mskfld,ntot)
      i=i+ntot
      call copy(Axmsk(i),mskfld,ntot)
      i=i+ntot
      if (if3d) then
        call copy(Axmsk(i),mskfld,ntot)
        i=i+ntot
      endif

      if (nekarn_ifpr) then
        ntotp=nx2*ny2*nz2*nelv
        call rone(Axmsk(i),ntotp)
        i=i+ntotp
      endif  

      if (IFFSI) then
        if (nid.eq.0) then       
          call copy(Axmsk(i),1.0,1)
          i=i+1
          call copy(Axmsk(i),1.0,1)
          i=i+1
        endif
      elseif (IFNLFSI) then
        if (nid.eq.0) then
          call copy(Axmsk(i),1.0,1)
          i=i+1
          call copy(Axmsk(i),1.0,1)
          i=i+1
        endif      
      endif

      return
      end subroutine set_arnoldi_msk
!----------------------------------------------------------------------

      subroutine getAx

!     Put action of the matrix into the vector Ax

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'
      include 'NLFSI'
      include 'NEK_ARNOLDI'

      integer i,ntot,ntotp
      
      ntot=nx1*ny1*nz1*nelv
      i=1
      call col3(Ax(i),vxp,Axmsk(i),ntot)
      i=i+ntot
      call col3(Ax(i),vyp,Axmsk(i),ntot)
      i=i+ntot
      if (if3d) then
        call col3(Ax(i),vzp,Axmsk(i),ntot)
        i=i+ntot
      endif

      if (nekarn_ifpr) then
        ntotp=nx2*ny2*nz2*nelv
        call copy(Ax(i),prp,ntotp)
        i=i+ntotp
      endif  

      if (IFFSI) then
        if (nid.eq.0) then       
          call col3(Ax(i),eta,Axmsk(i),1)
          i=i+1
          call col3(Ax(i),etav,Axmsk(i),1)
          i=i+1
        endif
      elseif (IFNLFSI) then
        if (nid.eq.0) then
          call col3(Ax(i),psi,Axmsk(i),1)
          i=i+1
          call col3(Ax(i),psiv,Axmsk(i),1)
          i=i+1
        endif      
      endif

      call copy(Axtmp,Ax,i-1)

      return
      end subroutine getAx
!---------------------------------------------------------------------- 

      subroutine setflds

!     assign velocity/fsi fields from Ax entries

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'ADJOINT_DEF'
      include 'ADJOINT'
      include 'FSI'
      include 'NLFSI'
      include 'NEK_ARNOLDI'

      integer i,ntot,ntotp
      
      ntot=nx1*ny1*nz1*nelv
      i=1
      call col3(vxp,Ax(i),Axmsk(i),ntot)
      i=i+ntot
      call col3(vyp,Ax(i),Axmsk(i),ntot)
      i=i+ntot
      if (if3d) then
        call col3(vzp,Ax(i),Axmsk(i),ntot)
        i=i+ntot
      endif

      ntotp=nx2*ny2*nz2*nelv
      if (nekarn_ifpr) then
        call copy(prp,Ax(i),ntotp)
        i=i+ntotp
      else
        call rzero(prp,ntotp)  
      endif  

      if (IFFSI) then
        if (nid.eq.0) then       
          call copy(eta,Ax(i),1)
          i=i+1
          call copy(etav,Ax(i),1)
          i=i+1
        endif

        call bcast(eta,wdsize)
        call bcast(etav,wdsize)

        eta_s  = eta
        etav_s = etav

      elseif (IFNLFSI) then
        if (nid.eq.0) then
          call copy(psi,Ax(i),1)
          i=i+1
          call copy(psiv,Ax(i),1)
          i=i+1
        endif
        call bcast(psi,wdsize)
        call bcast(psiv,wdsize)

        psi_s  = psi
        psiv_s = psiv

      endif


      return
      end subroutine setflds
!---------------------------------------------------------------------- 

      real function nek_innerprod(x,y,n)

      implicit none

      integer n
      real x(1),y(1)
      real glsc2

      nek_innerprod=glsc2(x,y,n)

      end function nek_innerprod
!---------------------------------------------------------------------- 

      real function nek_wt_innerprod(x,y,n)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEK_ARNOLDI'

      integer n
      real x(1),y(1)
      real glsc3

      nek_wt_innerprod=glsc3(x,y,ArWt,n)


      end function nek_wt_innerprod
!---------------------------------------------------------------------- 

      subroutine set_arnoldi_weight

!     Set weights for the Arnoldi vector

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'MASS_DEF'
      include 'MASS'        ! BM1
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'
      include 'NLFSI'
      include 'NEK_ARNOLDI'

      integer i,ntot,ntotp
      
      ntot=nx1*ny1*nz1*nelv
      i=1
      call copy(ArWt(i),BM1,ntot)
      i=i+ntot
      call copy(ArWt(i),BM1,ntot)
      i=i+ntot
      if (if3d) then
        call copy(ArWt(i),BM1,ntot)
        i=i+ntot
      endif

      
      if (nekarn_ifpr) then
        ntotp=nx2*ny2*nz2*nelv
!        call copy(ArWt(i),BM2,ntotp)

!       If we use the semi-norm then
!       weight of pressure is zero        
        call rzero(ArWt(i),ntotp)
        i=i+ntotp
      endif  

!      if (IFFSI) then
!        if (nid.eq.0) then       
!          call copy(ArWt(i),fsi_eta_wt,1)
!          i=i+1
!          call copy(ArWt(i),fsi_etav_wt,1)
!          i=i+1
!        endif
!      elseif (IFNLFSI) then
!        if (nid.eq.0) then
!          call copy(ArWt(i),nlfsi_psi_wt,1)
!          i=i+1
!          call copy(ArWt(i),nlfsi_psiv_wt,1)
!          i=i+1
!        endif 
!      endif

!     Size of krylov vector
      vlen = i-1
!      if (nio.eq.0) write(6,*) 'Arnoldi vector length=',vlen

      return
      end subroutine set_arnoldi_weight
!---------------------------------------------------------------------- 

      subroutine set_arnoldi_wtone

!     Set weights for the Arnoldi vector to identity

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'FSI'
      include 'NLFSI'
      include 'NEK_ARNOLDI'

      integer i,ntot,ntotp
      
      ntot=nx1*ny1*nz1*nelv
      i=1
      call rone(ArWt(i),ntot)
      i=i+ntot
      call rone(ArWt(i),ntot)
      i=i+ntot
      if (if3d) then
        call rone(ArWt(i),ntot)
        i=i+ntot
      endif

      if (nekarn_ifpr) then
        ntotp=nx2*ny2*nz2*nelv
!        call rone(ArWt(i),ntotp)

!       If we use the semi-norm then
!       weight of pressure is zero        
        call rzero(ArWt(i),ntotp)
       
        i=i+ntotp
      endif  

      if (IFFSI) then
        if (nid.eq.0) then       
          call copy(ArWt(i),1.0,1)
          i=i+1
          call copy(ArWt(i),1.0,1)
          i=i+1
        endif
      elseif (IFNLFSI) then
        if (nid.eq.0) then
          call copy(ArWt(i),1.0,1)
          i=i+1
          call copy(ArWt(i),1.0,1)
          i=i+1
        endif 
      endif

      return
      end subroutine set_arnoldi_wtone
!---------------------------------------------------------------------- 

      subroutine outhessen

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEK_ARNOLDI'

      character outfmt*64
      integer i

      call blank(outfmt,64)
      write(outfmt,'(A10,I5,A14)') '(A6,1x,I5,',nkryl+1,'(E25.16E3,1x))'
      if (nio.eq.0) then
        write(6,outfmt) 'hessen',nkryl, (hessen(i,nkryl),i=1,nkryl+1)
      endif  

      return
      end subroutine outhessen
!----------------------------------------------------------------------

      subroutine qorthocheck

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEK_ARNOLDI'

      character outfmt*64
      integer i
      real r(arnkryl),rmax
      real nek_wt_innerprod
      real vlamax

      do i=1,nkryl
        r(i)=nek_wt_innerprod(Qortho(1,i),Qortho(1,nkryl),vlen)
      enddo

      call blank(outfmt,64)
      i=2
      write(outfmt,'(A10,I5,A13)') '(A6,1x,I5,',i,'(E17.8E3,1x))'
      if (nio.eq.0) then
        rmax=vlamax(r,nkryl-1)
        write(6,outfmt) 'qortho',nkryl,rmax,r(nkryl)
      endif  

      return
      end subroutine qorthocheck
!----------------------------------------------------------------------

      subroutine nek_arnoldi_finalize

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'NEK_ARNOLDI'

      integer i,ntot

      ntot = nx1*ny1*nz1*nelv

      call outfile_hessen
      if (nekarn_ifpr) then
        do i=1,nkryl
          call outpost(Qortho(1,i),Qortho(1+ntot,i),
     $          Qortho(1+2*ntot,i),Qortho(1+ndim*ntot,i),tp,'qor')
        enddo
      else
        do i=1,nkryl
          call outpost(Qortho(1,i),Qortho(1+ntot,i),
     $          Qortho(1+2*ntot,i),prp,tp,'qor')
        enddo
      endif

      call exitt 

      return
      end subroutine nek_arnoldi_finalize
!----------------------------------------------------------------------

      subroutine outfile_hessen

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEK_ARNOLDI'

      character outfmt*64
      integer i,j

      call blank(outfmt,64)
      write(outfmt,'(A1,I5,A14)') '(',northo,'(E25.16E3,1x))'

      if (nid.eq.0) then
        open(unit=10101,file='hessenberg',status='unknown',
     $      form='formatted')
        do i=1,northo+1
          write(10101,outfmt) (hessen(i,j),j=1,northo)
        enddo
        close(10101) 
      endif  

      return
      end subroutine outfile_hessen
!----------------------------------------------------------------------
      subroutine nek_stepper_restart

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'ADJOINT_DEF'
      include 'ADJOINT'
      include 'FSI'
      include 'NLFSI'
      include 'NEK_ARNOLDI'
      
      integer ntot
      integer ntot2

      ntot=nx1*ny1*nz1*nelv
      ntot2=nx2*ny2*nz2*nelv

      istep=0
      time=0.

!     set vxp,psi etc from Ax      
      call setflds

!      if (iffsi) then
!        if (ifadj) then
!          continue 
!        else
!          call fsi_meshv_pert
!        endif  
!      elseif (ifnlfsi) then
!        if (ifadj) then
!          call nlfsi_meshv_adj
!        else
!          call nlfsi_meshv_pert
!        endif  
!      endif

      ifuzawa = nek_uzawa

      return
      end subroutine nek_stepper_restart
!---------------------------------------------------------------------- 

      subroutine calc_prp
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
c
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'ADJOINT_DEF'
      include 'ADJOINT'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'         ! ifield
      include 'MASS_DEF'
      include 'MASS'

!     include 'CTIMER'

      real w1,w2,w3
      real dv1,dv2,dv3,dp
      real h1,h2,h2inv
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

      integer ntot1,ntot2
      integer intype
      integer ifld_old,jpold

      ntot1  = nx1*ny1*nz1*nelv
      ntot2  = nx2*ny2*nz2*nelv

      call rone(h2,ntot1)
      call rzero(h1,ntot1)
      call rone(h2inv,ntot1)

      ifld_old = ifield
      ifield = 1
      jpold = jp
      jp = 1
      call opzero(bfxp,bfyp,bfzp)
      if (ifadj) then
        call advabp_adjoint
      else
        call advabp
      endif  

!     Remove Mass matrix 
      call opcol2(bfxp,bfyp,bfzp,binvm1,binvm1,binvm1)            
      call opdiv(dp,bfxp,bfyp,bfzp)

!      call chsign(dp,ntot2)
      call ortho(dp)

!     How do I verify that this change is correct?      
      intype = 1
      istep=99999
      if (nid.eq.0) write(6,*) 'Recalculating pressure'
      call esolver (dp,h1,h2,h2inv,intype)
      istep=0
      jp=jpold
      ifield=ifld_old

      call copy(prp,dp,ntot2)

      return
      end subroutine calc_prp
c------------------------------------------------------------------------








