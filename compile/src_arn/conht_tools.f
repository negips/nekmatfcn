!=======================================================================
!     Adam Peplinski; 2015.10.09
!     Tools for conjugated heat transfer to treat velocity and 
!     temperature together
!     Parameters used by this set of subroutines:
!     CONHT:
!     CHCST_SC, CHCFF_V, CHCFF_T - velocity and temperature scaling
!                                  factors
!     CHGRX, CHGRY, CHGRZ - gravitational acceleration
!     CHRA, CHPR - Rayleight and Prandtl number
!=======================================================================
!***********************************************************************
!     read parameters CONHT
      subroutine cht_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'CONHT'           ! CHCST_SC,CHCFF_V,CHCFF_T,CHGR[XYZ]

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr, len
      real rtmp(6)

!     namelists
      namelist /CONHT/ CHCST_SC,CHCFF_V,CHCFF_T,CHGRX,CHGRY,CHGRZ
!-----------------------------------------------------------------------
!     default values
      CHCST_SC = 3.36558d0
      CHCFF_V = 0.5
      CHCFF_T = 0.5
      CHGRX = 0.0
      CHGRY = 1.0
      CHGRZ = 0.0
!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=CONHT,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading CONHT parameters.$')

!     broadcast data
      if (NID.eq.0) then
         rtmp(1) = CHCST_SC
         rtmp(2) = CHCFF_V
         rtmp(3) = CHCFF_T
         rtmp(4) = CHGRX
         rtmp(5) = CHGRY
         rtmp(6) = CHGRZ
      endif
      len = 6*WDSIZE
      call bcast(rtmp,len)
      if (NID.ne.0) then
         CHCST_SC = rtmp(1)
         CHCFF_V = rtmp(2)
         CHCFF_T = rtmp(3)
         CHGRX = rtmp(4)
         CHGRY = rtmp(5)
         CHGRZ = rtmp(6)
      endif

      return
      end
!***********************************************************************
!     write parameters CONHT
      subroutine cht_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'CONHT'           ! CHCST_SC,CHCFF_V,CHCFF_T,CHGR[XYZ]

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /CONHT/ CHCST_SC,CHCFF_V,CHCFF_T,CHGRX,CHGRY,CHGRZ
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=CONHT,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing CONHT parameters.$')

      return
      end
!***********************************************************************
!     init parameters and variables in CONHT and ADJOINT
      subroutine cht_init()
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'INPUT_DEF'
      include 'INPUT'           ! PARAM
      include 'SOLN_DEF'
      include 'SOLN'            ! T
      include 'ADJOINT_DEF'
      include 'ADJOINT'         ! G_ADJ, BETA_B, DTD[XYZ]
      include 'CONHT'           ! CHCST_SC,CHCFF_V,CHCFF_T,CHGR[XYZ]

!     local variables
      integer ierr
!-----------------------------------------------------------------------
!     Rayleight and Prandtl numbers
      CHRA = abs(PARAM(2))
      CHPR = abs(PARAM(1))
      BETA_B = CHPR

!     gravity
      G_ADJ(1) = CHGRX
      G_ADJ(2) = CHGRY
      G_ADJ(3) = CHGRZ

!     temperature gradient
      call gradm1(DTDX,DTDY,DTDZ,T)

      return
      end
!***********************************************************************
!     calcualte forcing ralted to conjugated heat transfer
      subroutine cht_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFHEAT, CPFLD
      include 'PARALLEL_DEF'
      include 'PARALLEL'        ! GLLEL
      include 'TSTEP_DEF'
      include 'TSTEP'           ! IFIELD
      include 'SOLN_DEF'
      include 'SOLN'            ! JP, T, TP
      include 'ADJOINT_DEF'
      include 'ADJOINT'         ! IFADJ, G_ADJ, DTD[XYZ]
      include 'CONHT'           ! CHGR[XYZ]

!     argument list
      real ffx, ffy, ffz
      integer ix,iy,iz,ieg

!     local variables
      integer e, ip
      integer rtmp
!-----------------------------------------------------------------------
      if (IFHEAT) then
         e=GLLEL(ieg)
         if (JP.eq.0) then
            rtmp = T(ix,iy,iz,e,IFIELD)/CPFLD(1,2)
            ffx = ffx + CHGRX*rtmp
            ffy = ffy + CHGRY*rtmp
            if (IF3D) ffz = ffz + CHGRZ*rtmp
         else
            ip=ix+NX1*(iy-1+NY1*(iz-1+NZ1*(e-1)))  
            if (.not.IFADJ) then
               rtmp = TP(ip,IFIELD,JP)/CPFLD(1,2)
               ffx = ffx + G_ADJ(1)*rtmp
               ffy = ffy + G_ADJ(2)*rtmp
               if (IF3D) ffz = ffz + G_ADJ(3)*rtmp
            else
               ffx = ffx - DTDX(ip)*TP(ip,IFIELD,JP)
               ffy = ffy - DTDY(ip)*TP(ip,IFIELD,JP)
               if (IF3D) ffz = ffz - DTDZ(ip)*TP(ip,IFIELD,JP)
            end if
         end if
      endif

      return
      end
!***********************************************************************
!     set cpfld
      subroutine cht_cpfld_set()
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'INPUT_DEF'
      include 'INPUT'           ! CPFLD, PARAM
      include 'ADJOINT_DEF'
      include 'ADJOINT'         ! IFADJ
      include 'CONHT'           ! CHRA, CHPR

!     local variables
      integer ierr
!-----------------------------------------------------------------------
      if (IFHEAT) then
         if (IFADJ) then
            CPFLD(1,1)=CHPR/sqrt(CHRA)
            CPFLD(1,2)=1.0

            CPFLD(2,1)=1.0/sqrt(CHRA)
            CPFLD(2,2)=1.0
         else
            CPFLD(1,1)=1.0/sqrt(CHRA)
            CPFLD(1,2)=1.0/CHPR

            CPFLD(2,1)=1.0/sqrt(CHRA)
            CPFLD(2,2)=1.0
         endif
      else
         if (PARAM(2).lt.0.0) then
            CPFLD(1,1) = -1.0/PARAM(2)
         else
            CPFLD(1,1) = PARAM(2)
         endif

         if (PARAM(1).lt.0.0) then
            CPFLD(1,2) = -1.0/PARAM(1)
         else
            CPFLD(1,2) = PARAM(1)
         endif
      endif

      return
      end
!***********************************************************************
!     copy vectors
      subroutine cht_opcopy (a1,a2,a3,a4,b1,b2,b3,b4)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call copy(a1,b1,ntotv)
         call copy(a2,b2,ntotv)
         if(IF3D) call copy(a3,b3,ntotv)
      endif
      if (IFHEAT) call copy(a4,b4,ntott)

      return
      end
!***********************************************************************
!     subtract vectors A = A-B
      subroutine cht_opsub2 (a1,a2,a3,a4,b1,b2,b3,b4)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call sub2(a1,b1,ntotv)
         call sub2(a2,b2,ntotv)
         if(IF3D) call sub2(a3,b3,ntotv)
      endif
      if (IFHEAT) call sub2(a4,b4,ntott)

      return
      end
!***********************************************************************
!     subtract vectors A = B-C
      subroutine cht_opsub3 (a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
      real c1(1),c2(1),c3(1),c4(1)
!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call sub3(a1,b1,c1,ntotv)
         call sub3(a2,b2,c2,ntotv)
         if(IF3D) call sub3(a3,b3,c3,ntotv)
      endif
      if (IFHEAT) call sub3(a4,b4,c4,ntott)
      return
      end
!***********************************************************************
!     multiply vector by constant A = c*A
      subroutine cht_opcmult (a1,a2,a3,a4,const)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
!     argument list
      real a1(1),a2(1),a3(1),a4(1)
      real const
!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call cmult(a1,const,ntotv)
         call cmult(a2,const,ntotv)
         if(IF3D) call cmult(a3,const,ntotv)
      endif
      if (IFHEAT) call cmult(a4,const,ntott)
      return
      end
!***********************************************************************
!     multiply vector by constant A = c*A with separate const. for
!     velocity and temperature
      subroutine cht_opcmult2c (a1,a2,a3,a4,const1, const2)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
!     argument list
      real a1(1),a2(1),a3(1),a4(1)
      real const1, const2
!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call cmult(a1,const1,ntotv)
         call cmult(a2,const1,ntotv)
         if(IF3D) call cmult(a3,const1,ntotv)
      endif
      if (IFHEAT) call cmult(a4,const2,ntott)
      return
      end
!***********************************************************************
!     add vectors A = A+B
      subroutine cht_opadd2 (a1,a2,a3,a4,b1,b2,b3,b4)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call add2(a1,b1,ntotv)
         call add2(a2,b2,ntotv)
         if(IF3D) call add2(a3,b3,ntotv)
      endif
      if (IFHEAT) call add2(a4,b4,ntott)

      return
      end
!***********************************************************************
!     zero vectors
      subroutine cht_oprzero (a1,a2,a3,a4)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
!     argument list
      real a1(1),a2(1),a3(1),a4(1)
!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call rzero(a1,ntotv)
         call rzero(a2,ntotv)
         if(IF3D) call rzero(a3,ntotv)
      endif
      if (IFHEAT) call rzero(a4,ntott)
      return
      end
!***********************************************************************
!     vector summation with scaling A = A+c*B
      subroutine cht_opadd2cm (a1,a2,a3,a4,b1,b2,b3,b4,c)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
      real c
!     local variables
      integer ntotv, ntott
      integer i
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         if (IF3D) then
            do i=1,ntotv
               a1(i) = a1(i) + b1(i)*c
               a2(i) = a2(i) + b2(i)*c
               a3(i) = a3(i) + b3(i)*c
            enddo
         else
            do i=1,ntotv
               a1(i) = a1(i) + b1(i)*c
               a2(i) = a2(i) + b2(i)*c
            enddo
         endif
      endif
      if (IFHEAT) then
         do i=1,ntott
            a4(i) = a4(i) + b4(i)*c
         enddo
      endif
      return
      end
!***********************************************************************
!     vector subtraction with scaling A = A-c*B
      subroutine cht_opsub2cm (a1,a2,a3,a4,b1,b2,b3,b4,c)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
      real c
!     local variables
      integer ntotv, ntott
      integer i
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         if (IF3D) then
            do i=1,ntotv
               a1(i) = a1(i) - b1(i)*c
               a2(i) = a2(i) - b2(i)*c
               a3(i) = a3(i) - b3(i)*c
            enddo
         else
            do i=1,ntotv
               a1(i) = a1(i) - b1(i)*c
               a2(i) = a2(i) - b2(i)*c
            enddo
         endif
      endif
      if (IFHEAT) then
         do i=1,ntott
            a4(i) = a4(i) - b4(i)*c
         enddo
      endif
      return
      end
!***********************************************************************
!     global scalar
      real function cht_glsc2_wt (b1,b2,b3,b4,x1,x2,x3,x4,wt)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT_DEF'
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
      include 'MASS_DEF'
      include 'MASS'            ! VOLVM1, VOLTM1
      include 'CONHT'           ! CHCST_SC,CHCFF_V,CHCFF_T
!     argument list
      real b1(1),b2(1),b3(1),b4(1),x1(1),x2(1),x3(1),x4(1),wt(1)
!     local variables
      integer ntotv, ntott
      real sum, f1, f2
      integer i
!     functions
      real glsum
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

!     scaling factor velocity vs temperature
!     veorsion for newton
!      f1 = coeff_v
!      f2 = coeff_T
!     version for oic
      f1=CHCFF_V/VOLVM1
      f2=CHCFF_T*CHCST_SC/VOLTM1

      sum = 0.
      if (IFFLOW) then          !if vel
         if (IFHEAT) then       !if temp & vel
            if (IF3D) then
               do i=1,ntotv
                  sum = sum + wt(i)*(f1*(b1(i)*x1(i)+b2(i)*x2(i)
     &                 +b3(i)*x3(i))+f2*b4(i)*x4(i))
               end do
            else
               do i=1,ntotv
                  sum =sum + wt(i)*(f1*(b1(i)*x1(i)+b2(i)*x2(i))
     &                 +f2*b4(i)*x4(i))
               end do
            end if

!     for conjugate heat transfer
            if (ntott.gt.ntotv) then
               do i=ntotv+1,ntott
                  sum = sum + wt(i)*f2*b4(i)*x4(i)
               end do
            end if
        else                   !just vel
           if (IF3D) then
              do i=1,ntotv
                 sum = sum + wt(i)*f1*(b1(i)*x1(i)+
     $                b2(i)*x2(i)+b3(i)*x3(i))
              end do
           else
              do i=1,ntotv
                 sum = sum + wt(i)*f1*(b1(i)*x1(i)+b2(i)*x2(i))
              end do
           end if
        end if
      else                      !just temp
         if (IFHEAT) then
            do i=1,ntott
               sum = sum + wt(i)*(f2*b4(i)*x4(i))
            end do
         end if
      end if
      
      cht_glsc2_wt = glsum(sum,1)
      
      return
      end
!***********************************************************************
!     global scalar
      subroutine cht_weight_fun (lvx,lvy,lvz,lt,coeff)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'MASS_DEF'
      include 'MASS'            ! VOLVM1, VOLTM1
      include 'CONHT'           ! CHCST_SC,CHCFF_V,CHCFF_T
!     argument list
      real lvx(1),lvy(1),lvz(1),lt(1)
      real coeff
!     local variables
      real f1, f2
!-----------------------------------------------------------------------
      f1=CHCFF_V/VOLVM1/coeff
      f2=CHCFF_T*CHCST_SC/VOLTM1/coeff

c     rescale
      call cht_opcmult2c (lvx,lvy,lvz,lt,f1,f2)

      return
      end
