!====================================================================== 
!     Description:Correction for outflow boundary condition 
!                 for the Adjoint problem.
!     Author:     Guillaume Chauvat 
!====================================================================== 
     
      subroutine adjonbc(H2)

      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'
      include 'MASS_DEF'
      include 'MASS'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'ADJOINT_DEF'
      include 'ADJOINT'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'GEOM_DEF'
      include 'GEOM'
      
      real H2(lx1,ly1,lz1,1)
      character cb*3
      integer nfaces,nxyz,nel,ntot,iel,iface,ieg
      integer kx1,kx2,ky1,ky2,kz1,kz2,ix,iy,iz,ia
      logical ifalgn,ifnorx,ifnory,ifnorz
      real robcx,robcy,robcz,robfact
      
      nfaces=2*ndim
      nxyz  =nx1*ny1*nz1
      nel   =nelfld(ifield)
      ntot  =nxyz*nel
c
c     Add diagonal terms to the matrix for adjoint O/o and ON/on
c     boundary conditions. Does nothing for direct equations.
c      
      if (ifadj) then
c
c     check which faces have O/o/ON/on conditions and modify the H2 matrix
c     matrix accordingly. The equations must include a Robin term related
c     to the base flow velocity normal to the wall.
c         
         do iel=1,nel
            do iface=1,nfaces
               ieg=lglel(iel)
               cb = cbc(iface,iel,ifield)
               if (cb.eq.'O  '.or.cb.eq.'o  '
     &              .or.cb.eq.'ON '.or.cb.eq.'on ') then
                  ia = 0
c     
c     Loop on the points of the face
                  call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
                  do iz=kz1,kz2
                     do iy=ky1,ky2
                        do ix=kx1,kx2
                           ia = ia + 1
                           robfact=1./bm1(ix,iy,iz,iel)
                           robfact=robfact*area(ia,1,iface,iel)
c     add an U_n coefficient to H2
                           H2(ix,iy,iz,iel)=H2(ix,iy,iz,iel)+
     &                          area(ia,1,iface,iel)/bm1(ix,iy,iz,iel)*
     &                          (unx(ia,1,iface,iel)*vx(ix,iy,iz,iel)
     &                          +uny(ia,1,iface,iel)*vy(ix,iy,iz,iel)
     &                          +unz(ia,1,iface,iel)*vz(ix,iy,iz,iel))
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end if
         
      end subroutine adjonbc
!----------------------------------------------------------------------
      subroutine adjonbc_2(H2)

      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'
      include 'MASS_DEF'
      include 'MASS'            ! bm1
!      include 'PARALLEL_DEF'
!      include 'PARALLEL'
      include 'SOLN_DEF'
      include 'SOLN'            ! vx,vy,vz
      include 'ADJOINT_DEF'
      include 'ADJOINT'         ! ifadj
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ifield
      include 'INPUT_DEF'
      include 'INPUT'           ! ndim
      include 'GEOM_DEF'
      include 'GEOM'            ! area,unx,uny,unz
      include 'TOPOL_DEF'
      include 'TOPOL'           ! eface1
      
      real h2(lx1,ly1,lz1,lelv)
      character cb*3
      integer nfaces,nel,iel,iface,ieg
      integer ia,ifc
      integer j1,js1,jf1,jskip1
      integer j2,js2,jf2,jskip2
      
      nfaces=2*ndim
      nel   =nelfld(ifield)
c
c     Add diagonal terms to the matrix for adjoint O/o and ON/on
c     boundary conditions. Does nothing for direct equations.
c      
      if (ifadj) then
c       check which faces have O/o/ON/on conditions and modify the H2 matrix
c       matrix accordingly. The equations must include a Robin term related
c       to the base flow velocity normal to the wall.

        call dsset(nx1,ny1,nz1)
        do iel=1,nel
          do ifc=1,nfaces
            iface  = eface1(ifc)
            js1    = skpdat(1,iface)
            jf1    = skpdat(2,iface)
            jskip1 = skpdat(3,iface)
            js2    = skpdat(4,iface)
            jf2    = skpdat(5,iface)
            jskip2 = skpdat(6,iface)

!            ieg=lglel(iel)
            cb = cbc(ifc,iel,ifield)
            if (cb.eq.'O  '.or.cb.eq.'o  '
     $          .or.cb.eq.'ON '.or.cb.eq.'on ') then
              ia    = 0
              do j2=js2,jf2,jskip2
                do j1=js1,jf1,jskip1
                   ia = ia+1
                   if (if3d) then     
                     H2(j1,j2,1,iel)=H2(j1,j2,1,iel)+
     $                    area(ia,1,iface,iel)/bm1(j1,j2,1,iel)*
     $                    (unx(ia,1,iface,iel)*vx(j1,j2,1,iel)
     $                    +uny(ia,1,iface,iel)*vy(j1,j2,1,iel)
     $                    +unz(ia,1,iface,iel)*vz(j1,j2,1,iel))
                   else
                     H2(j1,j2,1,iel)=H2(j1,j2,1,iel)+
     $                    area(ia,1,iface,iel)/bm1(j1,j2,1,iel)*
     $                    (unx(ia,1,iface,iel)*vx(j1,j2,1,iel)
     $                    +uny(ia,1,iface,iel)*vy(j1,j2,1,iel))
                   endif          
                enddo
              enddo
            endif       ! cb
          enddo         ! nfaces
        enddo           ! iel
      endif             ! ifadj
         
      end subroutine adjonbc_2
!---------------------------------------------------------------------- 


 
