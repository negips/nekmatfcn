!     List of statistics calculation variables


      subroutine STAT_STD_VARS

c============================================================  
!     notation
!     U , V, W - velocity components
!     P - pressure
!     U,X = dU/dX - velocity derivative

!     average U
      lnvar = lnvar + 1
      npos = lnvar
c$$$!     for testing
c$$$      itmp = NX1*NY1*NZ1
c$$$      do i=1,NELV
c$$$         rtmp = STAT_LMAP(i)
c$$$         call cfill(slvel(1,1,1,i,1),rtmp,itmp)
c$$$      enddo
c$$$!     for testing; end
      call stat_compute_1Dav1(slvel(1,1,1,1,1),npos,alpha,beta)

!    lnvar = 1
c------------------------------------------------------------ 

!     average V
      lnvar = lnvar + 1
      npos = lnvar
c$$$!     for testing
c$$$      itmp = NX1*NY1*NZ1
c$$$      do i=1,NELV
c$$$         rtmp = STAT_GMAP(STAT_LMAP(i))
c$$$         call cfill(slvel(1,1,1,i,2),rtmp,itmp)
c$$$      enddo
c$$$!     for testing; end
      call stat_compute_1Dav1(slvel(1,1,1,1,2),npos,alpha,beta)

c------------------------------------------------------------ 

!     average W
      lnvar = lnvar + 1
      npos = lnvar
c$$$!     for testing
c$$$      itmp = NX1*NY1*NZ1
c$$$      do i=1,NELV
c$$$         rtmp = STAT_OWN(STAT_LMAP(i))
c$$$         call cfill(slvel(1,1,1,i,3),rtmp,itmp)
c$$$      enddo
c$$$!     for testing; end
      call stat_compute_1Dav1(slvel(1,1,1,1,3),npos,alpha,beta)

c------------------------------------------------------------ 

!     average P
      lnvar = lnvar + 1
      npos = lnvar
c$$$!     for testing
c$$$      itmp = NX1*NY1*NZ1
c$$$      do i=1,NELV
c$$$         rtmp = NID
c$$$         call cfill(slp(1,1,1,i),rtmp,itmp)
c$$$      enddo
c$$$!     for testing; end
      call stat_compute_1Dav1(slp(1,1,1,1),npos,alpha,beta)

!     get the rest

c============================================================  
! variables below not tested.

!    <uu>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slvel(1,1,1,1,1),slvel(1,1,1,1,1),
     $         npos,alpha,beta)
!npos==5
      
c------------------------------------------------------------ 
!    <vv>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slvel(1,1,1,1,2),slvel(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!    <ww>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slvel(1,1,1,1,3),slvel(1,1,1,1,3),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!    <pp>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),slp(1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!    <uv>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slvel(1,1,1,1,1),slvel(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!    <vw>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slvel(1,1,1,1,2),slvel(1,1,1,1,3),
     $         npos,alpha,beta)
!npos==10
      
c------------------------------------------------------------ 
!    <uw>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slvel(1,1,1,1,1),slvel(1,1,1,1,3),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!    <pu>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),slvel(1,1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!    <pv>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),slvel(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!    <pw>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),slvel(1,1,1,1,3),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!    Variable: pdu/dx. Notation: <p(u_x)>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),dudx(1,1,1,1,1),
     $         npos,alpha,beta)
!npos==15

c------------------------------------------------------------ 
!     <p(u_y)>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),dudx(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <p(u_z)>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),dudx(1,1,1,1,3),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <p(v_x)>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),dvdx(1,1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <p(v_y)>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),dvdx(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <p(v_z)>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),dvdx(1,1,1,1,3),
     $         npos,alpha,beta)
!npos==20

c------------------------------------------------------------ 
!     <p(w_x)>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),dwdx(1,1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <p(w_y)>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),dwdx(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <p(w_z)>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slp(1,1,1,1),dwdx(1,1,1,1,3),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
c------------------------------------------------------------ 

!     UU,VV,WW
      itmp = LX1*LY1*LZ1*LELT*LDIM
      call col3(tmpvel(1,1,1,1,1),slvel(1,1,1,1,1),slvel(1,1,1,1,1),
     $         itmp)
c------------------------------------------------------------ 
c------------------------------------------------------------ 

!     <uuu>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slvel(1,1,1,1,1),tmpvel(1,1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <vvv>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slvel(1,1,1,1,2),tmpvel(1,1,1,1,2),
     $         npos,alpha,beta)
!npos==25

c------------------------------------------------------------ 
!     <www>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(slvel(1,1,1,1,3),tmpvel(1,1,1,1,3),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <uuv>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(tmpvel(1,1,1,1,1),slvel(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <uuw>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(tmpvel(1,1,1,1,1),slvel(1,1,1,1,3),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <vvu>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(tmpvel(1,1,1,1,2),slvel(1,1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <vvw>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(tmpvel(1,1,1,1,2),slvel(1,1,1,1,3),
     $         npos,alpha,beta)
!npos==30

c------------------------------------------------------------ 
!     <wwu>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(tmpvel(1,1,1,1,3),slvel(1,1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <wwv>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(tmpvel(1,1,1,1,3),slvel(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <ppp>
      lnvar = lnvar + 1
      npos = lnvar

!      copy pp to tmppr 
      itmp = LX1*LY1*LZ1*LELT
      call col3(tmppr(1,1,1,1),slp(1,1,1,1),slp(1,1,1,1),
     $         itmp) 
     
      call stat_compute_1Dav2(tmppr(1,1,1,1),slp(1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <pppp>
      lnvar = lnvar + 1
      npos = lnvar

!     tmppr==pp
      call stat_compute_1Dav2(tmppr(1,1,1,1),tmppr(1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <uvw>
      lnvar = lnvar + 1
      npos = lnvar
      
!     copy uv to tmppr (don't need pp anymore) 
      itmp = LX1*LY1*LZ1*LELT
      call col3(tmppr(1,1,1,1),slvel(1,1,1,1,1),slvel(1,1,1,1,2),
     $         itmp) 
      
      call stat_compute_1Dav2(tmppr(1,1,1,1),slvel(1,1,1,1,3),
     $         npos,alpha,beta)
!npos==35

c------------------------------------------------------------ 
!     <uuuu>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <vvvv>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(tmpvel(1,1,1,1,2),tmpvel(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <wwww>
      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(tmpvel(1,1,1,1,3),tmpvel(1,1,1,1,3),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <uuuv>
      lnvar = lnvar + 1
      npos = lnvar

!     tmppr = uv            
      call stat_compute_1Dav2(tmpvel(1,1,1,1,1),tmppr(1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     <uuvv>
      lnvar = lnvar + 1
      npos = lnvar

      call stat_compute_1Dav2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),
     $         npos,alpha,beta)
!npos==40

c------------------------------------------------------------ 
!     <uvvv>
      lnvar = lnvar + 1
      npos = lnvar

!     tmppr = uv      
     
      call stat_compute_1Dav2(tmppr(1,1,1,1),tmpvel(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     Variable: (du/dx)^2 + (du/dy)^2 + (du/dz)^2.
!     Notation: e11

      itmp = LX1*LY1*LZ1*LELT*LDIM
      call col3(tmpvel(1,1,1,1,1),dudx(1,1,1,1,1),dudx(1,1,1,1,1),
     $         itmp)

      itmp = LX1*LY1*LZ1*LELT
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),itmp)
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,3),itmp)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav1(tmpvel(1,1,1,1,1),npos,alpha,beta)

c------------------------------------------------------------ 
!     (dv/dx)^2 + (dv/dy)^2 + (dv/dz)^2.
!     e22

      itmp = LX1*LY1*LZ1*LELT*LDIM
      call col3(tmpvel(1,1,1,1,1),dvdx(1,1,1,1,1),dvdx(1,1,1,1,1),
     $         itmp)

      itmp = LX1*LY1*LZ1*LELT
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),itmp)
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,3),itmp)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav1(tmpvel(1,1,1,1,1),npos,alpha,beta)

c------------------------------------------------------------ 
!     (dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2.
!     e33

      itmp = LX1*LY1*LZ1*LELT*LDIM
      call col3(tmpvel(1,1,1,1,1),dwdx(1,1,1,1,1),dwdx(1,1,1,1,1),
     $         itmp)

      itmp = LX1*LY1*LZ1*LELT
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),itmp)
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,3),itmp)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav1(tmpvel(1,1,1,1,1),npos,alpha,beta)

c------------------------------------------------------------ 
!     (du/dx)*(dv/dx) + (du/dy)*(dv/dy) + (du/dz)*(dv/dz).
!     e12

      itmp = LX1*LY1*LZ1*LELT*LDIM
      call col3(tmpvel(1,1,1,1,1),dudx(1,1,1,1,1),dvdx(1,1,1,1,1),
     $         itmp)

      itmp = LX1*LY1*LZ1*LELT
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),itmp)
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,3),itmp)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav1(tmpvel(1,1,1,1,1),npos,alpha,beta)
!npos==45

c------------------------------------------------------------ 
!     (du/dx)*(dw/dx) + (du/dy)*(dw/dy) + (du/dz)*(dw/dz).
!     e13

      itmp = LX1*LY1*LZ1*LELT*LDIM
      call col3(tmpvel(1,1,1,1,1),dudx(1,1,1,1,1),dwdx(1,1,1,1,1),
     $         itmp)

      itmp = LX1*LY1*LZ1*LELT
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),itmp)
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,3),itmp)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav1(tmpvel(1,1,1,1,1),npos,alpha,beta)

c------------------------------------------------------------ 
!     (dv/dx)*(dw/dx) + (dv/dy)*(dv/dy) + (dv/dz)*(dw/dz).
!     e23

      itmp = LX1*LY1*LZ1*LELT*LDIM
      call col3(tmpvel(1,1,1,1,1),dvdx(1,1,1,1,1),dwdx(1,1,1,1,1),
     $         itmp)

      itmp = LX1*LY1*LZ1*LELT
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),itmp)
      call add2(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,3),itmp)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav1(tmpvel(1,1,1,1,1),npos,alpha,beta)

c------------------------------------------------------------ 
!     (du/dx)*(du/dx)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dudx(1,1,1,1,1),dudx(1,1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     (du/dy)*(du/dy)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dudx(1,1,1,1,2),dudx(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     (du/dx)*(du/dy)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dudx(1,1,1,1,1),dudx(1,1,1,1,2),
     $         npos,alpha,beta)
!npos==50

c------------------------------------------------------------ 
!     (dv/dx)*(dv/dx)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dvdx(1,1,1,1,1),dvdx(1,1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     (dv/dy)*(dv/dy)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dvdx(1,1,1,1,2),dvdx(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     (dv/dx)*(dv/dy)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dvdx(1,1,1,1,1),dvdx(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     (dw/dx)*(dw/dx)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dwdx(1,1,1,1,1),dwdx(1,1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     (dw/dy)*(dw/dy)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dwdx(1,1,1,1,2),dwdx(1,1,1,1,2),
     $         npos,alpha,beta)
!npos==55

c------------------------------------------------------------ 
!     (dw/dx)*(dw/dy)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dwdx(1,1,1,1,1),dwdx(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     (du/dx)*(dv/dx)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dudx(1,1,1,1,1),dvdx(1,1,1,1,1),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     (du/dy)*(dv/dy)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dudx(1,1,1,1,2),dvdx(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     (du/dx)*(dv/dy)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dudx(1,1,1,1,1),dvdx(1,1,1,1,2),
     $         npos,alpha,beta)

c------------------------------------------------------------ 
!     (du/dy)*(dv/dx)

      lnvar = lnvar + 1
      npos = lnvar
      
      call stat_compute_1Dav2(dudx(1,1,1,1,2),dvdx(1,1,1,1,1),
     $         npos,alpha,beta)
!npos==60

c============================================================  

      return
      end subroutine STAT_STD_VARS

