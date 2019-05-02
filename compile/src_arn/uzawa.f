      subroutine incomprn_uzawa (ux,uy,uz,up)
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
c            3.  Mattias, 2015-01-22: This is a modified version of 
c                incomprn that enables computation of the pressure with
c                Uzawa method (without pressure projections), i.e. 
c                E=DH^(-1)D^T in addition to E=(dt/bd(1))DB^(-1)D^T. 
c                See e.g. Fischer (1997), JCP.
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

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

      parameter(nset = 1 + lbelv/lelv)
      common /orthov/ pset(lx2*ly2*lz2*lelv*mxprev,nset)
      common /orthbi/ nprv(2)
      logical ifprjp

      ifprjp=.false.    ! Project out previous pressure solutions?
      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0.and.(.not.IFUZAWA))
     $     ifprjp=.true.

      if (icalld.eq.0) tpres=0.0
      icalld = icalld+1
      npres  = icalld
      etime1 = dnekclock()

      ntot1  = nx1*ny1*nz1*nelv
      ntot2  = nx2*ny2*nz2*nelv

      call opdiv   (dp,ux,uy,uz)

      if (IFUZAWA) then
         intype = -1

         intloc = -1
         call sethlm  (h1,h2,intloc)
         call rzero   (h2inv,ntot1)

         call chsign(dp,ntot2)
      else
         intype = 1

         call rzero   (h1,ntot1)
         call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
         call invers2 (h2inv,h2,ntot1)

         bdti = -bd(1)/dt
         call cmult   (dp,bdti,ntot2)
      endif

      call add2col2(dp,bm2,usrdiv,ntot2) ! User-defined divergence.

      call ortho   (dp)

      i = 1 + ifield/ifldmhd
      if (ifprjp)   call setrhsp  (dp,h1,h2,h2inv,pset(1,i),nprv(i))
      if (IFUZAWA) then
                    call esolver  (dp,h1,h2,h2inv,intype)
      else
                    scaledt = dt/bd(1)
                    scaledi = 1./scaledt
                    call cmult(dp,scaledt,ntot2)        ! scale for tol
                    call esolver  (dp,h1,h2,h2inv,intype)
                    call cmult(dp,scaledi,ntot2)
      endif
      if (ifprjp)   call gensolnp (dp,h1,h2,h2inv,pset(1,i),nprv(i))

      call add2(up,dp,ntot2)

      call opgradt  (w1 ,w2 ,w3 ,dp)
      if (IFUZAWA) then
         call ophinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h1 ,h2 ,tolhv ,nmxh)
         call opadd2  (ux ,uy ,uz ,dv1,dv2,dv3)
      else
         call opbinv   (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
         dtb  = dt/bd(1)
         call opadd2cm (ux ,uy ,uz ,dv1,dv2,dv3, dtb )
      endif

      if (ifmhd)  call chkptol  ! to avoid repetition

      tpres=tpres+(dnekclock()-etime1)

      return
      end
c-----------------------------------------------------------------------
      subroutine incomprp_uzawa (ux,uy,uz,up)
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
c     Comments (Mattias, 2014-12-23):
c     1.  This routine has been modified relative to incomprp to enable 
c         the Uzawa method, i.e. E=DH^(-1)D^T in addition to 
c         E=(dt/bd(1))DB^(-1)D^T. See e.g. Fischer (1997), JCP
c
c     2.  For Uzawa method no pressure projection is implemented!
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

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
      COMMON /SCRCH/ PREXTR(LX2,LY2,LZ2,LELV)
      logical ifprjp
c
      if (icalld.eq.0) tpres=0.0
      icalld=icalld+1
      npres=icalld
c
      ntot1  = nx1*ny1*nz1*nelv
      ntot2  = nx2*ny2*nz2*nelv

      if (IFUZAWA) then
         intype = -1

         intloc = -1
         call sethlm  (h1,h2,intloc)
         call rzero   (h2inv,ntot1)
      else
         intype = 1
         dtbd   = bd(1)/dt

         call rzero   (h1,ntot1)
         call cmult2  (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)
         call invers2 (h2inv,h2,ntot1)
      endif

      call opdiv   (dp,ux,uy,uz)
      call chsign  (dp,ntot2)
      call ortho   (dp)


C******************************************************************


      ifprjp=.false.    ! project out previous pressure solutions?

      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0.and.(.not.IFUZAWA)) 
     $     ifprjp=.true.

      ! Most likely, the following can be commented out. (pff, 1/6/2010)
      if (npert.gt.1.or.ifbase)            ifprjp=.false.

      if (ifprjp)   call setrhs  (dp,h1,h2,h2inv)
                    call esolver (dp,h1,h2,h2inv,intype)
      if (ifprjp)   call gensoln (dp,h1,h2,h2inv)


C******************************************************************

      call opgradt (w1 ,w2 ,w3 ,dp)
      if (IFUZAWA) then
         call ophinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h1 ,h2 ,tolhv ,nmxh)
      else
         call opbinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
      endif
      call opadd2  (ux ,uy ,uz ,dv1,dv2,dv3)
c
      call extrapprp(prextr)
      call lagpresp
      call add3(up,prextr,dp,ntot2)
c
      return
      end
c------------------------------------------------------------------------
