
c=============================================================================
      subroutine qgconf
c-----------------------------------------------------------------------------
c interaction (cut Pomeron) configuration:
c b - impact parameter,
c xa(1-iap,3), xb(1-iat,3) - proj. and targ. nucleon coordinates,
c iddp(1-iap), iddt(1-iat) - proj. and targ. nucleon diffractive eigenstates,
c icona(1-iap) - connection for proj. nucleons (0 if too far from the target),
c iconab(1-iap,1-iat) - connection for proj.-targ. nucleons (0 if too far from
c each other),
c nwp, nwt - numbers of wounded proj. and targ. nucleons (inelastic or diff.),
c iwp(1-iap), iwt(1-iat) - indexes for wounded proj. and targ. nucleons
c (0 - intact, 1 - inel., 2,3 - diffr., -1 - recoiled from diffraction),
c ncola(1-iap), ncolb(1-iat) - index for inel.-wounded proj. and targ. nucleons,
c nbpom  - total number of Pomeron blocks,
c ias(k) (ibs(k)) - index of the proj. (targ.) nucleon for k-th Pomeron block,
c bbpom(k) - squared impact parameter (between proj. and targ.) for k-th block,
c vvxpom(k) - relative strenth of A-screening corrections for k-th block,
c nqs(k) - number of single Pomerons in k-th block (without cut 3P-vertexes),
c npompr(k) - number of proj. leg Pomerons in k-th block,
c npomtg(k) - number of targ. leg Pomerons in k-th block,
c npomin(k) - number of interm. Pomerons (between 2 3P-vertexes) in k-th block,
c xpopin(n,k) - LC momentum of the upper 3P-vertex for n-th interm. Pomeron
c in k-th block,
c xpomin(n,k) - LC momentum of the lower 3P-vertex for n-th interm. Pomeron
c in k-th block,
c nnpr(i,k) - proj. participant index for i-th single Pomeron in k-th block,
c nntg(i,k) - targ. participant index for i-th single Pomeron in k-th block,
c ilpr(i,k) - proj. index for i-th proj. leg Pomeron in k-th block,
c iltg(i,k) - proj. index for i-th targ. leg Pomeron in k-th block,
c lnpr(i,k) - proj. participant index for i-th proj. leg Pomeron in k-th block,
c lntg(i,k) - targ. participant index for i-th targ. leg Pomeron in k-th block,
c lqa(ip) - number of cut Pomerons connected to ip-th proj. nucleon (hadron),
c lqb(it) - number of cut Pomerons connected to it-th targ. nucleon (hadron),
c nbpi(n,ip) - block index for n-th Pomeron connected to ip-th proj. nucleon,
c nbti(n,it) - block index for n-th Pomeron connected to it-th targ. nucleon,
c idnpi(n,ip) - type of n-th Pomeron (0 - single, 1 - leg) connected to ip-th
c proj. nucleon,
c idnti(n,it) - type of n-th Pomeron (0 - single, 1 - leg) connected to it-th
c targ. nucleon,
c nppi(n,ip) - index in the block of n-th Pomeron connected to ip-th proj.
c nucleon (for single Pomerons),
c npti(n,it) - index in the block of n-th Pomeron connected to it-th targ.
c nucleon (for single Pomerons),
c nlpi(n,ip) - index in the block of n-th Pomeron connected to ip-th proj.
c nucleon (for leg Pomerons),
c nlti(n,it) - index in the block of n-th Pomeron connected to it-th targ.
c nucleon (for leg Pomerons),
c iprcn(ip) - index of the recoiled targ. nucleon for ip-th proj. nucleon
c (undergoing diffraction),
c itgcn(it) - index of the recoiled proj. nucleon for it-th targ. nucleon
c (undergoing diffraction),
c bpompr(n,ip) - squared impact parameter for n-th leg Pomeron connected
c to ip-th proj. nucleon,
c bpomtg(n,it) - squared impact parameter for n-th leg Pomeron connected
c to it-th targ. nucleon,
c vvxpr(n,ip) - relative strenth of A-screening corrections for n-th leg
c Pomeron connected to ip-th proj. nucleon,
c vvxtg(n,it) - relative strenth of A-screening corrections for n-th leg
c Pomeron connected to it-th targ. nucleon,
c xpompr(n,ip) - LC momentum of the 3P-vertex for n-th leg Pomeron connected
c to ip-th proj. nucleon,
c xpomtg(n,it) - LC momentum of the 3P-vertex for n-th leg Pomeron connected
c to it-th targ. nucleon
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
!            208---> 209 (KK)
      parameter(iapmax=209,npbmax=1000,npnmax=900,npmax=900,legmax=900)
      dimension xas(iapmax,3),vabs(2),vabsi(2,iapmax),wdifi(iapmax)
     *,vpac(iapmax),vtac(iapmax),xpomip(npmax),xpomim(npmax)
     *,vvxim(npmax),bpomim(npmax),xpompi(legmax),xpomti(legmax)
     *,vvxpi(legmax),vvxti(legmax),bpompi(legmax),bpomti(legmax)
     *,ipompi(legmax),ipomti(legmax),ncola(iapmax),ncolb(iapmax)
     *,wdp(2,iapmax),wdt(2,iapmax),wabs(2,2),xrapmin(100),xrapmax(100)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr4/  ey0(3)
      common /qgarr6/  pi,bm,amws
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),b
      common /qgarr9/  iwp(iapmax),iwt(iapmax),lqa(iapmax),lqb(iapmax)
     *,iprcn(iapmax),itgcn(iapmax),ias(npbmax),ibs(npbmax),nqs(npbmax)
     *,npompr(npbmax),npomtg(npbmax),npomin(npbmax),nnpr(npmax,npbmax)
     *,nntg(npmax,npbmax),ilpr(legmax,npbmax),iltg(legmax,npbmax)
     *,lnpr(legmax,npbmax),lntg(legmax,npbmax)
     *,nbpi(npnmax,iapmax),nbti(npnmax,iapmax),idnpi(npnmax,iapmax)
     *,idnti(npnmax,iapmax),nppi(npnmax,iapmax),npti(npnmax,iapmax)
     *,nlpi(npnmax,iapmax),nlti(npnmax,iapmax)
      common /qgarr10/ am(7),ammu
      common /qgarr11/ b10
      common /qgarr12/ nsp
      common /qgarr13/ nsf,iaf(iapmax)
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,delh,sgap
      common /qgarr23/ bbpom(npbmax),vvxpom(npbmax)
     *,bpompr(npnmax,iapmax),bpomtg(npnmax,iapmax)
     *,vvxpr(npnmax,iapmax),vvxtg(npnmax,iapmax)
     *,xpompr(npnmax,iapmax),xpomtg(npnmax,iapmax)
     *,xpopin(npmax,npbmax),xpomin(npmax,npbmax),vvxin(npmax,npbmax)
     *,bpomin(npmax,npbmax)
      common /qgarr43/ moniou
      common /qgarr46/ iconab(iapmax,iapmax),icona(iapmax)
     *,iconb(iapmax)
      common /qgarr55/ nwt,nwp       !N of wounded targ.(proj.) nucleons
      common /qgarr56/ nspec,nspect  !N of spectators targ.(proj.) nucleons
      common /qgdebug/  debug
      common /qgsIInex1/xan(iapmax,3),xbn(iapmax,3) !used to link with nexus
     *,bqgs,bmaxqgs,bmaxnex,bminnex
      common/jdiff/bdiff,jdiff     !for external use: impact parameter
                                   !for diffraction, diffraction type
ctp from epos
      integer ng1evt,ng2evt,ikoevt
      real    rglevt,sglevt,eglevt,fglevt,typevt
      common/c2evt/ng1evt,ng2evt,rglevt,sglevt,eglevt,fglevt,ikoevt
     *,typevt            !in epos.inc

      external qgran

      if(debug.ge.1)write (moniou,201)
      nsp=0
      nsf=0
      nsp0=nsp

c initialization
1     continue
      do i=1,ia(1)
       iddp(i)=1+int(qgran(b10)+cc(2,icz)) !diffractive eigenstates for proj.
      enddo
      do i=1,ia(2)
       iddt(i)=1+int(qgran(b10)+cc(2,2))   !diffractive eigenstates for targ.
      enddo

c-------------------------------------------------
c squared impact parameter is sampled uniformly (b**2<bm**2)
      b=bm*dsqrt(qgran(b10))
      if(debug.ge.1)write (moniou,202)b

      if(bmaxnex.ge.0.d0)then              !used to link with nexus
       b1=bminnex
       b2=min(bm,bmaxnex)
       if(b1.gt.b2)stop'bmin > bmax in qgsjet'
       b=dsqrt(b1*b1+(b2*b2-b1*b1)*qgran(b10))
       bqgs=b
      endif

c-------------------------------------------------
c nuclear configurations
      if(debug.ge.1)write (moniou,203)
      if(ia(1).gt.1)then          !projectile nucleon coordinates
       call qggea(ia(1),xa,1)     !xa(n,i), i=1,2,3 - x,y,z for n-th nucleon
      else
       do i=1,3
        xa(1,i)=0.d0              !projectile hadron
       enddo
      endif
      if(ia(2).gt.1)then          !target nucleon coordinates
       call qggea(ia(2),xb,2)     !xb(n,i), i=1,2,3 - x,y,z for n-th nucleon
      else
       do i=1,3
        xb(1,i)=0.d0              !target proton
       enddo
      endif

c-------------------------------------------------
c check connections
      if(debug.ge.1)write (moniou,204)
      do it=1,ia(2)
       iconb(it)=0
      enddo

      do ip=1,ia(1)
       icdp=iddp(ip)
       icona(ip)=0
       do it=1,ia(2)
        icdt=iddt(it)
        bbp=(xa(ip,1)+b-xb(it,1))**2+(xa(ip,2)-xb(it,2))**2
        vv1p=qgpomi(scm,bbp,0.d0,0.d0,0.d0,icdp,icdt,icz,1)
        if(vv1p.gt.1.d-3)then
         if(debug.ge.2)write (moniou,205)ip,it
         iconab(ip,it)=1
         icona(ip)=1
         iconb(it)=1
         if(debug.ge.2)write (moniou,206)ip
         if(debug.ge.2)write (moniou,207)it
        else
         iconab(ip,it)=0
        endif
       enddo
      enddo

      nrej=0
2     nrej=nrej+1
      if(debug.ge.2)write (moniou,208)nrej
      if(nrej.gt.10)then
       if(debug.ge.1)write (moniou,209)
       goto 1
      endif
      nsp=nsp0
      nbpom=0
      nwp=0
      nwt=0
      do i=1,ia(1)
       lqa(i)=0
       iwp(i)=0
       ncola(i)=0
       wdp(1,i)=0.d0
       wdp(2,i)=0.d0
      enddo
      do i=1,ia(2)
       lqb(i)=0
       iwt(i)=0
       ncolb(i)=0
       wdt(1,i)=0.d0
       wdt(2,i)=0.d0
      enddo
      nqs(1)=0
      npomin(1)=0
      npompr(1)=0
      npomtg(1)=0

c-------------------------------------------------
c Pomeron configuration
      if(debug.ge.1)write (moniou,210)
      do 4 ip=1,ia(1)             !loop over all projectile nucleons
       if(debug.ge.2)write (moniou,211)ip
       if(icona(ip).eq.0)goto 4
       x=xa(ip,1)+b               !proj. x is shifted by the impact parameter b
       y=xa(ip,2)
       icdp=iddp(ip)              !diffr. eigenstate for ip

       do 3 it=1,ia(2)            !loop over all target nucleons
        if(debug.ge.2)write (moniou,212)it
        if(iconab(ip,it).eq.0)goto 3
        icdt=iddt(it)                         !diffr. eigenstate for it
        bbp=(x-xb(it,1))**2+(y-xb(it,2))**2   !distance squared between ip, it

c calculate nuclear screening factors for "middle point" -> eikonals
        xpomr=1.d0/dsqrt(scm)
        xxp=.5d0*(x+xb(it,1))
        yyp=.5d0*(y+xb(it,2))
        call qgfdf(xxp,yyp,xpomr,vpac,vtac,vvx,vvxp,vvxt,vvxpl,vvxtl
     *  ,ip,it)
        vv=qgpomi(scm,bbp,vvx,vvxp,vvxt,icdp,icdt,icz,1)        !total eikonal
        vv1p=min(vv,qgpomi(scm,bbp,vvx,vvxp,vvxt,icdp,icdt,icz,2)) !1P-eikonal
        if(debug.ge.2)write (moniou,213)vv,vv1p

        if(qgran(b10).gt.1.d0-exp(-2.d0*vv))goto 3 !1.-exp(-2*vv) - probability
                                                   !for inelastic interaction
        iwt(it)=1
        iwp(ip)=1
        ncola(ip)=ncola(ip)+1                   !N of binary collisions for ip
        ncolb(it)=ncolb(it)+1                   !N of binary collisions for it

        n=npgen(2.d0*vv,1,50) !number of elem. inter. for (ip-it) collision
        nbpom=nbpom+1         !new Pomeron block
        if(nbpom.gt.npbmax)then
         goto 2
        endif
        ias(nbpom)=ip         !proj. index for current elementary interaction
        ibs(nbpom)=it         !targ. index for current elementary interaction
        bbpom(nbpom)=bbp      !distance squared between ip, it
        vvxpom(nbpom)=1.d0-(1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt)
        if(debug.ge.2)write (moniou,214)nbpom,ip,it,n

        nqs(nbpom)=0
        npomin(nbpom)=0
        npompr(nbpom)=0
        npomtg(nbpom)=0
        do i=1,n
         if(qgran(b10).lt.vv1p/vv.or.scm.le.sgap**2)then  !single Pomeron
          if(debug.ge.2)write (moniou,215)i
          np=nqs(nbpom)+1
          if(np.gt.legmax)then
           goto 2
          endif
          nqs(nbpom)=np                  !update Pomeron number in the block
          l0=lqa(ip)+1
          if(l0.gt.npnmax)then
           goto 2
          endif
          lqa(ip)=l0                     !update number of connections for proj.
          nnpr(np,nbpom)=l0              !index for connected proj. participant
          nbpi(l0,ip)=nbpom
          idnpi(l0,ip)=0
          nppi(l0,ip)=np
          l0=lqb(it)+1
          if(l0.gt.npnmax)then
           goto 2
          endif
          lqb(it)=l0
          nntg(np,nbpom)=l0              !index for connected targ. participant
          nbti(l0,it)=nbpom
          idnti(l0,it)=0
          npti(l0,it)=np

         else                            !multi-Pomeron vertex
          if(debug.ge.2)write (moniou,219)
          call qg3pdf(vvxpi,vvxti,xpompi,xpomti,bpompi,bpomti,xpomip
     *    ,xpomim,vvxim,bpomim,npompi,npomti,npin,ipompi,ipomti
     *    ,wdp,wdt,ip,it,iret)
          if(iret.ne.0)goto 2

          if(npin.ne.0)then
           if(debug.ge.2)write (moniou,220)i,npin
           npomin(nbpom)=npomin(nbpom)+npin
           if(npomin(nbpom).gt.npmax)then
            goto 2
           endif
           do l=1,npin
            l1=npomin(nbpom)+l-npin
            xpopin(l1,nbpom)=xpomip(l)
            xpomin(l1,nbpom)=xpomim(l)
            vvxin(l1,nbpom)=vvxim(l)
            bpomin(l1,nbpom)=bpomim(l)
           enddo
          endif
          if(npompi.ne.0)then
           if(debug.ge.2)write (moniou,221)i,npompi
           do m=1,npompi
            np=npompr(nbpom)+1
            if(np.gt.legmax)then
             goto 2
            endif
            npompr(nbpom)=np
            ipp=ipompi(m)
            iwp(ipp)=1
            ilpr(np,nbpom)=ipp
            l0=lqa(ipp)+1
            if(l0.gt.npnmax)then
             goto 2
            endif
            lqa(ipp)=l0
            lnpr(np,nbpom)=l0
            nbpi(l0,ipp)=nbpom
            idnpi(l0,ipp)=1
            nlpi(l0,ipp)=np
            vvxpr(l0,ipp)=vvxpi(m)
            xpompr(l0,ipp)=1.d0/xpompi(m)/scm
            bpompr(l0,ipp)=bpompi(m)
           enddo
          endif
          if(npomti.ne.0)then
           if(debug.ge.2)write (moniou,222)i,npomti
           do m=1,npomti
            np=npomtg(nbpom)+1
            if(np.gt.legmax)then
             goto 2
            endif
            npomtg(nbpom)=np
            itt=ipomti(m)
            iwt(itt)=1
            iltg(np,nbpom)=itt
            l0=lqb(itt)+1
            if(l0.gt.npnmax)then
             goto 2
            endif
            lqb(itt)=l0
            lntg(np,nbpom)=l0
            nbti(l0,itt)=nbpom
            idnti(l0,itt)=1
            nlti(l0,itt)=np
            vvxtg(l0,itt)=vvxti(m)
            xpomtg(l0,itt)=xpomti(m)
            bpomtg(l0,itt)=bpomti(m)
           enddo
          endif
         endif
        enddo                   !end of Pomeron loop
3      continue                 !end of it-loop
4     continue                  !end of ip-loop

c-------------------------------------------------
c   diffraction (hadron-hadron case)
      if(ia(1).eq.1.and.ia(2).eq.1.and.iwp(1).eq.0.and.iwt(1).eq.0)then
       wel=0.d0
       winel=0.d0
       do icdp=1,2
       do icdt=1,2
        vv=qgpomi(scm,b*b,0.d0,0.d0,0.d0,icdp,icdt,icz,1)   !total eikonal
        wabs(icdp,icdt)=exp(-vv)
        wel=wel+cc(icdp,icz)*cc(icdt,2)*wabs(icdp,icdt)
        winel=winel+cc(icdp,icz)*cc(icdt,2)*wabs(icdp,icdt)**2
       enddo
       enddo
       if(qgran(b10).le.wel**2/winel)then
        if(debug.ge.1)write (moniou,231)
        goto 1
       endif

       wdifp=cc(1,icz)*cc(2,icz)*(cc(1,2)**2*(wabs(1,1)-wabs(2,1))**2
     * +cc(2,2)**2*(wabs(1,2)-wabs(2,2))**2+2.d0*cc(1,2)*cc(2,2)
     * *(wabs(1,1)-wabs(2,1))*(wabs(1,2)-wabs(2,2)))
       wdift=cc(1,2)*cc(2,2)*(cc(1,icz)**2*(wabs(1,1)-wabs(1,2))**2
     * +cc(2,icz)**2*(wabs(2,1)-wabs(2,2))**2+2.d0*cc(1,icz)*cc(2,icz)
     * *(wabs(1,1)-wabs(1,2))*(wabs(2,1)-wabs(2,2)))
       wdifd=cc(1,icz)*cc(2,icz)*cc(1,2)*cc(2,2)
     * *(wabs(1,1)+wabs(2,2)-wabs(1,2)-wabs(2,1))**2
       aks=(wdifp+wdift+wdifd)*qgran(b10)
       if(aks.lt.wdifp)then
        nwp=nwp+1
        iwp(1)=2
        iprcn(1)=1
        iwt(1)=-1
       elseif(aks.lt.wdifp+wdift)then
        nwt=nwt+1
        iwt(1)=2
        itgcn(1)=1
        iwp(1)=-1
       else
        nwp=nwp+1
        nwt=nwt+1
        iwp(1)=2
        iwt(1)=2
        iprcn(1)=1
        itgcn(1)=1
       endif
       goto 9
      endif

c-------------------------------------------------
c   diffraction (hadron-nucleus & nucleus-nucleus)
      do ip=1,ia(1)             !loop over all projectile nucleons
       x=xa(ip,1)+b             !proj. x is shifted by b
       y=xa(ip,2)
       if(iwp(ip).ne.0)then
        nwp=nwp+1               !one more wounded proj. nucleon
        if(lqa(ip).eq.0.and.(wdp(1,ip).ne.0.d0.or.wdp(2,ip).ne.0.d0))
     *  then
         icdps=iddp(ip)
         xpomr=1.d0/dsqrt(scm)
         do it=1,ia(2)
          if(iconab(ip,it).ne.0)then
            bbp=(x-xb(it,1))**2+(y-xb(it,2))**2
           xxp=.5d0*(x+xb(it,1))
           yyp=.5d0*(y+xb(it,2))
           icdt=iddt(it)
           do icdp=1,2
            iddp(ip)=icdp
            call qgfdf(xxp,yyp,xpomr,vpac,vtac
     *      ,vvx,vvxp,vvxt,vvxpl,vvxtl,ip,it)
            vv=qgpomi(scm,bbp,vvx,vvxp,vvxt,icdp,icdt,icz,1)   !total eikonal
            wdp(icdp,ip)=wdp(icdp,ip)*exp(-vv)
           enddo
          endif
         enddo
         iddp(ip)=icdps
         wdifr=cc(1,icz)*cc(2,icz)*(wdp(1,ip)-wdp(2,ip))**2
     *   /(cc(1,icz)*wdp(1,ip)**2+cc(2,icz)*wdp(2,ip)**2)
         if(qgran(b10).lt.wdifr)iwp(ip)=3                     !LMD excitation
        endif

       elseif(icona(ip).ne.0)then
        if(debug.ge.2)write (moniou,223)ip
        vabs(1)=0.d0
        vabs(2)=0.d0
        icdps=iddp(ip)
        do it=1,ia(2)
          bbp=(x-xb(it,1))**2+(y-xb(it,2))**2
         icdt=iddt(it)
         do icdp=1,2
          if(iconab(ip,it).eq.0)then
           vabsi(icdp,it)=0.d0
          else
           iddp(ip)=icdp
           xpomr=1.d0/dsqrt(scm)
           xxp=.5d0*(x+xb(it,1))
           yyp=.5d0*(y+xb(it,2))
           call qgfdf(xxp,yyp,xpomr,vpac,vtac,vvx,vvxp,vvxt,vvxpl,vvxtl
     *     ,ip,it)
           vv=qgpomi(scm,bbp,vvx,vvxp,vvxt,icdp,icdt,icz,1)   !total eikonal
           vabsi(icdp,it)=vv
           vabs(icdp)=vabs(icdp)+vv
          endif
         enddo
        enddo
        iddp(ip)=icdps
        wdifr=cc(1,icz)*cc(2,icz)*(exp(-vabs(1))-exp(-vabs(2)))**2
     *  /(cc(1,icz)*exp(-2.d0*vabs(1))+cc(2,icz)*exp(-2.d0*vabs(2)))

        if(qgran(b10).lt.wdifr)then       !projectile diffraction
         wdift=0.d0
         do it=1,ia(2)
          if(iwt(it).ne.-1)then
           wdifi(it)=cc(1,icz)*cc(2,icz)*(exp(-vabsi(1,it))
     *     -exp(-vabsi(2,it)))**2/(cc(1,icz)*exp(-2.d0*vabsi(1,it))
     *     +cc(2,icz)*exp(-2.d0*vabsi(2,it)))
           wdift=wdift+wdifi(it)
          else
           wdifi(it)=0.d0
          endif
         enddo
         if(wdift.ne.0.d0)then
          nwp=nwp+1
          iwp(ip)=2
          aks=qgran(b10)*wdift
          do it=1,ia(2)
           aks=aks-wdifi(it)
           if(aks.lt.0.d0)goto 5
          enddo
5          continue
          iprcn(ip)=it
          if(iwt(it).eq.0)iwt(it)=-1
          if(debug.ge.2)write (moniou,224)ip,it
         endif
        endif
       endif
      enddo                            !end of ip-loop

      do 8 it=1,ia(2)                     !check target diffraction
       if(iwt(it).gt.0)then
        nwt=nwt+1                         !one more wounded targ. nucleon
        if(lqb(it).eq.0.and.(wdt(1,it).ne.0.d0.or.wdt(2,it).ne.0.d0))
     *  then
         icdts=iddt(it)
         xpomr=1.d0/dsqrt(scm)
         do ip=1,ia(1)
          if(iconab(ip,it).ne.0)then
           bbp=(xa(ip,1)+b-xb(it,1))**2+(xa(ip,2)-xb(it,2))**2
           xxp=.5d0*(xa(ip,1)+b+xb(it,1))
           yyp=.5d0*(xa(ip,2)+xb(it,2))
           icdp=iddp(ip)
           do icdt=1,2
            iddt(it)=icdt
            call qgfdf(xxp,yyp,xpomr,vpac,vtac
     *      ,vvx,vvxp,vvxt,vvxpl,vvxtl,ip,it)
            vv=qgpomi(scm,bbp,vvx,vvxp,vvxt,icdp,icdt,icz,1)   !total eikonal
             wdt(icdt,it)=wdt(icdt,it)*exp(-vv)
           enddo
          endif
         enddo
         iddt(it)=icdts
         wdifr=cc(1,2)*cc(2,2)*(wdt(1,it)-wdt(2,it))**2
     *   /(cc(1,2)*wdt(1,it)**2+cc(2,2)*wdt(2,it)**2)
         if(qgran(b10).lt.wdifr)iwt(it)=3
        endif

       elseif(iconb(it).ne.0)then
        if(debug.ge.2)write (moniou,225)it
        vabs(1)=0.d0
        vabs(2)=0.d0
        icdts=iddt(it)
        do ip=1,ia(1)
         bbp=(xa(ip,1)+b-xb(it,1))**2+(xa(ip,2)-xb(it,2))**2
         icdp=iddp(ip)
         do icdt=1,2
          if(iconab(ip,it).eq.0)then
           vabsi(icdt,ip)=0.d0
          else
           iddt(it)=icdt
           xpomr=1.d0/dsqrt(scm)
           xxp=.5d0*(xa(ip,1)+b+xb(it,1))
           yyp=.5d0*(xa(ip,2)+xb(it,2))
           call qgfdf(xxp,yyp,xpomr,vpac,vtac,vvx,vvxp,vvxt,vvxpl,vvxtl
     *     ,ip,it)
           vv=qgpomi(scm,bbp,vvx,vvxp,vvxt,icdp,icdt,icz,1)   !total eikonal
           vabsi(icdt,ip)=vv
           vabs(icdt)=vabs(icdt)+vv
          endif
         enddo
        enddo
        iddt(it)=icdts
        wdifr=cc(1,2)*cc(2,2)*(exp(-vabs(1))-exp(-vabs(2)))**2
     *  /(cc(1,2)*exp(-2.d0*vabs(1))+cc(2,2)*exp(-2.d0*vabs(2)))

        if(qgran(b10).lt.wdifr)then       !target diffraction
         wdift=0.d0
         do ip=1,ia(1)
          if(iwp(ip).eq.-1)then
           wdifi(ip)=0.d0
          else
           if(iwp(ip).eq.2)then
            itt=iprcn(ip)
            if(itt.eq.it)goto 7
            if(iwt(itt).eq.2)then
             wdifi(ip)=0.d0
             goto 6
            endif
           endif
           wdifi(ip)=cc(1,2)*cc(2,2)*(exp(-vabsi(1,ip))
     *     -exp(-vabsi(2,ip)))**2/(cc(1,2)*exp(-2.d0*vabsi(1,ip))
     *     +cc(2,2)*exp(-2.d0*vabsi(2,ip)))
          endif
6          wdift=wdift+wdifi(ip)
         enddo
         if(wdift.eq.0.d0)goto 8
         nwt=nwt+1
         iwt(it)=2
         aks=qgran(b10)*wdift
         do ip=1,ia(1)
          aks=aks-wdifi(ip)
          if(aks.lt.0.d0)goto 7
         enddo
7         continue
         itgcn(it)=ip
         if(debug.ge.2)write (moniou,226)it,ip
         if(iwp(ip).eq.0)then
          iwp(ip)=-1
         elseif(iwp(ip).eq.2)then
          itt=iprcn(ip)
          iprcn(ip)=it
          if(itt.ne.it.and.iwt(itt).eq.-1)iwt(itt)=0
         endif
        endif
       endif
8     continue
c check diffractive cross sections          !so060413-beg
9     jdiff=0                             !non-diffractive
      nqst=0
      nint=0
      if(nbpom.ne.0)then
       do i=1,nbpom
        nqst=nqst+nqs(i)
        nint=nint+npomin(i)
       enddo
      endif
      if((nwp.ne.0.or.nwt.ne.0).and.nqst.eq.0)then   !not elastic nor ND
       lqat=0
       do ip=1,ia(1)
        lqat=lqat+lqa(ip)
       enddo
       lqbt=0
       do it=1,ia(2)
        lqbt=lqbt+lqb(it)
       enddo
       iwpt=0
       do ip=1,ia(1)
        if(iwp(ip).eq.1)then
         iwpt=1
         goto 10
        elseif(iwp(ip).ge.2)then
         iwpt=2
        endif
       enddo
10     continue
       iwtt=0
       do it=1,ia(2)
        if(iwt(it).eq.1)then
         iwtt=1
         goto 11
        elseif(iwt(it).ge.2)then
         iwtt=2
        endif
       enddo
11     continue
       if(lqat.eq.0.and.lqbt.eq.0)then
        if(nbpom.eq.0.or.nint.eq.0)then
         if(iwpt.eq.2.and.iwtt.ne.2)then
          jdiff=6                         !SD(LM)-proj
         elseif(iwpt.ne.2.and.iwtt.eq.2)then
          jdiff=7                         !SD(LM)-targ
         elseif(iwpt.eq.2.and.iwtt.eq.2)then
          jdiff=8                         !DD(LM)
         else
          goto 14
         endif
        else
         if(iwpt.ne.2.and.iwtt.ne.2)then
          jdiff=9                         !CD(DPE)
         else
          jdiff=10                        !CD+LMD
         endif
        endif
       elseif(lqat.gt.0.and.lqbt.eq.0.and.iwtt.ne.2)then
        jdiff=1                          !SD(HM)-proj
       elseif(lqat.eq.0.and.lqbt.gt.0.and.iwpt.ne.2)then
        jdiff=2                          !SD(HM)-targ
       elseif(lqat.gt.0.and.lqbt.eq.0.and.iwtt.eq.2)then
        jdiff=3                          !DD(LHM)-proj
       elseif(lqat.eq.0.and.lqbt.gt.0.and.iwpt.eq.2)then
        jdiff=4                          !DD(LHM)-targ

       elseif(lqat.gt.0.and.lqbt.gt.0)then
        if(nbpom.eq.0)stop'problem with nbpom!!!'
        xrapmax(1)=1.d0
        xrapmin(1)=1.d0/scm
        do ibpom=1,nbpom
         if(npompr(ibpom).gt.0)then
          do i=1,npompr(ibpom)
           ip=ilpr(i,ibpom)
           lpom=lnpr(i,ibpom)
           xrapmax(1)=min(xrapmax(1),1.d0/xpompr(lpom,ip)/scm)
          enddo
         endif
         if(npomtg(ibpom).gt.0)then
          do i=1,npomtg(ibpom)
           it=iltg(i,ibpom)
           lpom=lntg(i,ibpom)
           xrapmin(1)=max(xrapmin(1),xpomtg(lpom,it))
          enddo
         endif
        enddo
        if(xrapmin(1).gt..999d0*xrapmax(1))goto 14
        nraps=1
12      if(nraps.gt.90)stop'nraps>90'
        do ibpom=1,nbpom
         if(npomin(ibpom).gt.0)then
          do i=1,npomin(ibpom)
           if(nraps.eq.1)then
            if(1.d0/scm/xpomin(i,ibpom).lt..999d0*xrapmax(1)
     *      .and.xpopin(i,ibpom).gt.1.001d0*xrapmin(1))then !rap-gaps changed
             if(1.d0/scm/xpomin(i,ibpom).lt.1.001d0*xrapmin(1)
     *       .and.xpopin(i,ibpom).gt..999d0*xrapmax(1))then !no rap-gap (filled)
               goto 14
             elseif(xpopin(i,ibpom).gt..999d0*xrapmax(1))then
              xrapmax(1)=1.d0/scm/xpomin(i,ibpom)
             elseif(1.d0/scm/xpomin(i,ibpom).lt.1.001d0*xrapmin(1))then
              xrapmin(1)=xpopin(i,ibpom)
             else
              xrapmin(2)=xrapmin(1)
              xrapmin(1)=xpopin(i,ibpom)
              xrapmax(2)=1.d0/scm/xpomin(i,ibpom)
              nraps=2
              goto 12
             endif
            endif
           else
            if(1.d0/scm/xpomin(i,ibpom).lt..999d0*xrapmax(1)
     *      .and.xpopin(i,ibpom).gt.1.001d0*xrapmin(nraps))then !rap-gaps changed
             if(1.d0/scm/xpomin(i,ibpom).lt.1.001d0*xrapmin(nraps)
     *       .and.xpopin(i,ibpom).gt..999d0*xrapmax(1))then !no rap-gaps (filled)
              goto 14
             else
              do irap=1,nraps
               if(xpopin(i,ibpom).gt..999d0*xrapmax(irap).and.1.d0/scm
     *         /xpomin(i,ibpom).lt.1.001d0*xrapmin(irap))then !gap filled
                if(irap.lt.nraps)then
                 do j=irap,nraps-1
                  xrapmax(j)=xrapmax(j+1)
                  xrapmin(j)=xrapmin(j+1)
                 enddo
                endif
                nraps=nraps-1
                goto 12
               elseif(xpopin(i,ibpom).gt..999d0*xrapmax(irap))then
                xrapmax(irap)=min(1.d0/scm/xpomin(i,ibpom)
     *          ,xrapmax(irap))
               elseif(1.d0/scm/xpomin(i,ibpom)
     *         .lt.1.001d0*xrapmin(irap))then
                xrapmin(irap)=max(xpopin(i,ibpom),xrapmin(irap))
               elseif(1.d0/scm/xpomin(i,ibpom).gt.xrapmin(irap)
     *         .and.xpopin(i,ibpom).lt.xrapmax(irap))then
                xrapmin(irap)=max(xpopin(i,ibpom),xrapmin(irap))
                if(irap.lt.nraps)then
                 do j=1,nraps-irap
                  xrapmax(nraps-j+2)=xrapmax(nraps-j+1)
                  xrapmin(nraps-j+2)=xrapmin(nraps-j+1)
                 enddo
                endif
                xrapmin(irap+1)=xrapmin(irap)
                xrapmin(irap)=xpopin(i,ibpom)
                xrapmax(irap+1)=1.d0/scm/xpomin(i,ibpom)
                nraps=nraps+1
                goto 12
               endif
              enddo                       !end of irap-loop
             endif
            endif
           endif
          enddo                           !end of npin-loop
         endif
        enddo                             !end of ibpom-loop
        jdiff=5                          !DD(HM)
       endif
      endif                              !end of diffr. check
14    bdiff=b

ctp define collision type
      typevt=0                      !no interaction
      if(nwp.gt.0.or.nwt.gt.0)then         !so060413-end
       if(jdiff.eq.0)then                                  !ND (no rap-gaps)
        typevt=1
       elseif(jdiff.eq.8.or.jdiff.eq.10.or.
     *       (jdiff.gt.2.and.jdiff.lt.6))then !DD + (CD+LMD)
        typevt=2
       elseif(jdiff.eq.1.or.jdiff.eq.6)then                  !SD pro
        typevt=4
       elseif(jdiff.eq.2.or.jdiff.eq.7)then                  !SD tar
        typevt=-4
       elseif(jdiff.eq.9)then                                !CD
        typevt=3
       else
        stop'problem with typevt!'
       endif
      endif


c form projectile spectator part
      if(debug.ge.1)write (moniou,227)
      nspec=0
      do ip=1,ia(1)
       if(iwp(ip).eq.0)then
        if(debug.ge.2)write (moniou,228)ip
        nspec=nspec+1
        do l=1,3
         xas(nspec,l)=xa(ip,l)
        enddo
       endif
      enddo

      nspect=0
      do it=1,ia(2)
       if(iwt(it).eq.0)nspect=nspect+1
      enddo

c inelastic interaction: energy sharing and particle production
      if(nwp.ne.0.or.nwt.ne.0)then
       if(ia(1).eq.nspec.or.ia(2).eq.nspect)stop'ia(1)=nspec!!!'
       if(debug.ge.1)write (moniou,229)

       call qgsha(nbpom,ncola,ncolb,iret)
       if(iret.ne.0)goto 1
       if(nsp.le.nsp0+2)then
        if(debug.ge.1)write (moniou,230)
        goto 1
       endif
      else                                 !no interaction
       if(debug.ge.1)write (moniou,231)
       goto 1
      endif
      if(debug.ge.1)write (moniou,232)nsp

c fragmentation of the projectile spectator part
      if(debug.ge.1)write (moniou,233)
      call qgfrgm(nspec,xas)
      if(debug.ge.1)write (moniou,234)nsf
      if(debug.ge.1)write (moniou,235)

201   format(2x,'qgconf - configuration of the interaction')
202   format(2x,'qgconf: impact parameter b=',e10.3,' fm')
203   format(2x,'qgconf: nuclear configurations')
204   format(2x,'qgconf: check connections')
205   format(2x,'qgconf: ',i3,'-th proj. nucleon may interact with '
     *,i3,'-th target nucleon')
206   format(2x,'qgconf: ',i3,'-th projectile nucleon may interact')
207   format(2x,'qgconf: ',i3,'-th target nucleon may interact')
208   format(2x,'qgconf: ',i3,'-th rejection,'
     *,' redo Pomeron configuration')
209   format(2x,'qgconf: too many rejections,'
     *,' redo nuclear configuartions')
210   format(2x,'qgconf: Pomeron configuration')
211   format(2x,'qgconf: check ',i3,'-th projectile nucleon')
212   format(2x,'qgconf: interaction with ',i3,'-th target nucleon?')
213   format(2x,'qgconf: eikonals - total: ',e10.3,2x,'single: ',e10.3)
214   format(2x,'qgconf: ',i4,'-th Pomeron block connected to ',i3
     *,'-th proj. nucleon and'/4x,i3,'-th targ. nucleon;'
     *,' number of element. processes in the block: ',i3)
215   format(2x,'qgconf: ',i3
     *,'-th process in the block is single cut Pomeron')
219   format(2x,'qgconf: configuration of multi-Pomeron vertexes')
220   format(2x,'qgconf: ',i3,'-th process in the block contains '
     *,i3,' interm. Pomerons')
221   format(2x,'qgconf: ',i3,'-th process in the block contains '
     *,i3,' proj. legs')
222   format(2x,'qgconf: ',i3,'-th process in the block contains '
     *,i3,' targ. legs')
223   format(2x,'qgconf: check diffraction for ',i3,'-th proj. nucleon')
224   format(2x,'qgconf: diffr. of ',i3,'-th proj. nucleon,'
     *,' recoil of ',i3,'-th targ. nucleon')
225   format(2x,'qgconf: check diffraction for ',i3,'-th targ. nucleon')
226   format(2x,'qgconf: diffr. of ',i3,'-th targ. nucleon,'
     *,' recoil of ',i3,'-th proj. nucleon')
227   format(2x,'qgconf: projectile spectator part')
228   format(2x,'qgconf: ',i3,'-th proj. nucleon stays idle')
229   format(2x,'qgconf: inelastic interaction: energy sharing'
     *,' and particle production')
230   format(2x,'qgconf: no particle produced - rejection')
231   format(2x,'qgconf: no interaction - rejection')
232   format(2x,'qgconf: ',i5,' particles have been produced')
233   format(2x,'qgconf: fragmentation of the proj. spectator part')
234   format(2x,'qgconf: ',i3,' proj. fragments have been produced')
235   format(2x,'qgconf - end')
      return
      end
