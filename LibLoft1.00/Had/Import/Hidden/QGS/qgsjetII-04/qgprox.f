      subroutine qgprox(imode)
c-------------------------------------------------------------------------
c qgprox - propose Pomeron end LC momenta
c imod = 0 - to define normalization
c imod = 1 - propose values according to x^delf * (1 - sum_i x_i)^ahl
c-------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=209,npbmax=1000,npnmax=900,npmax=900,legmax=900)
!      parameter(iapmax=208,npbmax=1000,npnmax=900,npmax=900,legmax=900)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr9/  iwp(iapmax),iwt(iapmax),lqa(iapmax),lqb(iapmax)
     *,iprcn(iapmax),itgcn(iapmax),ias(npbmax),ibs(npbmax),nqs(npbmax)
     *,npompr(npbmax),npomtg(npbmax),npomin(npbmax),nnpr(npmax,npbmax)
     *,nntg(npmax,npbmax),ilpr(legmax,npbmax),iltg(legmax,npbmax)
     *,lnpr(legmax,npbmax),lntg(legmax,npbmax)
     *,nbpi(npnmax,iapmax),nbti(npnmax,iapmax),idnpi(npnmax,iapmax)
     *,idnti(npnmax,iapmax),nppi(npnmax,iapmax),npti(npnmax,iapmax)
     *,nlpi(npnmax,iapmax),nlti(npnmax,iapmax)
      common /qgarr11/ b10
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,delh,sgap
      common /qgarr19/ ahl(3)
      common /qgarr23/ bbpom(npbmax),vvxpom(npbmax)
     *,bpompr(npnmax,iapmax),bpomtg(npnmax,iapmax)
     *,vvxpr(npnmax,iapmax),vvxtg(npnmax,iapmax)
     *,xpompr(npnmax,iapmax),xpomtg(npnmax,iapmax)
     *,xpopin(npmax,npbmax),xpomin(npmax,npbmax),vvxin(npmax,npbmax)
     *,bpomin(npmax,npbmax)
      common /qgarr40/ xppr(npnmax,iapmax),xmtg(npnmax,iapmax)
      common /qgarr43/ moniou
      common /qgdebug/  debug
      external qgran

      if(debug.ge.3)write (moniou,201)imode

      delf=dels
      if(imode.eq.0)then                    !0-configuration (for normalization)
       do ip=1,ia(1)                        !loop over proj. nucleons
        if(lqa(ip).ne.0)then
         do n=1,lqa(ip)                     !loop over proj. constituents
          if(idnpi(n,ip).eq.0)then
           xppr(n,ip)=1.d0/wp0              !LC+ for single Pomeron
          else
           xppr(n,ip)=1.d0/xpompr(n,ip)/scm !LC+ for leg Pomeron
          endif
          enddo
        endif
       enddo
       do it=1,ia(2)                        !loop over targ. nucleons
        if(lqb(it).ne.0)then
         do n=1,lqb(it)                     !loop over targ. constituents
          if(idnti(n,it).eq.0)then
           xmtg(n,it)=1.d0/wm0              !LC- for single Pomeron
          else
           xmtg(n,it)=1.d0/xpomtg(n,it)/scm !LC- for leg Pomeron
          endif
         enddo
        endif
       enddo

      else                                  !proposed configuration
       do ip=1,ia(1)                        !loop over proj. nucleons
        if(lqa(ip).ne.0)then
         xpt=1.d0
         do n=1,lqa(ip)                     !loop over proj. constituents
          nrej=0
          alfl=ahl(icz)+(lqa(ip)-n)*(1.d0+delf)
c          if(icz.eq.2)alfl=alfl-float(lqa(ip)-1)/lqa(ip)  !baryon "junction"
!!!!!!!!!
!          write(0,*) ' 1.d0 + delf=', 1.d0 + delf, ' delf =',delf
!          write(0,*) ' alfl =', alfl, ' lqa(ip) =', lqa(ip), ' n=',n
!          write(0,*) 'ip=',ip, ' ahl(icz)=', ahl(icz), ' icz=',icz
!!!!!!
          gb0=(1.d0-.11d0**(1.d0/(1.d0+delf)))**alfl
     *    *exp(alfl*(1.d0+delf)*.11d0)*2.d0
!!!!!!
!          write(0,*) ' gb0=', gb0
!!!!!!!!!!!
1         continue
c proposal functions are chosen depending on the parameters
c to assure an efficient procedure
          if(delf.ge.0.d0.and.alfl.ge.0.d0
     *    .or.delf.lt.0.d0.and.alfl.le.0.d0)then
           up=1.d0-qgran(b10)**(1.d0/(1.d0+delf))
           if(1.d0-up.lt.1.d-20)goto 1
           tp=1.d0-up**(1.d0/(1.d0+alfl))
!!!!!!!!!!!!!
!           write(0,*) '  1- up= ', 1-up
!!!!!!!!!!!!!!
           gb=(tp/(1.d0-up))**delf
          elseif(delf.lt.0.d0.and.alfl.gt.0.d0)then
           up=-log(1.d0-qgran(b10)*(1.d0-exp(-alfl*(1.d0+delf))))
     *     /alfl/(1.d0+delf)
           tp=up**(1.d0/(1.d0+delf))
           gb=(1.d0-tp)**alfl*exp(alfl*(1.d0+delf)*up)/gb0
          else
           tp=1.d0-qgran(b10)**(1.d0/(1.d0+alfl))
           gb=tp**delf
          endif
          if(qgran(b10).gt.gb)then
           nrej=nrej+1
           goto 1
          endif
          xppr(n,ip)=tp*xpt                 !proposed LC+ for the constituent
          xpt=xpt-xppr(n,ip)                !LC+ of the remnant
          enddo
        endif
       enddo

       do it=1,ia(2)                        !loop over targ. nucleons
        if(lqb(it).ne.0)then
         xmt=1.d0
         do n=1,lqb(it)                     !loop over targ. constituents
          nrej=0
          alfl=ahl(2)+(lqb(it)-n)*(1.d0+delf)
c     *    -float(lqb(it)-1)/lqb(it)                       !baryon "junction"
          gb0=(1.d0-.11d0**(1.d0/(1.d0+delf)))**alfl
     *    *exp(alfl*(1.d0+delf)*.11d0)*2.d0
2         continue
          if(delf.ge.0.d0.and.alfl.ge.0.d0
     *    .or.delf.lt.0.d0.and.alfl.le.0.d0)then
           up=1.d0-qgran(b10)**(1.d0/(1.d0+delf))
           if(1.d0-up.lt.1.d-20)goto 2
           tp=1.d0-up**(1.d0/(1.d0+alfl))
           gb=(tp/(1.d0-up))**delf
          elseif(delf.lt.0.d0.and.alfl.gt.0.d0)then
           up=-log(1.d0-qgran(b10)*(1.d0-exp(-alfl*(1.d0+delf))))
     *     /alfl/(1.d0+delf)
           tp=up**(1.d0/(1.d0+delf))
           gb=(1.d0-tp)**alfl*exp(alfl*(1.d0+delf)*up)/gb0
          else
           tp=1.d0-qgran(b10)**(1.d0/(1.d0+alfl))
           gb=tp**delf
          endif
          if(qgran(b10).gt.gb)then
           nrej=nrej+1
           goto 2
          endif
          if(qgran(b10).gt.gb)goto 2
          xmtg(n,it)=tp*xmt                 !proposed LC- for the constituent
          xmt=xmt-xmtg(n,it)                !LC- of the remnant
          enddo
        endif
       enddo
      endif
      if(debug.ge.4)write (moniou,202)

201   format(2x,'qgprox - propose Pomeron end LC momenta, imode=',i2)
202   format(2x,'qgprox - end')
      return
      end
      subroutine qghot(wpp,wpm,b,vvx,nva,nvb,izp,izt,icdp,icdt,icz,iqq
     *,jpt)
c---------------------------------------------------------------------------
c qghot - semi-hard process
c wpp,wpm   - LC momenta for the constituent partons,
c b         - impact parameter for the semi-hard Pomeron,
c izp, izt  - types of proj. and targ. remnants,
c icdp,icdt - proj. and targ.  diffractive eigenstates,
c iqq - type of the semi-hard process: 0 - gg, 1 - q_vg, 2 - gq_v, 3 - q_vq_v
c jpt=0 - single Pomeron,
c jpt=1 - proj. leg Pomeron,
c jpt=2 - targ. leg Pomeron
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      character*2 tyq
      parameter(njmax=50000)
      dimension ept(4),ep3(4),ey(3),ebal(4),
     *qmin(2),wp(2),iqc(2),iqp(2),nqc(2),ncc(2,2),
     *qv1(30,50),zv1(30,50),qm1(30,50),iqv1(30,50),
     *ldau1(30,49),lpar1(30,50),
     *qv2(30,50),zv2(30,50),qm2(30,50),iqv2(30,50),
     *ldau2(30,49),lpar2(30,50)
      parameter(iapmax=209,npbmax=1000,npnmax=900,npmax=900,legmax=900)
      common /qgarr2/  scm,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr8/  wwm,be(4),dc(5),deta,almpt,ptdif,ptndi
      common /qgarr10/ am(7),ammu
      common /qgarr11/ b10
      common /qgarr12/ nsp
      common /qgarr15/ fp(3),rq(2,3),cd(2,3),gsoft(3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr26/ factk,fqscal
      common /qgarr37/ eqj(4,njmax),iqj(njmax),ncj(2,njmax),nj
      common /qgarr42/ tyq(16)
      common /qgarr43/ moniou
      common /qgarr51/ epsxmn
      common /qgdebug/ debug
      external qgran

      if(debug.ge.1)write (moniou,201)iqq,wpp,wpm,izp,izt,icdp,icdt
     *,icz,jpt,nj

      wwgg=0.d0
      wwqg=0.d0
      wwgq=0.d0
      wwqq=0.d0
      wpi=0.d0
      wmi=0.d0
      sjqg=0.d0
      sjqq=0.d0
      sea1=0.d0
      sea2=0.d0
      glu1=0.d0
      glu2=0.d0
      nj0=nj                       !store number of final partons
      nsp0=nsp                     !store number of final particles

1     sy=wpp*wpm  !energy squared for semi-hard inter. (including preevolution)
      nj=nj0
      nsp=nsp0
      s2min=4.d0*fqscal*qt0       !threshold energy
      if(sy.lt.s2min)stop'qghot: sy<s2min!!!'

      if(iqq.eq.3)then             !q_vq_v-ladder
       wpi=wpp                     !LC+ for the hard interaction
       wmi=wpm                     !LC- for the hard interaction
      else

c-------------------------------------------------
c normalization of acceptance
       xmin=s2min/sy
       iq=(iqq+1)/2+1              !auxilliary type of parton (1 - g, 2 - q(q~))
       sj=qgjit(qt0,qt0,sy,1,iq)   !inclusive parton-parton cross-sections
       if(iqq.eq.0)then
        gb0=-dlog(xmin)*(1.d0-dsqrt(xmin))**(2.d0*betp)*sj
       else
        gb0=(1.d0-xmin)**betp*sj
       endif
       if(jpt.eq.0)then            !single Pomeron
        if(iqq.eq.0)then
         rp0=(rq(icdp,icz)+rq(icdt,2)+alfp*dlog(scm/s2min))
     *   *4.d0*.0389d0
         gb0=gb0/(rq(icdp,icz)+rq(icdt,2)+alfp*dlog(scm/sy))
     *   *exp(-b*b/rp0)
        elseif(iqq.eq.1)then
         rp0=(rq(icdp,icz)+rq(icdt,2)+alfp*dlog(wpp*wm0/s2min))
     *   *4.d0*.0389d0
         gb0=gb0/(rq(icdp,icz)+rq(icdt,2)+alfp*dlog(wm0/wpm))
     *   *exp(-b*b/rp0)
        elseif(iqq.eq.2)then
         rp0=(rq(icdp,icz)+rq(icdt,2)+alfp*dlog(wpm*wp0/s2min))
     *   *4.d0*.0389d0
         gb0=gb0/(rq(icdp,icz)+rq(icdt,2)+alfp*dlog(wp0/wpp))
     *   *exp(-b*b/rp0)
        endif
       elseif(jpt.eq.1)then        !proj. leg Pomeron
        if(iqq.eq.0)then
         rp0=(rq(icdp,icz)+alfp*dlog(wp0*wpm/s2min))*4.d0*.0389d0
         gb0=gb0/(rq(icdp,icz)+alfp*dlog(wp0/wpp))*exp(-b*b/rp0)
        elseif(iqq.eq.1)then
         rp0=(rq(icdp,icz)+alfp*dlog(sy/s2min))*4.d0*.0389d0
         gb0=gb0/rq(icdp,icz)*exp(-b*b/rp0)
        endif
       elseif(jpt.eq.2)then        !targ. leg Pomeron
        if(iqq.eq.0)then
         rp0=(rq(icdt,2)+alfp*dlog(wm0*wpp/s2min))*4.d0*.0389d0
         gb0=gb0/(rq(icdt,2)+alfp*dlog(wm0/wpm))*exp(-b*b/rp0)
        elseif(iqq.eq.2)then
         rp0=(rq(icdt,2)+alfp*dlog(sy/s2min))*4.d0*.0389d0
         gb0=gb0/rq(icdt,2)*exp(-b*b/rp0)
        endif
       endif

c-------------------------------------------------
c sharing of LC momenta between soft preevolution and hard ladder
2      zpm=(1.d0-qgran(b10)*(1.d0-xmin**(delh-dels)))
     * **(1.d0/(delh-dels))
       sjqq=qgjit(qt0,qt0,zpm*sy,2,2)  !inclusive qq cross-section
       sjqg=qgjit(qt0,qt0,zpm*sy,1,2)  !inclusive qg cross-section
       sjgg=qgjit(qt0,qt0,zpm*sy,1,1)  !inclusive gg cross-section

       if(iqq.eq.0)then              !gg-ladder
        xp=zpm**qgran(b10)           !LC+ momentum share
        xm=zpm/xp                    !LC- momentum share
        wpi=wpp*xp                   !LC+ for the hard interaction
        wmi=wpm*xm                   !LC- for the hard interaction
        if(jpt.eq.0)then             !single Pomeron
         rp1=(rq(icdp,icz)+alfp*dlog(wp0/wpi))*4.d0*.0389d0
         rp2=(rq(icdt,2)+alfp*dlog(wm0/wmi))*4.d0*.0389d0
         rp=rp1*rp2/(rp1+rp2)
         z=qgran(b10)
         phi=pi*qgran(b10)
         b0=dsqrt(-rp*dlog(z))
         bb1=(b*rp1/(rp1+rp2)+b0*cos(phi))**2+(b0*sin(phi))**2
         bb2=(b*rp2/(rp1+rp2)-b0*cos(phi))**2+(b0*sin(phi))**2

         xpomr=wpi/wp0
         if(xpomr*sgap.ge.1.d0.or.xpomr*scm.le.sgap)then
          vvx1=0.d0
         else
          v1pnu0=qgfani(1.d0/xpomr,bb1,vvx,0.d0,0.d0,icdp,icz,1)
          v1tnu0=qgfani(xpomr*scm,bb2,vvx,0.d0,0.d0,icdt,2,1)
          nn=0
21        nn=nn+1
          vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx)
          vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx)
          v1pnu=qgfani(1.d0/xpomr,bb1,vvxp,0.d0,0.d0,icdp,icz,1)
          v1tnu=qgfani(xpomr*scm,bb2,vvxt,0.d0,0.d0,icdt,2,1)
          if((abs(v1pnu0-v1pnu).gt.1.d-1.or.abs(v1tnu0-v1tnu).gt.1.d-1)
     *    .and.nn.lt.100)then
           v1pnu0=v1pnu
           v1tnu0=v1tnu
           goto 21
          endif
          vvx1=1.d0-exp(-v1tnu)*(1.d0-vvx)
         endif

         xpomr=wm0/wmi/scm
         if(xpomr*sgap.ge.1.d0.or.xpomr*scm.le.sgap)then
          vvx2=0.d0
         else
          v1pnu0=qgfani(1.d0/xpomr,bb1,vvx,0.d0,0.d0,icdp,icz,1)
          v1tnu0=qgfani(xpomr*scm,bb2,vvx,0.d0,0.d0,icdt,2,1)
          nn=0
22        nn=nn+1
          vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx)
          vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx)
          v1pnu=qgfani(1.d0/xpomr,bb1,vvxp,0.d0,0.d0,icdp,icz,1)
          v1tnu=qgfani(xpomr*scm,bb2,vvxt,0.d0,0.d0,icdt,2,1)
          if((abs(v1pnu0-v1pnu).gt.1.d-1.or.abs(v1tnu0-v1tnu).gt.1.d-1)
     *    .and.nn.lt.100)then
           v1pnu0=v1pnu
           v1tnu0=v1tnu
           goto 22
          endif
          vvx2=1.d0-exp(-v1pnu)*(1.d0-vvx)
         endif

         glu1=qglegc(1.d0/xp,wpp/wp0,bb1,vvx1,icdp,icz,12) !upper gluon PDF
         sea1=qglegc(1.d0/xp,wpp/wp0,bb1,vvx1,icdp,icz,13) !upper quark PDF
         glu2=qglegc(1.d0/xm,wpm/wm0,bb2,vvx2,icdt,2,12)   !lower gluon PDF
         sea2=qglegc(1.d0/xm,wpm/wm0,bb2,vvx2,icdt,2,13)   !lower quark PDF
        elseif(jpt.eq.1)then                         !proj. leg Pomeron
         rp1=(rq(icdp,icz)+alfp*dlog(wp0/wpi))*4.d0*.0389d0
         rp2=-alfp*dlog(xm)*4.d0*.0389d0
         rp=rp1*rp2/(rp1+rp2)
         z=qgran(b10)
         phi=pi*qgran(b10)
         b0=dsqrt(-rp*dlog(z))
         bb1=(b*rp1/(rp1+rp2)+b0*cos(phi))**2+(b0*sin(phi))**2
         bb2=(b*rp2/(rp1+rp2)-b0*cos(phi))**2+(b0*sin(phi))**2

         glu1=qglegc(1.d0/xp,wpp/wp0,bb1,vvx,icdp,icz,12) !upper gluon PDF
         sea1=qglegc(1.d0/xp,wpp/wp0,bb1,vvx,icdp,icz,13) !upper quark PDF
         glu2=qgppdi(xm,0)
         sea2=qgppdi(xm,1)
        elseif(jpt.eq.2)then                         !proj. leg Pomeron
         rp1=(rq(icdt,2)+alfp*dlog(wm0/wmi))*4.d0*.0389d0
         rp2=-alfp*dlog(xp)*4.d0*.0389d0
         rp=rp1*rp2/(rp1+rp2)
         z=qgran(b10)
         phi=pi*qgran(b10)
         b0=dsqrt(-rp*dlog(z))
         bb1=(b*rp1/(rp1+rp2)+b0*cos(phi))**2+(b0*sin(phi))**2
         bb2=(b*rp2/(rp1+rp2)-b0*cos(phi))**2+(b0*sin(phi))**2

         glu1=qglegc(1.d0/xm,wpm/wm0,bb1,vvx,icdt,2,12) !upper gluon PDF
         sea1=qglegc(1.d0/xm,wpm/wm0,bb1,vvx,icdt,2,13) !upper quark PDF
         glu2=qgppdi(xp,0)
         sea2=qgppdi(xp,1)
        endif
        wwgg=glu1*glu2*sjgg
        wwqg=sea1*glu2*sjqg
        wwgq=glu1*sea2*sjqg
        wwqq=sea1*sea2*sjqq
        gbyj=-dlog(zpm)*(wwgg+wwqg+wwgq+wwqq)
        if(jpt.eq.0)then
         rh=rq(icdp,icz)+rq(icdt,2)-alfp*dlog(zpm*sy/scm)
        elseif(jpt.eq.1)then
         rh=rq(icdp,icz)-alfp*dlog(wpp/wp0*zpm)
        elseif(jpt.eq.2)then
         rh=rq(icdt,2)-alfp*dlog(wpm/wm0*zpm)
        else
         rh=0.d0
         stop 'Should not happen in qghot'
        endif
        gbyj=gbyj/rh*exp(-b*b/(4.d0*.0389d0*rh))

       else                          !q_vg-(gq_v-)ladder
        if(iqq.eq.1)then             !q_vg-ladder
         wpi=wpp
         wmi=wpm*zpm
         xm=zpm
         if(jpt.eq.0)then            !single Pomeron
          rp1=rq(icdp,icz)*4.d0*.0389d0
          rp2=(rq(icdt,2)+alfp*dlog(wm0/wmi))*4.d0*.0389d0
          rp=rp1*rp2/(rp1+rp2)
          z=qgran(b10)
          phi=pi*qgran(b10)
          b0=dsqrt(-rp*dlog(z))
          bb1=(b*rp1/(rp1+rp2)+b0*cos(phi))**2+(b0*sin(phi))**2
          bb2=(b*rp2/(rp1+rp2)-b0*cos(phi))**2+(b0*sin(phi))**2

          xpomr=wm0/wmi/scm
          if(xpomr*sgap.ge.1.d0.or.xpomr*scm.le.sgap)then
           vvx2=0.d0
          else
           v1pnu0=qgfani(1.d0/xpomr,bb1,vvx,0.d0,0.d0,icdp,icz,1)
           v1tnu0=qgfani(xpomr*scm,bb2,vvx,0.d0,0.d0,icdt,2,1)
           nn=0
23         nn=nn+1
           vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx)
           vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx)
           v1pnu=qgfani(1.d0/xpomr,bb1,vvxp,0.d0,0.d0,icdp,icz,1)
           v1tnu=qgfani(xpomr*scm,bb2,vvxt,0.d0,0.d0,icdt,2,1)
           if((abs(v1pnu0-v1pnu).gt.1.d-1.or.abs(v1tnu0-v1tnu).gt.1.d-1)
     *     .and.nn.lt.100)then
            v1pnu0=v1pnu
            v1tnu0=v1tnu
            goto 23
           endif
           vvx2=1.d0-exp(-v1pnu)*(1.d0-vvx)
          endif

          glu2=qglegc(1.d0/xm,wpm/wm0,bb2,vvx2,icdt,2,12) !upper gluon PDF
          sea2=qglegc(1.d0/xm,wpm/wm0,bb2,vvx2,icdt,2,13) !upper quark PDF
          wwqg=glu2*sjqg
          wwqq=sea2*sjqq
         else                        !leg Pomeron
          wwqg=qgppdi(xm,0)*sjqg
          wwqq=qgppdi(xm,1)*sjqq
         endif
        elseif(iqq.eq.2)then         !gq_v-ladder
         wpi=wpp*zpm
         wmi=wpm
         xp=zpm
         if(jpt.eq.0)then            !single Pomeron
          rp1=(rq(icdp,icz)+alfp*dlog(wp0/wpi))*4.d0*.0389d0
          rp2=rq(icdt,2)*4.d0*.0389d0
          rp=rp1*rp2/(rp1+rp2)
          z=qgran(b10)
          phi=pi*qgran(b10)
          b0=dsqrt(-rp*dlog(z))
          bb1=(b*rp1/(rp1+rp2)+b0*cos(phi))**2+(b0*sin(phi))**2
          bb2=(b*rp2/(rp1+rp2)-b0*cos(phi))**2+(b0*sin(phi))**2

          xpomr=wpi/wp0
          if(xpomr*sgap.ge.1.d0.or.xpomr*scm.le.sgap)then
           vvx1=0.d0
          else
           v1pnu0=qgfani(1.d0/xpomr,bb1,vvx,0.d0,0.d0,icdp,icz,1)
           v1tnu0=qgfani(xpomr*scm,bb2,vvx,0.d0,0.d0,icdt,2,1)
           nn=0
24         nn=nn+1
           vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx)
           vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx)
           v1pnu=qgfani(1.d0/xpomr,bb1,vvxp,0.d0,0.d0,icdp,icz,1)
           v1tnu=qgfani(xpomr*scm,bb2,vvxt,0.d0,0.d0,icdt,2,1)
           if((abs(v1pnu0-v1pnu).gt.1.d-1.or.abs(v1tnu0-v1tnu).gt.1.d-1)
     *     .and.nn.lt.100)then
            v1pnu0=v1pnu
            v1tnu0=v1tnu
            goto 24
           endif
           vvx1=1.d0-exp(-v1tnu)*(1.d0-vvx)
          endif

          glu1=qglegc(1.d0/xp,wpp/wp0,bb1,vvx1,icdp,icz,12) !upper gluon PDF
          sea1=qglegc(1.d0/xp,wpp/wp0,bb1,vvx1,icdp,icz,13) !upper quark PDF
          wwqg=glu1*sjqg
          wwqq=sea1*sjqq
         else                        !leg Pomeron
          wwqg=qgppdi(xp,0)*sjqg
          wwqq=qgppdi(xp,1)*sjqq
         endif
        endif
        gbyj=wwqg+wwqq
        if(jpt.eq.0)then
         if(iqq.eq.1)then
          rh=rq(icdp,icz)+rq(icdt,2)-alfp*dlog(wpm/wm0*zpm)
         else
          rh=rq(icdp,icz)+rq(icdt,2)-alfp*dlog(wpp/wp0*zpm)
         endif
        elseif(jpt.eq.1)then
         rh=rq(icdp,icz)-alfp*dlog(zpm)
        elseif(jpt.eq.2)then
         rh=rq(icdt,2)-alfp*dlog(zpm)
        else
         rh=0.d0
         stop 'Should not happen in qghot'
        endif
        gbyj=gbyj/rh*exp(-b*b/(4.d0*.0389d0*rh))
       endif

       gbyj=gbyj/gb0/zpm**delh
       if(qgran(b10).gt.gbyj)goto 2
      endif
      if(debug.ge.2)write (moniou,202)wpi*wmi

11    wpi1=wpi
      wmi1=wmi
      wpq=0.d0
      wmq=0.d0
      nj=nj0                     !initialization for the number of final partons
      rrr=qgran(b10)
      jqq=0                                  !gg-ladder
      if(iqq.eq.1.or.iqq.eq.2)then
       if(rrr.lt.wwqq/(wwqg+wwqq))jqq=1      !q_vq_s-laddder
      elseif(iqq.eq.0)then
       if(rrr.lt.wwqg/(wwgg+wwqg+wwgq+wwqq))then
        jqq=1                                !q_sg-ladder
       elseif(rrr.lt.(wwqg+wwgq)/(wwgg+wwqg+wwgq+wwqq))then
        jqq=2                                !gq_s-ladder
       elseif(rrr.lt.(wwqg+wwgq+wwqq)/(wwgg+wwqg+wwgq+wwqq))then
        jqq=3                                !q_sq_s-ladder
       endif
      endif

c-------------------------------------------------
c parton types for the ladder legs and for the leading jets
c iqc(1) - flavor for the upper quark (0 in case of gluon),
c iqc(2) - the same for the lower one
      if(iqq.ne.0.and.iqq.ne.2)then          !q_v from the proj.
       call qgvdef(izp,ic1,ic2,icz)          !leading state flavor
       iqc(1)=ic1                            !upper leg parton
       nj=nj+1
       if(nj.gt.njmax)stop'increase njmax!!!'
       nva=nj
       iqj(nj)=ic2                           !leading jet parton
       ncc(1,1)=nj                           !color connection with leading jet
       ncc(2,1)=0
      else                                   !g(q_s) from the proj.
       nj=nj+1
       if(nj.gt.njmax)stop'increase njmax!!!'
       if(qgran(b10).lt.dc(2))then
        iqj(nj)=-4
       else
        iqj(nj)=-int(2.d0*qgran(b10)+1.d0)
       endif
       iqj(nj+1)=-iqj(nj)
       wp1=wpp-wpi
       wp2=wp1*qgran(b10)
       wp1=wp1-wp2
       eqj(1,nj)=.5d0*wp1
       eqj(2,nj)=eqj(1,nj)
       eqj(3,nj)=0.d0
       eqj(4,nj)=0.d0
       eqj(1,nj+1)=.5d0*wp2
       eqj(2,nj+1)=eqj(1,nj+1)
       eqj(3,nj+1)=0.d0
       eqj(4,nj+1)=0.d0
       if(jqq.eq.0.or.iqq.eq.0.and.jqq.eq.2)then
        iqc(1)=0
        ncc(1,1)=nj
        ncc(2,1)=nj+1
        nj=nj+1
        if(nj.gt.njmax)stop'increase njmax!!!'
       else
        if(qgran(b10).lt..3333d0)then
         iqc(1)=3*(2.d0*int(.5d0+qgran(b10))-1.d0)
        else
         iqc(1)=int(2.d0*qgran(b10)+1.d0)
     *   *(2.d0*int(.5d0+qgran(b10))-1.d0)
        endif
12      zg=xp+qgran(b10)*(1.d0-xp)           !gluon splitting into qq~
        if(qgran(b10).gt.zg**dels*((1.d0-xp/zg)/ (1.d0-xp))**betp)
     *  goto 12
        xg=xp/zg
        wpq0=wpp*(xg-xp)
        wmq=1.d0/wpq0
        wmi1=wmi1-wmq
        if(wmi1*wpi1.le.s2min)goto 11
        nj=nj+2
        if(nj.gt.njmax)stop'increase njmax!!!'
        iqj(nj)=-iqc(1)
        if(iabs(iqc(1)).eq.3)iqj(nj)=iqj(nj)*4/3
        eqj(1,nj)=.5d0*wmq
        eqj(2,nj)=-.5d0*wmq
        eqj(3,nj)=0.d0
        eqj(4,nj)=0.d0
        if(iqc(1).gt.0)then
         ncj(1,nj)=nj-1
         ncj(1,nj-1)=nj
         ncj(2,nj)=0
         ncj(2,nj-1)=0
         ncc(1,1)=nj-2
         ncc(2,1)=0
        else
         ncj(1,nj)=nj-2
         ncj(1,nj-2)=nj
         ncj(2,nj)=0
         ncj(2,nj-2)=0
         ncc(1,1)=nj-1
         ncc(2,1)=0
        endif
       endif
      endif

      if((iqq-2)*(iqq-3)*(iqq-4).eq.0)then     !q_v from the targ.
       call qgvdef(izt,ic1,ic2,2)              !leading state flavor
       iqc(2)=ic1                              !lower leg parton
       nj=nj+1
       if(nj.gt.njmax)stop'increase njmax!!!'
       nvb=nj
       iqj(nj)=ic2
       ncc(1,2)=nj
       ncc(2,2)=0
      else
       nj=nj+1
       if(nj.gt.njmax)stop'increase njmax!!!'
       if(qgran(b10).lt.dc(2))then
        iqj(nj)=-4
       else
        iqj(nj)=-int(2.d0*qgran(b10)+1.d0)
       endif
       iqj(nj+1)=-iqj(nj)
       wm1=wpm-wmi
       wm2=wm1*qgran(b10)
       wm1=wm1-wm2
       eqj(1,nj)=.5d0*wm1
       eqj(2,nj)=-eqj(1,nj)
       eqj(3,nj)=0.d0
       eqj(4,nj)=0.d0
       eqj(1,nj+1)=.5d0*wm2
       eqj(2,nj+1)=-eqj(1,nj+1)
       eqj(3,nj+1)=0.d0
       eqj(4,nj+1)=0.d0
       if(jqq.eq.0.or.iqq.eq.0.and.jqq.eq.1)then
        iqc(2)=0
        ncc(1,2)=nj
        ncc(2,2)=nj+1
        nj=nj+1
        if(nj.gt.njmax)stop'increase njmax!!!'
       else
        if(qgran(b10).lt..3333d0)then
         iqc(2)=3*(2.d0*int(.5d0+qgran(b10))-1.d0)
        else
         iqc(2)=int(2.d0*qgran(b10)+1.d0)
     *   *(2.d0*int(.5d0+qgran(b10))-1.d0)
        endif
14      zg=xm+qgran(b10)*(1.d0-xm)           !gluon splitting into qq~
        if(qgran(b10).gt.zg**dels*((1.d0-xm/zg)/ (1.d0-xm))**betp)
     *  goto 14
        xg=xm/zg
        wmq0=wpm*(xg-xm)
        wpq=1.d0/wmq0
        wpi1=wpi1-wpq
        if(wmi1*wpi1.le.s2min)goto 11
        nj=nj+2
        if(nj.gt.njmax)stop'increase njmax!!!'
        iqj(nj)=-iqc(2)
        if(iabs(iqc(2)).eq.3)iqj(nj)=iqj(nj)*4/3
        eqj(1,nj)=.5d0*wpq
        eqj(2,nj)=.5d0*wpq
        eqj(3,nj)=0.d0
        eqj(4,nj)=0.d0
        if(iqc(2).gt.0)then
         ncj(1,nj)=nj-1
         ncj(1,nj-1)=nj
         ncj(2,nj)=0
         ncj(2,nj-1)=0
         ncc(1,2)=nj-2
         ncc(2,2)=0
        else
         ncj(1,nj)=nj-2
         ncj(1,nj-2)=nj
         ncj(2,nj)=0
         ncj(2,nj-2)=0
         ncc(1,2)=nj-1
         ncc(2,2)=0
        endif
       endif
      endif

      if(jqq.ne.0)then
       if(iqq.ne.0.or.iqq.eq.0.and.jqq.eq.3)then
        sjqq1=qgjit(qt0,qt0,wpi1*wmi1,2,2)
        gbs=sjqq1/sjqq
       else
        sjqg1=qgjit(qt0,qt0,wpi1*wmi1,1,2)
        gbs=sjqg1/sjqg
       endif
       if(qgran(b10).gt.gbs)goto 11
      endif
      wpi=wpi1
      wmi=wmi1

      ept(1)=.5d0*(wpi+wmi)      !ladder 4-momentum
      ept(2)=.5d0*(wpi-wmi)
      ept(3)=0.d0
      ept(4)=0.d0
      qmin(1)=qt0                !q^2 cutoff for the upper leg
      qmin(2)=qt0                !q^2 cutoff for the downer leg
      qminn=max(qmin(1),qmin(2)) !overall q^2 cutoff
      si=qgnrm(ept)
      jini=1
      jj=int(1.5d0+qgran(b10)) !1st parton at upper (jj=1) or downer (jj=2) leg

3     continue

      aaa=qgnrm(ept)             !ladder mass squared
      if(debug.ge.3)write (moniou,203)si,iqc,ept,aaa

      pt2=ept(3)**2+ept(4)**2
      pt=dsqrt(pt2)
      ww=si+pt2

      iqp(1)=min(1,iabs(iqc(1)))+1
      iqp(2)=min(1,iabs(iqc(2)))+1
      wp(1)=ept(1)+ept(2)                 !LC+ for the ladder
      wp(2)=ept(1)-ept(2)                 !LC- for the ladder
      s2min=4.d0*fqscal*qminn   !minimal energy squared for 2-parton production
      if(jini.eq.1)then                   !general ladder
       sj=qgjit(qmin(jj),qmin(3-jj),si,iqp(jj),iqp(3-jj))   !total ladder contribution
       sj1=qgjit1(qmin(3-jj),qmin(jj),si,iqp(3-jj),iqp(jj)) !one-way ordered
       sjb=qgbit(qmin(1),qmin(2),si,iqp(1),iqp(2))          !born contribution
       aks=qgran(b10)
       if(aks.lt.sjb/sj)then
        goto 6      !born process sampled
       elseif(aks.lt.sj1/sj)then       !change to one-way ordered ladder
        jj=3-jj
        sj=sj1
        jini=0
       endif
      else                                !one-way ordered ladder
       sj=qgjit1(qmin(jj),qmin(3-jj),si,iqp(jj),iqp(3-jj)) !one-way ordered
       sjb=qgbit(qmin(1),qmin(2),si,iqp(1),iqp(2))         !born contribution
       if(qgran(b10).lt.sjb/sj)goto 6      !born process sampled
      endif
      wwmin=(s2min+qmin(jj)+pt2-2.d0*pt*dsqrt(qmin(jj)*epsxmn))
     */(1.d0-epsxmn)           !minimal energy squared for 3-parton production

      if(debug.ge.3)write (moniou,204)s2min,wwmin,sj,sjb

      if(ww.lt.1.1d0*wwmin)goto 6         !energy too low -> born process

      xxx=pt*dsqrt(qmin(jj))/ww
      xmin=(s2min+qmin(jj)+pt2)/ww
      xmin=xmin-2.d0*xxx*(xxx+dsqrt(xxx**2+1.d0-xmin))

      xmax=1.d0-epsxmn
      if(debug.ge.3)write (moniou,205)xmin,xmax

      qqmax=(pt*dsqrt(epsxmn)+dsqrt(max(0.d0,pt2*epsxmn
     *+(1.d0+4.d0*fqscal)*(xmax*ww-pt2))))/(1.d0+4.d0*fqscal)
      qqmin=qmin(jj)        !minimal parton virtuality in the current rung
      if(debug.ge.3)write (moniou,206)qqmin,qqmax

      qm0=qqmin
      xm0=xmax
      s2max=xm0*ww

      if(jini.eq.1)then
       sj0=qgjit(qm0,qmin(3-jj),s2max,1,iqp(3-jj))*qgfap(xm0,iqp(jj),1)
     * +qgjit(qm0,qmin(3-jj),s2max,2,iqp(3-jj))*qgfap(xm0,iqp(jj),2)
      else
       sj0=qgjit1(qm0,qmin(3-jj),s2max,1,iqp(3-jj))
     * *qgfap(xm0,iqp(jj),1)
     * +qgjit1(qm0,qmin(3-jj),s2max,2,iqp(3-jj))*qgfap(xm0,iqp(jj),2)
      endif

      gb0=sj0*qm0*qgalf(qm0/alm)*qgsudx(qm0,iqp(jj)) *4.5d0  !normal. of accept.
      if(xm0.le..5d0)then
       gb0=gb0*xm0**(1.d0-delh)
      else
       gb0=gb0*(1.d0-xm0)*2.d0**delh
      endif
      if(debug.ge.3)write (moniou,208)xm0,xmin,xmax,gb0

      xmin2=max(.5d0,xmin)
      xmin1=xmin**delh
      xmax1=min(xmax,.5d0)**delh
      if(xmin.ge..5d0)then                             !choose proposal function
       djl=1.d0
      elseif(xmax.lt..5d0)then
       djl=0.d0
      else
       djl=1.d0/(1.d0+((2.d0*xmin)**delh-1.d0)/delh
     * /dlog(2.d0*(1.d0-xmax)))
      endif

c-------------------------------------------------
c propose x, q^2
4     continue
      if(qgran(b10).gt.djl)then
       x=(xmin1+qgran(b10)*(xmax1-xmin1))**(1.d0/delh) !parton LC share
      else
       x=1.d0-(1.d0-xmin2)*((1.d0-xmax)/(1.d0-xmin2))**qgran(b10)
      endif
      qq=qqmin/(1.d0+qgran(b10)*(qqmin/qqmax-1.d0))    !parton virtuality
      qt2=qq*(1.d0-x)                                  !parton p_t^2
      if(debug.ge.4)write (moniou,209)qq,qqmin,qqmax,x,qt2

      if(qq.gt.qminn)then                  !update virtuality cutoff
       qmin2=qq
      else
       qmin2=qminn
      endif
      qt=dsqrt(qt2)
      call qgcs(c,s)
      ep3(3)=qt*c                          !final parton p_x, p_y
      ep3(4)=qt*s
      pt2new=(ept(3)-ep3(3))**2+(ept(4)-ep3(4))**2!p_t^2 for the remained ladder
      s2min2=max(s2min,4.d0*fqscal*qmin2)  !new ladder kinematic limit
      s2=x*ww-qt2*x/(1.d0-x)-pt2new        !mass squared for the remained ladder
      if(s2.lt.s2min2)goto 4           !ladder mass below threshold -> rejection

      if(jini.eq.1)then                    !weights for g- and q-legs
       sj1=qgjit(qq,qmin(3-jj),s2,1,iqp(3-jj))*qgfap(x,iqp(jj),1)
       sj2=qgjit(qq,qmin(3-jj),s2,2,iqp(3-jj))*qgfap(x,iqp(jj),2)
      else
       sj1=qgjit1(qq,qmin(3-jj),s2,1,iqp(3-jj))*qgfap(x,iqp(jj),1)
       sj2=qgjit1(qq,qmin(3-jj),s2,2,iqp(3-jj))*qgfap(x,iqp(jj),2)
      endif
      gb7=(sj1+sj2)*qgalf(qq/alm)*qq*qgsudx(qq,iqp(jj))/gb0  /2.d0
                               !acceptance probability for x and q**2 simulation
      if(x.le..5d0)then
       gb7=gb7*x**(1.d0-delh)
      else
       gb7=gb7*(1.d0-x)*2.d0**delh
      endif
      if(debug.ge.4)write (moniou,210)gb7,s2,sj1,sj2,jj,jini
      if(qgran(b10).gt.gb7)goto 4          !rejection

c-------------------------------------------------
c define color flow for the emitted jet; perform final state emission
      nqc(2)=0
      if(qgran(b10).lt.sj1/(sj1+sj2))then         !new gluon-leg ladder
       if(iqc(jj).eq.0)then                       !g -> gg
        jt=1
        jq=int(1.5d0+qgran(b10))
        nqc(1)=ncc(jq,jj)                         !color connection for the jet
        nqc(2)=0
       else                                       !q -> qg
        jt=2
        if(iqc(jj).gt.0)then                      !orientation of color flow
         jq=1
        else
         jq=2
        endif
        nqc(1)=0
        ncc(jq,jj)=ncc(1,jj)                      !color connection for the jet
       endif
       iq1=iqc(jj)                                !jet flavor (type)
       iqc(jj)=0                                  !new ladder leg flavor (type)

      else                                        !new quark-leg ladder
       if(iqc(jj).ne.0)then                       !q -> gq
        iq1=0
        jt=3
        if(iqc(jj).gt.0)then                      !orientation of color flow
         jq=1
        else
         jq=2
        endif
        nqc(1)=ncc(1,jj)                          !color connection for the jet
        nqc(2)=0

       else                                       !g -> qq~
        jq=int(1.5d0+qgran(b10))                  !orientation of color flow
        iq1=int(3.d0*qgran(b10)+1.d0)*(3-2*jq)    !jet flavor (type)
        iqc(jj)=-iq1                              !new ladder leg flavor (type)
        jt=4
        nqc(1)=ncc(jq,jj)                         !color connections for the jet
        ncc(1,jj)=ncc(3-jq,jj)
       endif
      endif
      if(debug.ge.3)write (moniou,211)jt

      call qgcjet(qt2,iq1,qv1,zv1,qm1,iqv1,ldau1,lpar1,jq) !final state emission
      si=x*ww-(qt2+qm1(1,1))*x/(1.d0-x)-pt2new  !mass squared for the new ladder
      if(si.gt.s2min2)then
       iq=min(1,iabs(iqc(jj)))+1
       if(jini.eq.1)then
        gb=qgjit(qq,qmin(3-jj),si,iq,iqp(3-jj))
     *  /qgjit(qq,qmin(3-jj),s2,iq,iqp(3-jj))
       else
        gb=qgjit1(qq,qmin(3-jj),si,iq,iqp(3-jj))
     *  /qgjit1(qq,qmin(3-jj),s2,iq,iqp(3-jj))
       endif
       if(qgran(b10).gt.gb)goto 1        !jet mass correction for the acceptance
      else                                        !below threshold -> rejection
       goto 1
      endif

      wp3=wp(jj)*(1.d0-x)
      wm3=(qt2+qm1(1,1))/wp3
      ep3(1)=.5d0*(wp3+wm3)                       !jet 4-momentum
      ep3(2)=.5d0*(wp3-wm3)*(3-2*jj)
      call qgrec(ep3,nqc,qv1,zv1,qm1,iqv1,ldau1,lpar1,jq)
                               !reconstruction of 4-momenta of all final partons
c-------------------------------------------------
c define color connections for the new ladder
      if(jt.eq.1)then
       if(ncc(1,jj).eq.0.and.ncc(2,jj).eq.0)ncc(3-jq,jj)=nqc(1)
       ncc(jq,jj)=nqc(2)
      elseif(jt.eq.2)then
       ncc(3-jq,jj)=nqc(1)
      elseif(jt.eq.3)then
       ncc(1,jj)=nqc(2)
      elseif(jt.eq.4.and.ncc(1,jj).eq.0.and.ncc(2,jj).eq.0)then
       ncc(1,jj)=nqc(1)
      endif

      if(iabs(iq1).eq.3)then
       iqqq=8+iq1/3*4
      else
       iqqq=8+iq1
      endif
      if(debug.ge.3)write (moniou,212)tyq(iqqq),qt2,ep3
      do i=1,4
       ept(i)=ept(i)-ep3(i)                       !new ladder 4-momentum
      enddo
      qmin(jj)=qq                                 !new virtuality cutoffs
      qminn=qmin2
      goto 3                                      !consider next parton emission

c------------------------------------------------
c born process - last parton pair production in the ladder
6     continue
      if(debug.ge.2)write (moniou,214)si,qminn,iqc
      tmin=qminn*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qminn*fqscal/si)))
      qtmin=tmin*(1.d0-tmin/si)
      if(iqc(1).ne.0.or.iqc(2).ne.0)then
       gb0=tmin**2*qgalf(qtmin/fqscal/alm)**2
     * *qgfbor(si,tmin,iqc(1),iqc(2),1)    *1.1d0
      else
       gb0=.25d0*si**2*qgalf(qtmin/fqscal/alm)**2
     * *qgfbor(si,.5d0*si,iqc(1),iqc(2),1)
      endif
      gb0=gb0*qgsudx(qtmin/fqscal,iqp(1))*qgsudx(qtmin/fqscal,iqp(2))
                                                    !normalization of acceptance
      if(debug.ge.3)write (moniou,215)gb0

7     q2=tmin/(1.d0-qgran(b10)*(1.d0-2.d0*tmin/si))   !proposed q^2
      z=q2/si                                         !parton LC momentum share
      qt2=q2*(1.d0-z)                                 !parton p_t^2
      if(qgran(b10).lt..5d0)then
       jm=2
       tq=si-q2
      else
       jm=1
       tq=q2
      endif
      gb=q2**2*qgalf(qt2/fqscal/alm)**2*qgfbor(si,tq,iqc(1),iqc(2),1)
     **qgsudx(qt2/fqscal,iqp(1))*qgsudx(qt2/fqscal,iqp(2))/gb0
                                                      !acceptance probabilty
      if(debug.ge.4)write (moniou,216)gb,q2,z,qt2
      if(qgran(b10).gt.gb)goto 7                      !rejection

c-------------------------------------------------
c define color connections for the 1st emitted jet
      nqc(2)=0
      if(iqc(1).eq.0.and.iqc(2).eq.0)then             !gg-process
       jq=int(1.5d0+qgran(b10))                       !orientation of color flow
       nqc(1)=ncc(jq,jm)

       if(qgran(b10).lt..5d0)then
        jt=1                                          !gg -> gg
        nqc(2)=0
        njc1=ncc(3-jq,jm)                         !color connections for 1st jet
        njc2=ncc(jq,3-jm)
        if(ncc(1,1).eq.0.and.ncc(2,1).eq.0)then
         if(jm.eq.1)nqc(1)=njc2
        else
         if(iqj(njc1).ne.0)then
          ncj(1,njc1)=njc2
         else
          ncj(jq,njc1)=njc2
         endif
         if(iqj(njc2).ne.0)then
          ncj(1,njc2)=njc1
         else
          ncj(3-jq,njc2)=njc1
         endif
        endif
       else                                 !gg -> gg (inverse color connection)
        jt=2
        nqc(2)=ncc(3-jq,3-jm)
       endif

      elseif(iqc(1)*iqc(2).eq.0)then                  !qg -> qg
       if(iqc(1)+iqc(2).gt.0)then                     !orientation of color flow
        jq=1
       else
        jq=2
       endif
       if(qgran(b10).lt..5d0)then
        if(iqc(jm).eq.0)then
         jt=3
         nqc(1)=ncc(jq,jm)
         nqc(2)=0
         njc1=ncc(3-jq,jm)
         njc2=ncc(1,3-jm)
         if(ncc(1,jm).eq.0.and.ncc(2,jm).eq.0)then
          nqc(1)=njc2
         else
          if(iqj(njc1).ne.0)then
           ncj(1,njc1)=njc2
          else
           ncj(jq,njc1)=njc2
          endif
          if(iqj(njc2).ne.0)then
           ncj(1,njc2)=njc1
          else
           ncj(3-jq,njc2)=njc1
          endif
         endif
        else
         jt=4
         nqc(1)=0
         njc1=ncc(1,jm)
         njc2=ncc(3-jq,3-jm)
         if(njc2.ne.0)then
          if(iqj(njc1).ne.0)then
           ncj(1,njc1)=njc2
          else
           ncj(3-jq,njc1)=njc2
          endif
          if(iqj(njc2).ne.0)then
           ncj(1,njc2)=njc1
          else
           ncj(jq,njc2)=njc1
          endif
         endif
        endif
       else
        if(iqc(jm).eq.0)then
         jt=5
         nqc(2)=ncc(3-jq,jm)
         nqc(1)=ncc(1,3-jm)
        else
         jt=6
         nqc(1)=ncc(jq,3-jm)
        endif
       endif

      elseif(iqc(1)*iqc(2).gt.0)then                  !qq (q~q~) -> qq (q~q~)
       jt=7
       if(iqc(1).gt.0)then
        jq=1
       else
        jq=2
       endif
       nqc(1)=ncc(1,3-jm)
      else                                            !qq~ -> qq~
       jt=8
       if(iqc(jm).gt.0)then
        jq=1
       else
        jq=2
       endif
       nqc(1)=0
       njc1=ncc(1,jm)
       njc2=ncc(1,3-jm)
       if(iqj(njc1).ne.0)then
        ncj(1,njc1)=njc2
       else
        ncj(3-jq,njc1)=njc2
       endif
       if(iqj(njc2).ne.0)then
        ncj(1,njc2)=njc1
       else
        ncj(jq,njc2)=njc1
       endif
      endif
      if(jt.ne.8)then
       jq2=jq
      else
       jq2=3-jq
      endif
      if(debug.ge.3)write (moniou,211)jt
      call qgcjet(qt2,iqc(jm),qv1,zv1,qm1,iqv1,ldau1,lpar1,jq)!final state emis.
      call qgcjet(qt2,iqc(3-jm),qv2,zv2,qm2,iqv2,ldau2,lpar2,jq2)
      amt1=qt2+qm1(1,1)
      amt2=qt2+qm2(1,1)
      if(dsqrt(si).gt.dsqrt(amt1)+dsqrt(amt2))then
       z=qgtwd(si,amt1,amt2)
      else
       if(debug.ge.4)write (moniou,217)dsqrt(si),dsqrt(amt1),dsqrt(amt2)
       goto 1                                      !below threshold -> rejection
      endif

      call qgdeft(si,ept,ey)
      wp3=z*dsqrt(si)
      wm3=(qt2+qm1(1,1))/wp3
      ep3(1)=.5d0*(wp3+wm3)                        !1st jet 4-momentum
      ep3(2)=.5d0*(wp3-wm3)
      qt=dsqrt(qt2)
      call qgcs(c,s)
      ep3(3)=qt*c
      ep3(4)=qt*s

      call qgtran(ep3,ey,1)
      call qgrec(ep3,nqc,qv1,zv1,qm1,iqv1,ldau1,lpar1,jq)
                               !reconstruction of 4-momenta of all final partons
      if(iabs(iqc(jm)).eq.3)then
       iqqq=8+iqc(jm)/3*4
      else
       iqqq=8+iqc(jm)
      endif
      if(debug.ge.3)write (moniou,212)tyq(iqqq),qt2,ep3

      wp3=(1.d0-z)*dsqrt(si)
      wm3=(qt2+qm2(1,1))/wp3
      ep3(1)=.5d0*(wp3+wm3)                        !2nd jet 4-momentum
      ep3(2)=.5d0*(wp3-wm3)
      ep3(3)=-qt*c
      ep3(4)=-qt*s
      call qgtran(ep3,ey,1)

c-------------------------------------------------
c define color connections for the 2nd emitted jet
      if(jt.eq.1)then
       nqc(1)=nqc(2)
       if(ncc(1,3-jm).eq.0.and.ncc(2,3-jm).eq.0)then
        nqc(2)=ncc(3-jq,jm)
       else
        nqc(2)=ncc(3-jq,3-jm)
       endif
      elseif(jt.eq.2)then
       if(ncc(1,1).eq.0.and.ncc(2,1).eq.0)then
        if(jm.eq.1)then
         nqc(2)=nqc(1)
         nqc(1)=ncc(jq,3-jm)
        else
         nqc(1)=nqc(2)
         nqc(2)=ncc(3-jq,jm)
        endif
       else
        nqc(2)=ncc(3-jq,jm)
        nqc(1)=ncc(jq,3-jm)
       endif
      elseif(jt.eq.3)then
       nqc(1)=nqc(2)
      elseif(jt.eq.4)then
       nqc(2)=nqc(1)
       if(ncc(1,1).eq.0.and.ncc(2,1).eq.0)then
        nqc(1)=ncc(1,jm)
       else
        nqc(1)=ncc(jq,3-jm)
       endif
      elseif(jt.eq.5)then
       if(ncc(1,jm).eq.0.and.ncc(2,jm).eq.0)then
        nqc(1)=nqc(2)
       else
        nqc(1)=ncc(jq,jm)
       endif
      elseif(jt.eq.6)then
       if(ncc(1,3-jm).eq.0.and.ncc(2,3-jm).eq.0)then
        nqc(2)=nqc(1)
       else
        nqc(2)=ncc(3-jq,3-jm)
       endif
       nqc(1)=ncc(1,jm)
      elseif(jt.eq.7)then
       nqc(1)=ncc(1,jm)
      endif
      call qgrec(ep3,nqc,qv2,zv2,qm2,iqv2,ldau2,lpar2,jq2)
                               !reconstruction of 4-momenta of all final partons
      if(iabs(iqc(3-jm)).eq.3)then
       iqqq=8+iqc(3-jm)/3*4
      else
       iqqq=8+iqc(3-jm)
      endif
      if(debug.ge.3)write (moniou,212)tyq(iqqq),qt2,ep3

      ebal(1)=.5d0*(wpp+wpm)                          !balans of 4-momentum
      ebal(2)=.5d0*(wpp-wpm)
      ebal(3)=0.d0
      ebal(4)=0.d0
      do i=nj0+1,nj
       if(iqq.eq.0.or.iqq.eq.1.and.i.ne.nva.or.iqq.eq.2
     * .and.i.ne.nvb.or.iqq.eq.3.and.i.ne.nva.and.i.ne.nvb)then
        do j=1,4
         ebal(j)=ebal(j)-eqj(j,i)
        enddo
       endif
      enddo
      if(debug.ge.2)write (moniou,218)nj
      if(debug.ge.5)write (moniou,219)ebal
      if(debug.ge.1)write (moniou,220)

201   format(2x,'qghot - semihard interaction:'/
     *4x,'type of the interaction - ',i2/
     *4x,'initial light cone momenta - ',2e10.3/
     *4x,'remnant types - ',2i3,2x,'diffr. eigenstates - ',2i2/
     *4x,'proj. class - ',i2,2x,'Pomeron type - ',i2/
     *4x,'initial number of final partons - ',i4)
202   format(2x,'qghot: mass squared for parton ladder - ',e10.3)
203   format(2x,'qghot: ',' mass squared for the laddder:',e10.3/
     *4x,'ladder end flavors:',2i3/4x,'ladder 5-momentum: ',5e10.3)
204   format(2x,'qghot: kinematic bounds s2min=',e10.3,
     *2x,'wwmin=',e10.3/4x,'jet cross section sj=',e10.3,
     *2x,'born cross section sjb=',e10.3)
205   format(2x,'qghot: xmin=',e10.3,2x,'xmax=',e10.3)
206   format(2x,'qghot: qqmin=',e10.3,2x,'qqmax=',e10.3)
208   format(2x,'qghot: xm0=',e10.3,2x,'xmin=',e10.3,2x,
     *'xmax=',e10.3,2x,'gb0=',e10.3)
209   format(2x,'qghot: qq=',e10.3,2x,'qqmin=',e10.3,2x,
     *'qqmax=',e10.3,2x,'x=',e10.3,2x,'qt2=',e10.3)
210   format(2x,'qghot: gb7=',e10.3,2x,'s2=',e10.3,2x,'sj1=',e10.3
     *,2x,'sj2=',e10.3,2x,'jj=',i2,2x,'jini=',i2)
211   format(2x,'qghot: colour connection jt=:',i1)
212   format(2x,'qghot: new jet flavor:',a2,
     *' pt squared for the jet:',e10.3/4x,'jet 4-momentum:',4e10.3)
214   format(2x,'qghot - highest virtuality subprocess in the ladder:'/
     *4x,'mass squared for the process:',e10.3/4x,'q^2-cutoff:',e10.3
     *,2x,'iqc=',2i3)
215   format(2x,'qghot - normalization of acceptance:',' gb0=',e10.3)
216   format(2x,'qghot - acceptance probabilty:'/
     *4x,'gb=',e10.3,2x,'q2=',e10.3,2x,'z=',e10.3,2x,'qt2=',e10.3)
217   format(2x,'qghot: ecm=',e10.3,2x,'mt1=',e10.3,2x,'mt2=',e10.3)
218   format(2x,'qghot: total number of jets - ',i4)
219   format(2x,'qghot: 4-momentum balans - ',4e10.3)
220   format(2x,'qghot - end')
      return
      end

c------------------------------------------------------------------------
      function npgen(vv,npmin,npmax)
c-----------------------------------------------------------------------
c npgen -  Poisson distribution
c vv    - average number
c npmin - minimal number
c npmax - maximal number
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /qgdebug/  debug
      EXTERNAL qgran

      if(npmin.eq.0)then
       aks=qgran(b10)
       vvn=exp(-vv)
       do n=1,npmax
         aks=aks-vvn
        if(aks.lt.0.d0)goto 1
         vvn=vvn*vv/dble(n)
       enddo
      elseif(npmin.eq.1)then
       aks=qgran(b10)*(1.d0-exp(-vv))
       vvn=exp(-vv)
       do n=1,npmax
         vvn=vvn*vv/dble(n)
         aks=aks-vvn
        if(aks.lt.0.d0)goto 2
       enddo
      elseif(npmin.eq.2)then
       aks=qgran(b10)*(1.d0-exp(-vv)*(1.d0+vv))
       vvn=vv*exp(-vv)
       do n=2,npmax
         vvn=vvn*vv/dble(n)
         aks=aks-vvn
        if(aks.lt.0.d0)goto 2
       enddo
      else
       stop'npgen'
      endif
1     n=n-1
2     npgen=n
      return
      end
