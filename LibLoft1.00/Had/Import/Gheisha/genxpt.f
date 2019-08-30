*cmz :          10/01/91  18.46.53  by  federico carminati
*cmz :  3.14/16 13/03/89  14.48.46  by  nick van eijndhoven (cern)
*-- author :
      subroutine genxpt(ippp,nfl,avern)
c
c *** generation of x- and pt- values for all produced particles ***
c *** nve 02-may-1988 cern geneva ***
c
c origin : h.fesefeldt 11-oct-1987
c
c a simple single variable description e d3s/dp3= f(q) with
c q**2 = (m*x)**2 + pt**2 is used. final state kinematic is produced
c by an ff-type iterative cascade method
c
      common/consts/ pi,twpi,pibtw,mp,mpi,mmu,mel,mkch,mk0,smp,smpi,
     $               smu,ct,ctkch,ctk0,
     $               ml0,msp,ms0,msm,mx0,mxm,ctl0,ctsp,ctsm,ctx0,ctxm,
     $               rmass(35),rcharg(35)
c
                     real mp,mpi,mmu,mel,mkch,mk0,
     *                    ml0,msp,ms0,msm,mx0,mxm
c
      common/curpar/weight(10),ddeltn,ifile,irun,nevt,nevent,shflag,
     *              ithst,ittot,itlst,ifrnd,tofcut,cmom(5),ceng(5),
     *              rs,s,enp(10),np,nm,nn,nr,no,nz,ipa(200),
     *              atno2,zno2
c
      common/result/xend,yend,zend,rca,rce,amas,nch,tof,px,py,pz,
     *              userw,intct,p,en,ek,amasq,deltn,itk,ntk,ipart,ind,
     *              lcalo,icel,sinl,cosl,sinp,cosp,
     *              xold,yold,zold,pold,pxold,pyold,pzold,
     *              xscat,yscat,zscat,pscat,pxscat,pyscat,pzscat
                    real nch,intct
c
      common/mat   / lmat,
     *               den(21),radlth(21),atno(21),zno(21),absl(21),
     *               cden(21),mden(21),x0den(21),x1den(21),rion(21),
     *               matid(21),matid1(21,24),parmat(21,10),
     *               ifrat,ifrac(21),frac1(21,10),den1(21,10),
     *               atno1(21,10),zno1(21,10)
c
      common/event / nsize,ncur,next,ntot,eve(1200)
c
      common/prntfl/inbcd,newbcd,inbin,newbin,npevt,nevtp,lprt,nprt(10)
                    logical lprt,nprt
c
      common/errcom/ ier(100)
c
      common /vecuty/ pv(10,200)
c
c
      common/genin /tecm,amass(18),npg,kgenev
      common/genout/pcm(5,18),wgt
c
c
      real maspar,lamb,nucsup
      dimension maspar(8),bp(8),ptex(8),c1par(5),g1par(5),tavai(2),
     $          side(200),iavai(2),binl(20),dndl(20),twsup(8),
     $          nucsup(6),psup(6),ipax(100)
      dimension rndm(3)
      data maspar/0.75,0.70,0.65,0.60,0.50,0.40,0.75,0.20/
      data     bp/3.50,3.50,3.50,6.00,5.00,4.00,3.50,3.50/
      data   ptex/1.70,1.70,1.50,1.70,1.40,1.20,1.70,1.20/
      data  c1par/0.6,0.6,0.35,0.15,0.10/
      data  g1par/2.6,2.6,1.80,1.30,1.20/
      data binl/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.11,1.25
     $         ,1.43,1.67,2.0,2.5,3.33,5.00,10.00/
      data twsup/1.,1.,0.7,0.5,0.3,0.2,0.1,0.0/
      data nucsup/1.00,0.7,0.5,0.4,0.35,0.3/
      data   psup/3.,6.,20.,50.,100.,1000./
c
c**
c**  for annihilation interactions introduce proper kinematics
c**
      call coranh(nihil,nfl)
c**
c**
c** check first mass-indices
c**
      ek=enp(5)
      en=enp(6)
      p=enp(7)
      s=enp(8)
      rs=enp(9)
      nt=0
      do 1 i=1,100
      if(ipa(i).eq.0) goto 1
      nt=nt+1
      ipa(nt)=ipa(i)
    1 continue
      call vzero(ipa(nt+1),200-nt)
      call ucopy(ipa(1),ipax(1),100)
c**
c** for low multiplicity use two-body resonance model or single/double
c** diffraction model (--> twoclu (--> twob (--> coscat)))
c**
      cfa=0.025*((atno2-1.)/120.)*exp(-(atno2-1.)/120.)
      if(nihil.gt.0) goto 200
      if(nt.ge.8) goto 200
      if(ek.lt.1.) goto 60
      call grndm(rndm,1)
      ran=rndm(1)
      if(ipart.ge.10.and.ipart.le.13.and.ran.lt.0.5) goto 200
      call grndm(rndm,1)
      ran=rndm(1)
      wsup=twsup(nt)
      if(ran.gt.wsup) goto 200
   60 call ucopy(ipax,ipa,100)
      call twoclu(ippp,nfl,avern)
      go to 9999
c**
c** set effective 4-momentum of primary particle
c**
  200 pv(1,199)=p*px
      pv(2,199)=p*py
      pv(3,199)=p*pz
      pv(4,199)=en
      pv(5,199)=amas
      pv(6,199)=nch
      pv(7,199)=tof
      pv(8,199)=ipart
      pv(9,199)=0.
      pv(10,199)=userw
      ier(49)=ier(49)+1
c**
c** some randomisation of order of final state particles
c**
      do 201 i=3,nt
      call grndm(rndm,1)
      ipx=ifix(3.+rndm(1)*(nt-2.))
      if(ipx.gt.nt) ipx=nt
      ipa1=ipa(ipx)
      ipa(ipx)=ipa(i)
  201 ipa(i)  =ipa1
c**
c** distribute in forward and backward hemisphere in cms
c**
      side(1)= 1.
      side(2)=-1.
      ntb=1
      targ=0.
      if(ipart.lt.10.or.ipart.gt.13) goto 53
      call grndm(rndm,1)
      if(rndm(1).lt.0.7) goto 53
      ipa1=ipa(1)
      ipa(1)=ipa(2)
      ipa(2)=ipa1
   53 lead=0
      if(ipart.lt.10.or.ipart.eq.14.or.ipart.eq.16) goto 532
      ipa1=abs(ipa(1))
      if(ipa1.lt.10.or.ipa1.eq.14.or.ipa1.eq.16) goto 531
      lead=ipa1
      goto 532
  531 ipa1=abs(ipa(2))
      if(ipa1.lt.10.or.ipa1.eq.14.or.ipa1.eq.16) goto 532
      lead=ipa1
  532 do 3 i=1,nt
      if(i.le.2) goto 54
      side(i)= 1.
      call grndm(rndm,1)
      if(rndm(1).lt.0.5) side(i)=-1.
      if(side(i).lt.-0.5) ntb=ntb+1
   54 continue
    3 continue
      tb=2.*ntb
      call grndm(rndm,1)
      if(rs.lt.(2.0+rndm(1))) tb=(2.*ntb+nt)/2.
c**
c** add particles from intranuclear cascade
c**
      afc=0.312+0.200*log(log(s))+s**1.5/6000.
      if(afc.gt.0.75) afc=0.75
      xtarg=afc*(atno2**0.33 -1.0)*tb
      if(xtarg.le.0.) xtarg=0.01
      call poisso(xtarg,ntarg)
      nt2=nt+ntarg
      if(nt2.le.180) goto 2
      nt2=180
      ntarg=nt2-nt
    2 continue
c      if (nprt(4)) write(newbcd,3001) ntarg,nt
      nt1=nt+1
      if(ntarg.eq.0) goto 51
c**
c** check number of extra nucleons and pions
c**
      do 881 ipx=1,6
      if(p.le.psup(ipx)) goto 882
  881 continue
      ipx=6
  882 do 4 i=nt1,nt2
      call grndm(rndm,1)
      ran=rndm(1)
      if(ran.lt.nucsup(ipx)) goto 52
      call grndm(rndm,1)
      ipa(i)=-(7+ifix(rndm(1)*3.0))
      goto 4
   52 ipa(i)=-16
      pnrat=1.-zno2/atno2
      call grndm(rndm,1)
      if(rndm(1).gt.pnrat) ipa(i)=-14
      targ=targ+1.
    4 side(i)=-2.
      nt=nt2
c**
c** choose masses and charges for all particles
c**
   51 do 5 i=1,nt
      ipa1=abs(ipa(i))
      pv(5,i)=rmass(ipa1)
      pv(6,i)=rcharg(ipa1)
      pv(7,i)=1.
      if(pv(5,i).lt.0.) pv(7,i)=-1.
      pv(5,i)=abs(pv(5,i))
    5 continue
c**
c** check available kinetic energy, in this model conservation of
c** kinetic energy in forward and backward hemisphere is assumed
c**
    6 if(nt.le.1) goto 60
      tavai(1)=rs/2.
      tavai(2)=(targ+1.)*rs/2.
      iavai(1)=0
      iavai(2)=0
      do 7 i=1,nt
      l=1
      if(side(i).lt.0.) l=2
      iavai(l)=iavai(l)+1
      tavai(l)=tavai(l)-abs(pv(5,i))
    7 continue
      nth=nt
      if(nth.gt.10) nth=10
      if(iavai(1).le.0) goto 60
      if(iavai(2).le.0) goto 60
      if(tavai(1).gt.0.) goto 11
      call grndm(rndm,1)
      iskip=ifix(rndm(1)*(iavai(1)-1))+1
      is=0
      do 10  i=1,nt
      ii=nt-i+1
      if(side(ii).lt.0.) goto 10
      is=is+1
      if(is.ne.iskip) goto 10
      if(ii.eq.nt) goto 9
      nt1=ii+1
      nt2=nt
      do 8 j=nt1,nt2
      ipa(j-1)=ipa(j)
      side(j-1)=side(j)
      do 71 k=1,10
   71 pv(k,j-1)=pv(k,j)
    8 continue
      goto 9
   10 continue
    9 ipa(nt)=0
      side(nt)=0.
      nt=nt-1
      goto 6
   11 if(tavai(2).gt.0.) goto 15
      call grndm(rndm,1)
      iskip=ifix(rndm(1)*(iavai(2)-1))+1
      is=0
      do 14  i=1,nt
      ii=nt-i+1
      if(side(ii).gt.0.) goto 14
      is=is+1
      if(is.ne.iskip) goto 14
      if(side(ii).lt.-1.5) ntarg=ntarg-1
      if(ntarg.lt.0) ntarg=0
      if(ii.eq.nt) goto 13
      nt1=ii+1
      nt2=nt
      do 12 j=nt1,nt2
      ipa(j-1)=ipa(j)
      side(j-1)=side(j)
      do 74 k=1,10
   74 pv(k,j-1)=pv(k,j)
   12 continue
      goto 13
   14 continue
   13 ipa(nt)=0
      side(nt)=0.
      nt=nt-1
      goto 6
   15 if(nt.le.1) goto 60
      if(nt.eq.180) goto 29
      nt1=nt+1
      nt2=180
      do 28 i=nt1,nt2
   28 ipa(i)=0
   29 continue
c**
c** now the preparation is finished.
c** define initial state vectors for lorentz transformations.
c**
      pv(1,181)=0.
      pv(2,181)=0.
      pv(3,181)=p
      pv(4,181)=sqrt(p*p+amas*amas)
      pv(5,181)=abs(amas)
      pv(1,182)=0.
      pv(2,182)=0.
      pv(3,182)=0.
      pv(4,182)=mp
      pv(5,182)=mp
      pv(1,184)=0.
      pv(2,184)=0.
      pv(3,184)=0.
      pv(4,184)=mp*(1.+targ)
      pv(5,184)=pv(4,184)
      pv(1,188)=0.
      pv(2,188)=0.
      pv(3,188)=0.
      pv(1,189)=1.
      pv(2,189)=0.
      pv(3,189)=0.
      call add(181,182,183)
      call add(184,181,184)
      call lor(181,183,181)
      call lor(182,183,182)
c**
c** main loop for 4-momentum generation , see pitha-report (aachen)
c** for a detailed description of the method.
c**
      call grndm(rndm,1)
      phi=rndm(1)*twpi
      ekin1=0.
      ekin2=0.
      do 39 j=1,10
      pv(j,185)=0.
   39 pv(j,186)=0.
      npg=0
      targ1=0.
      do 16 iii=1,nt
      i=nt-iii+1
      ipa1=abs(ipa(i))
c**
c** count number of backward nucleons
c**
      if(i.eq.2) goto 301
      if(side(i).lt.-1.5.and.ipa1.ge.14) goto 301
      goto 38
  301 npg=npg+1
      if(npg.gt.18) goto 38
      side(i)=-3.
      targ1=targ1+1.
      goto 16
   38 j=3
      if(ipa1.lt.14) j=2
      if(ipa1.lt.10) j=1
      if(i.le.2) j=j+3
      if(side(i).lt.-1.5) j=7
      if(j.eq.7.and.ipa1.ge.14) j=8
c**
c** set pt - and phi values, they are changed somewhat in the iteration
c** loop, set mass parameter for lambda fragmentation model
c**
      call grndm(rndm,1)
      ran=rndm(1)
      bpp=bp(j)
      bpe=ptex(j)
      pt2=-log(1.-ran)/bpp
      aspar=maspar(j)
      pt2=pt2**bpe
      pt =sqrt(pt2)
      if(pt.lt.0.001) pt=0.001
      pv(1,i)=pt*cos(phi)
      pv(2,i)=pt*sin(phi)
      pv(10,i)=pt
      binl(1)=0.
      rlmax=1./pv(10,i)
      do 73 j=2,20
   73 binl(j)=rlmax*(j-1)/19.
      et=pv(4,181)
      p0=pv(3,181)
      if(side(i).lt.0.) et=pv(4,182)
      if(side(i).lt.0.) p0=abs(pv(3,182))
      dndl(1)=0.
      ntrial=0
c**
c** start of big iteration loop
c**
   30 ntrial=ntrial+1
      if(ntrial.gt. 2) goto 169
      do 17 l=2,20
      dndl(l)=0.
      x=(binl(l)+binl(l-1))/2.
      if(pv(10,i).lt.0.001) pv(10,i)=0.001
      if(x.gt.1./pv(10,i)) goto 17
      dx=binl(l)-binl(l-1)
      dndl(l)=aspar/sqrt((1.+(aspar*x)**2)**3)
      dndl(l)=et*dndl(l)/sqrt((x*pv(10,i)*et)**2+pv(10,i)**2
     *                             +pv(5,i)**2)
      dndl(l)=dndl(l)*dx
   17 dndl(l)=dndl(l-1)+dndl(l)
      ntri=0
   31 call grndm(rndm,1)
      ran=rndm(1)*dndl(20)
      do 18 l=2,20
      if(ran.lt.dndl(l)) goto 19
   18 continue
c**
c** start of small iteration loop
c**
   19 ntri=ntri+1
      call grndm(rndm,1)
      ran=rndm(1)
      dx=binl(l)-binl(l-1)
      lamb=binl(l-1)+ran*dx/2.
      x=pv(10,i)*lamb
      if(x.gt.1.) x=1.
      x=x*side(i)/abs(side(i))
      pv(3,i)=x*et
      pv(4,i)=pv(3,i)**2+pv(10,i)**2+pv(5,i)**2
      pv(4,i)=sqrt(pv(4,i))
      if(side(i).lt.0.) goto 165
      if(i.gt.2) goto 20
      ekin=tavai(1)-ekin1
      call normal(ran)
      if(ekin.lt.0.) ekin=0.04*abs(ran)
      pv(4,i)=abs(pv(5,i))+ekin
      rnve=abs(pv(4,i)**2-pv(5,i)**2)
      pp=sqrt(rnve)
      call lengtx(i,pp1)
c
      if (pp1 .ge. 1.0e-6) go to 8000
      call grndm(rndm,2)
      rthnve=pi*rndm(1)
      phinve=twpi*rndm(2)
      pv(1,i)=pp*sin(rthnve)*cos(phinve)
      pv(2,i)=pp*sin(rthnve)*sin(phinve)
      pv(3,i)=pp*cos(rthnve)
      go to 8001
 8000 continue
      pv(1,i)=pv(1,i)*pp/pp1
      pv(2,i)=pv(2,i)*pp/pp1
      pv(3,i)=pv(3,i)*pp/pp1
 8001 continue
c
      call add(185,i,185)
      goto 16
   20 ekin=ekin1+pv(4,i)-abs(pv(5,i))
      if(ekin.lt.0.95*tavai(1)) goto 161
      if(ntri.gt. 5) goto 167
      pv(10,i)=pv(10,i)*0.9
      pv( 1,i)=pv( 1,i)*0.9
      pv( 2,i)=pv( 2,i)*0.9
      dndl(20)=dndl(20)*0.9
      if((tavai(2)-abs(pv(5,i))).lt.0.) goto 31
      side(i)=-1.
      tavai(1)=tavai(1)+abs(pv(5,i))
      tavai(2)=tavai(2)-abs(pv(5,i))
      goto 31
  161 call add(185,i,185)
      ekin1=ekin1+pv(4,i)-abs(pv(5,i))
      goto 163
  165 ekin=ekin2+pv(4,i)-abs(pv(5,i))
      xxx=0.95+0.05*targ/20.
      if(xxx.gt.0.999) x=0.999
      if(ekin.lt.xxx*tavai(2)) goto 166
      if(ntri.gt. 5) goto 167
      pv(10,i)=pv(10,i)*0.9
      pv( 1,i)=pv( 1,i)*0.9
      pv( 2,i)=pv( 2,i)*0.9
      dndl(20)=dndl(20)*0.9
      if((tavai(1)-abs(pv(5,i))).lt.0.) goto 31
      side(i)=+1.
      tavai(1)=tavai(1)-abs(pv(5,i))
      tavai(2)=tavai(2)+abs(pv(5,i))
      goto 31
  166 call add(186,i,186)
      ekin2=ekin2+pv(4,i)-abs(pv(5,i))
  163 call add(185,186,187)
      pv(3,187)=0.
      call ang(187,189,cost,phis)
      if(pv(2,187).lt.0.) phis=twpi-phis
      call normal(ran)
      ran=ran*pi/12.
      phi=phis+pi+ran
      if(phi.gt.twpi) phi=phi-twpi
      if(phi.lt.0.) phi=twpi-phi
      goto 16
c**
c** particle momentum zero, reduce kinetic energy of all other
c**
  167 ekin1=0.
      ekin2=0.
      do 162 j=1,10
      pv(j,185)=0.
  162 pv(j,186)=0.
      ii=i+1
      do 168 l=ii,nt
      if(abs(ipa(l)).ge.14.and.side(l).lt.0.) goto 168
      pv(4,l)=pv(4,l)*0.95+0.05*abs(pv(5,l))
      if(pv(4,l).lt.abs(pv(5,l))) pv(4,l)=abs(pv(5,l))
      rnve=abs(pv(4,l)**2-pv(5,l)**2)
      pp=sqrt(rnve)
      call lengtx(l,pp1)
c
      if (pp1 .ge. 1.0e-6) go to 8002
      call grndm(rndm,2)
      rthnve=pi*rndm(1)
      phinve=twpi*rndm(2)
      pv(1,l)=pp*sin(rthnve)*cos(phinve)
      pv(2,l)=pp*sin(rthnve)*sin(phinve)
      pv(3,l)=pp*cos(rthnve)
      go to 8003
 8002 continue
      pv(1,l)=pv(1,l)*pp/pp1
      pv(2,l)=pv(2,l)*pp/pp1
      pv(3,l)=pv(3,l)*pp/pp1
 8003 continue
c
      pv(10,l)=sqrt(pv(1,l)**2+pv(2,l)**2)
      if(side(l).lt.0.) goto 164
      ekin1=ekin1+pv(4,l)-abs(pv(5,l))
      call add(185,l,185)
      goto 168
  164 ekin2=ekin2+pv(4,l)-abs(pv(5,l))
      call add(186,l,186)
  168 continue
c *** next stmt. changed to prevent from infinite looping ***
c*************      goto 38
      go to 30
c**
c** skip particle, if not enough energy
c**
  169 ipa(i)=0
      do 170 j=1,10
  170 pv(j,i)=0.
      goto 163
   16 continue
      ntri=0
      ii=0
      do 320 i=1,nt
      if(ipa(i).eq.0) goto 320
      ii=ii+1
      ipa(ii)=ipa(i)
      side(ii)=side(i)
      do 321 j=1,10
  321 pv(j,ii)=pv(j,i)
  320 continue
      nt=ii
c**
c** backward nucleons produced with a cluster model
c**
      call lor(184,183,187)
      call sub(187,185,187)
      call sub(187,186,187)
      if(targ1.gt.1.5) goto 310
  322 i=2
      call normal(ran)
      ekin=tavai(2)-ekin2
      ekinm=rs/2.-mp
      if(ekin.gt.ekinm) ekin=ekinm
      call normal(ran)
      if(ekin.lt.0.04) ekin=0.04*abs(ran)
      pv(4,i)=abs(pv(5,i))+ekin
      rnve=abs(pv(4,i)**2-pv(5,i)**2)
      pp=sqrt(rnve)
      call lengtx(187,pp1)
c
      if (pp1 .ge. 1.0e-6) go to 8004
      call grndm(rndm,2)
      rthnve=pi*rndm(1)
      phinve=twpi*rndm(2)
      pv(1,i)=pp*sin(rthnve)*cos(phinve)
      pv(2,i)=pp*sin(rthnve)*sin(phinve)
      pv(3,i)=pp*cos(rthnve)
      go to 8005
 8004 continue
      pv(1,i)=pv(1,187)*pp/pp1
      pv(2,i)=pv(2,187)*pp/pp1
      pv(3,i)=pv(3,187)*pp/pp1
 8005 continue
c
      call add(186,i,186)
      goto 330
  310 itarg1=ifix(targ1+0.1)
      if(itarg1.gt.5) itarg1=5
      rmb0=0.
      npg=0
      do 311 i=1,nt
      if(side(i).gt.-2.5) goto 311
      npg=npg+1
      rmb0=rmb0+abs(pv(5,i))
  311 continue
      if(npg.lt.2) goto 322
      call grndm(rndm,1)
      ran=rndm(1)
      rmb=-log(1.-ran)
      gpar=g1par(itarg1)
      cpar=c1par(itarg1)
      rmb=rmb0+rmb**cpar/gpar
      pv(5,187)=rmb
      if(pv(5,187).gt.pv(4,187)) pv(5,187)=pv(4,187)
      rnve=abs(pv(4,187)**2-pv(5,187)**2)
      pp=sqrt(rnve)
      call lengtx(187,pp1)
c
      if (pp1 .ge. 1.0e-6) go to 8006
      call grndm(rndm,2)
      rthnve=pi*rndm(1)
      phinve=twpi*rndm(2)
      pv(1,187)=pp*sin(rthnve)*cos(phinve)
      pv(2,187)=pp*sin(rthnve)*sin(phinve)
      pv(3,187)=pp*cos(rthnve)
      go to 8007
 8006 continue
      pv(1,187)=pv(1,187)*pp/pp1
      pv(2,187)=pv(2,187)*pp/pp1
      pv(3,187)=pv(3,187)*pp/pp1
 8007 continue
c
      i=187
c      if (nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,5)
      pv(1,187)=-pv(1,187)
      pv(2,187)=-pv(2,187)
      pv(3,187)=-pv(3,187)
      kgenev=1
      tecm=pv(5,187)
      npg=0
      do 312 i=1,nt
      if(side(i).gt.-2.5)goto 312
      npg=npg+1
      amass(npg)=abs(pv(5,i))
  312 continue
      call phasp
      npg=0
      do 314 i=1,nt
      if(side(i).gt.-2.5) goto 314
      npg=npg+1
      pv(1,i)=pcm(1,npg)
      pv(2,i)=pcm(2,npg)
      pv(3,i)=pcm(3,npg)
      pv(4,i)=pcm(4,npg)
      call lor(i,187,i)
      call add(186,i,186)
  314 continue
330   continue
      call add(185,186,187)
      ekin1=pv(4,181)+pv(4,182)
      ekin2=pv(4,185)+pv(4,186)
      if(nprt(4)) then
      write(newbcd,2000) ekin1,ekin2
      i=181
      write(newbcd,2001) i,(pv(j,i),j=1,4)
      i=182
      write(newbcd,2001) i,(pv(j,i),j=1,4)
      i=185
      write(newbcd,2001) i,(pv(j,i),j=1,5)
      i=186
      write(newbcd,2001) i,(pv(j,i),j=1,5)
      do 37 i=1,nt
   37 write(newbcd,2001) i,(pv(j,i),j=1,10),ipa(i),side(i)
      endif !nprt(4).ne.0
c**
c** lorentz transformation in lab system
c**
   36 if(nt.le.2) goto 60
      targ=0.
      do 601 i=1,nt
      if(pv(5,i).gt.0.5) targ=targ+1.
      call lor(i,182,i)
  601 continue
      if(targ.lt.0.5) targ=1.
      if(lead.eq.0) goto 6085
      do 6081 i=1,nt
      if(abs(ipa(i)).eq.lead) goto 6085
 6081 continue
      i=1
      if(lead.ge.14.and.abs(ipa(2)).ge.14) i=2
      if(lead.lt.14.and.abs(ipa(2)).lt.14) i=2
      ipa(i)=lead
      ekin=pv(4,i)-abs(pv(5,i))
      pv(5,i)=rmass(lead)
      pv(7,i)=1.
      if(pv(5,i).lt.0.) pv(7,i)=-1.
      pv(5,i)=abs(pv(5,i))
      pv(6,i)=rcharg(lead)
      pv(4,i)=pv(5,i)+ekin
      call lengtx(i,pp)
      rnve=abs(pv(4,i)**2-pv(5,i)**2)
      pp1=sqrt(rnve)
      pv(1,i)=pp1*pv(1,i)/pp
      pv(2,i)=pp1*pv(2,i)/pp
      pv(3,i)=pp1*pv(3,i)/pp
 6085 kgenev=1
      pv(1,184)=0.
      pv(2,184)=0.
      pv(3,184)=p
      pv(4,184)=sqrt(p*p+amas*amas)
      pv(5,184)=abs(amas)
      ekin0=pv(4,184)-pv(5,184)
      pv(1,185)=0.
      pv(2,185)=0.
      pv(3,185)=0.
      pv(4,185)=mp*targ
      pv(5,185)=pv(4,185)
      ekin=pv(4,184)+pv(4,185)
      i=184
c      if (nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,5)
      i=185
c      if (nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,5)
      call add(184,185,186)
      call lor(184,186,184)
      call lor(185,186,185)
      tecm=pv(4,184)+pv(4,185)
      npg=nt
      pv(1,188)=0.
      pv(2,188)=0.
      pv(3,188)=0.
      pv(4,188)=0.
      pv(5,188)=0.
      ekin1=0.
      do 598 i=1,npg
c      if (nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,10),ipa(i),side(i)
      call add(188,i,188)
      ekin1=ekin1+pv(4,i)-pv(5,i)
      ekin=ekin-pv(5,i)
      if(i.gt.18) goto 598
      amass(i)=pv(5,i)
  598 continue
      if(npg.gt.18) goto 597
      call phasp
      ekin=0.
      do 599 i=1,npg
      pv(1,187)=pcm(1,i)
      pv(2,187)=pcm(2,i)
      pv(3,187)=pcm(3,i)
      pv(4,187)=pcm(4,i)
      pv(5,187)=amass(i)
      call lor(187,185,187)
  599 ekin=ekin+pv(4,187)-pv(5,187)
      call ang(188,184,cost,teta)
c      if (nprt(4)) write(newbcd,2003) teta,ekin0,ekin1,ekin
c**
c** make shure, that  kinetic energies are correct.
c** ekin= kinetic energy theoretically
c** ekin1= kinetic energy simulated
c**
  597 if(ekin1.eq.0.) goto 600
      pv(1,187)=0.
      pv(2,187)=0.
      pv(3,187)=0.
      pv(4,187)=0.
      pv(5,187)=0.
      wgt=ekin/ekin1
      ekin1=0.
      do 602 i=1,nt
      ekin=pv(4,i)-pv(5,i)
      ekin=ekin*wgt
      pv(4,i)=ekin+pv(5,i)
      rnve=abs(pv(4,i)**2-pv(5,i)**2)
      pp=sqrt(rnve)
      call lengtx(i,pp1)
c
      if (pp1 .ge. 1.0e-6) go to 8008
      call grndm(rndm,2)
      rthnve=pi*rndm(1)
      phinve=twpi*rndm(2)
      pv(1,i)=pp*sin(rthnve)*cos(phinve)
      pv(2,i)=pp*sin(rthnve)*sin(phinve)
      pv(3,i)=pp*cos(rthnve)
      go to 8009
 8008 continue
      pv(1,i)=pv(1,i)*pp/pp1
      pv(2,i)=pv(2,i)*pp/pp1
      pv(3,i)=pv(3,i)*pp/pp1
 8009 continue
c
      ekin1=ekin1+ekin
      call add(187,i,187)
  602 continue
      call ang(187,184,cost,teta)
c      if (nprt(4)) write(newbcd,2003) teta,ekin0,ekin1
c**
c** rotate in direction of z-axis, this does disturb in some way our
c** inclusive distributions, but it is nessacary for momentum conser-
c** vation.
c**
  600 pv(1,187)=0.
      pv(2,187)=0.
      pv(3,187)=0.
      pv(4,187)=0.
      pv(5,187)=0.
      do 596 i=1,nt
  596 call add(187,i,187)
c**
c** some smearing in transverse direction from fermi motion
c**
*          call rannor(ran1,ran2)
      call grndm(rndm,2)
      ry=rndm(1)
      rz=rndm(2)
      rx=6.283185*rz
      a1=sqrt(-2.*log(ry))
      ran1=a1*sin(rx)
      ran2=a1*cos(rx)
      pv(1,187)=pv(1,187)+ran1*0.020*targ
      pv(2,187)=pv(2,187)+ran2*0.020*targ
      call defs(184,187,188)
      pv(1,187)=0.
      pv(2,187)=0.
      pv(3,187)=0.
      pv(4,187)=0.
      pv(5,187)=0.
      do 595 i=1,nt
      call trac(i,188,i)
  595 call add(187,i,187)
      call ang(187,184,cost,teta)
c      if (nprt(4)) write(newbcd,2003) teta
c**
c** rotate in direction of primary particle, subtract binding energies
c** and make some further corrections if required (steep, steeq)
c**
      dekin=0.
      npions=0
      ek1=0.
      ek2=0.
      do 21 i=1,nt
      call defs1(i,199,i)
c      if (nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,10),ipa(i),side(i)
      if(atno2.lt.1.5) goto 21
      call lengtx(i,pp)
      ekin=pv(4,i)-abs(pv(5,i))
      call normal(ran)
      ekin=ekin-cfa*(1.+0.5*ran)
      if (ekin .lt. 1.0e-6) ekin=1.0e-6
      call steeq(xxh,i)
      dekin=dekin+ekin*(1.-xxh)
      ekin=ekin*xxh
      if(abs(ipa(i)).ge.7.and.abs(ipa(i)).le.9) npions=npions+1
      if(abs(ipa(i)).ge.7.and.abs(ipa(i)).le.9) ek1=ek1+ekin
      pp1=sqrt(ekin*(ekin+2.*abs(pv(5,i))))
      pv(4,i)=ekin+abs(pv(5,i))
c
      if (pp .ge. 1.0e-6) go to 8010
      call grndm(rndm,2)
      rthnve=pi*rndm(1)
      phinve=twpi*rndm(2)
      pv(1,i)=pp1*sin(rthnve)*cos(phinve)
      pv(2,i)=pp1*sin(rthnve)*sin(phinve)
      pv(3,i)=pp1*cos(rthnve)
      go to 8011
 8010 continue
      pv(1,i)=pv(1,i)*pp1/pp
      pv(2,i)=pv(2,i)*pp1/pp
      pv(3,i)=pv(3,i)*pp1/pp
 8011 continue
c
   21 continue
      if(ek1.eq.0.) goto 23
      if(npions.eq.0) goto 23
      dekin=1.+dekin/ek1
      do 22 i=1,nt
      if(abs(ipa(i)).lt.7.or.abs(ipa(i)).gt.9) goto 22
      call lengtx(i,pp)
      ekin=pv(4,i)-abs(pv(5,i))
      ekin=ekin*dekin
      if (ekin .lt. 1.0e-6) ekin=1.0e-6
      pp1=sqrt(ekin*(ekin+2.*abs(pv(5,i))))
      pv(4,i)=ekin+abs(pv(5,i))
c
      if (pp .ge. 1.0e-6) go to 8012
      call grndm(rndm,2)
      rthnve=pi*rndm(1)
      phinve=twpi*rndm(2)
      pv(1,i)=pp1*sin(rthnve)*cos(phinve)
      pv(2,i)=pp1*sin(rthnve)*sin(phinve)
      pv(3,i)=pp1*cos(rthnve)
      go to 8013
 8012 continue
      pv(1,i)=pv(1,i)*pp1/pp
      pv(2,i)=pv(2,i)*pp1/pp
      pv(3,i)=pv(3,i)*pp1/pp
 8013 continue
c
   22 continue
c**
c** add black track particles, the total number of particles produced
c** is restricted to 198, this may have influence on very high energy
c** first protons and neutrons
c**
   23 if(atno2.lt.1.5) goto 40
      call selfab(sprob)
      tex=enp(1)
      spall=targ
      if(tex.lt.0.001) goto 445
      black=(1.5+1.25*targ)*enp(1)/(enp(1)+enp(3))
      call poisso(black,nbl)
c      if (nprt(4)) write(newbcd,3003) nbl,tex
      if(ifix(targ)+nbl.gt.atno2) nbl=atno2-targ
      if(nt+nbl.gt.198) nbl=198-nt
      if(nbl.le.0) goto 445
      ekin=tex/nbl
      ekin2=0.
      call steep(xx)
      do 441 i=1,nbl
      call grndm(rndm,1)
      if(rndm(1).lt.sprob) goto 441
      if(nt.eq.198) goto 441
      if(ekin2.gt.tex) goto 443
      call grndm(rndm,1)
      ran1=rndm(1)
      call normal(ran2)
      ekin1=-ekin*log(ran1)-cfa*(1.+0.5*ran2)
      if(ekin1.lt.0.0) ekin1=-0.010*log(ran1)
      ekin1=ekin1*xx
      ekin2=ekin2+ekin1
      if(ekin2.gt.tex) ekin1=tex-(ekin2-ekin1)
      if (ekin1 .lt. 0.0) ekin1=1.0e-6
      ipa1=16
      pnrat=1.-zno2/atno2
      call grndm(rndm,3)
      if(rndm(1).gt.pnrat) ipa1=14
      nt=nt+1
      spall=spall+1.
      cost=-1.+rndm(2)*2.
      sint=sqrt(abs(1.-cost*cost))
      phi=twpi*rndm(3)
      ipa(nt)=-ipa1
      side(nt)=-4.
      pv(5,nt)=abs(rmass(ipa1))
      pv(6,nt)=rcharg(ipa1)
      pv(7,nt)=1.
      pv(4,nt)=ekin1+pv(5,nt)
      rnve=abs(pv(4,nt)**2-pv(5,nt)**2)
      pp=sqrt(rnve)
      pv(1,nt)=pp*sint*sin(phi)
      pv(2,nt)=pp*sint*cos(phi)
      pv(3,nt)=pp*cost
  441 continue
  443 if(atno2.lt.10.) goto 445
      if(ek.gt.2.0) goto 445
      ii=nt+1
      kk=0
      eka=ek
      if(eka.gt.1.) eka=eka*eka
      if(eka.lt.0.1) eka=0.1
      ika=3.6*exp((zno2**2/atno2-35.56)/6.45)/eka
      if(ika.le.0) go to 445
      do 444 i=1,nt
      ii=ii-1
      if(ipa(ii).ne.-14) goto 444
      ipa(ii)=-16
      ipa1  = 16
      pv(5,ii)=abs(rmass(ipa1))
      pv(6,ii)=rcharg(ipa1)
      kk=kk+1
      if(kk.gt.ika) goto 445
  444 continue
c**
c** then also deuterons, tritons and alphas
c**
  445 tex=enp(3)
      if(tex.lt.0.001) goto 40
      black=(1.5+1.25*targ)*enp(3)/(enp(1)+enp(3))
      call poisso(black,nbl)
      if(nt+nbl.gt.198) nbl=198-nt
      if(nbl.le.0) goto 40
      ekin=tex/nbl
      ekin2=0.
      call steep(xx)
c      if (nprt(4)) write(newbcd,3004) nbl,tex
      do 442 i=1,nbl
      call grndm(rndm,1)
      if(rndm(1).lt.sprob) goto 442
      if(nt.eq.198) goto 442
      if(ekin2.gt.tex) goto 40
      call grndm(rndm,1)
      ran1=rndm(1)
      call normal(ran2)
      ekin1=-ekin*log(ran1)-cfa*(1.+0.5*ran2)
      if(ekin1.lt.0.0) ekin1=-0.010*log(ran1)
      ekin1=ekin1*xx
      ekin2=ekin2+ekin1
      if(ekin2.gt.tex) ekin1=tex-(ekin2-ekin1)
      if (ekin1 .lt. 0.0) ekin1=1.0e-6
      call grndm(rndm,3)
      cost=-1.+rndm(1)*2.
      sint=sqrt(abs(1.-cost*cost))
      phi=twpi*rndm(2)
      ran=rndm(3)
      ipa(nt+1)=-30
      if(ran.gt.0.60) ipa(nt+1)=-31
      if(ran.gt.0.90) ipa(nt+1)=-32
      side(nt+1)=-4.
      pv(5,nt+1)=(abs(ipa(nt+1))-28)*mp
      spall=spall+pv(5,nt+1)*1.066
      if(spall.gt.atno2) goto 40
      nt=nt+1
      pv(6,nt)=1.
      if(ipa(nt).eq.-32) pv(6,nt)=2.
      pv(7,nt)=1.
      pv(4,nt)=pv(5,nt)+ekin1
      rnve=abs(pv(4,nt)**2-pv(5,nt)**2)
      pp=sqrt(rnve)
      pv(1,nt)=pp*sint*sin(phi)
      pv(2,nt)=pp*sint*cos(phi)
      pv(3,nt)=pp*cost
  442 continue
c**
c** store on event common
c**
   40 call grndm(rndm,1)
      if(rs.gt.(4.+rndm(1))) goto 42
      do 41 i=1,nt
      call lengtx(i,etb)
      if(etb.lt.p) goto 41
      etf=p
      pv(4,i)=sqrt(pv(5,i)**2+etf**2)
      etf=etf/etb
      pv(1,i)=pv(1,i)*etf
      pv(2,i)=pv(2,i)*etf
      pv(3,i)=pv(3,i)*etf
   41 continue
   42 ekin=pv(4,200)-abs(pv(5,200))
      ekin1=pv(4,199)-abs(pv(5,199))
      ekin2=0.
      call tdelay(tof1)
      call grndm(rndm,1)
      ran=rndm(1)
      tof=tof-tof1*log(ran)
      do 44 i=1,nt
      ekin2=ekin2+pv(4,i)-abs(pv(5,i))
      if(pv(7,i).lt.0.) pv(5,i)=-pv(5,i)
      pv(7,i)=tof
      pv(8,i)=abs(ipa(i))
      pv(9,i)=0.
   44 pv(10,i)=0.
c      if (nprt(4)) write(newbcd,2006) nt,ekin,enp(1),enp(3),ekin1,ekin2
      intct=intct+1.
      nmode=4
      if(spall.lt.0.5.and.atno2.gt.1.5) nmode=14
      call setcur(nt)
      ntk=ntk+1
      if(nt.eq.1) go to 9999
      do 50 ii=2,nt
      i=ii-1
      if(ntot.lt.nsize/12) goto 43
      go to 9999
   43 call settrk(i)
   50 continue
c
 2002 format(' *genxpt* production of final state kinematic after ',i3,
     $ ' trials.  kinetic energies ',2f6.2,' out of ',2f6.2)
 2000 format(' *genxpt* cms parameters of final state particles,',
     $ ' energies in initial and final state ',2f6.2)
 2001 format(' *genxpt* track',2x,i3,2x,10f8.3,2x,i3,2x,f4.0)
 2003 format(' *genxpt* teta,ekin0,ekin1,ekin ',4f10.4)
 2006 format(' *genxpt* comp.',1x,i5,1x,5f7.2)
 3001 format(' *genxpt* nuclear excitation',i5,
     $ ' particles produced in addition  to ',i5,' normal particles')
 3002 format(' *genxpt* available energies ',2f10.4,
     $ ' for ',2i3,' particles in beam/target fragm. region',
     $ ' with ipa/side array '/
     $ 1h ,5x,10(i3,2x,f3.0,4x))
 3003 format(' *genxpt* ',i3,' black track particles produced',
     $ ' with total kinetic energy of ',f8.3,' gev')
 3004 format(' *genxpt* ',i5,' heavy fragments produced',
     $ ' with total energy of',f8.4,' gev')
c
 9999 continue
c
      return
      end
