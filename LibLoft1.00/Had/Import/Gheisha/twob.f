*cmz :          10/01/91  18.53.21  by  federico carminati
*cmz :  3.14/16 01/11/90  16.16.00  by  rene brun
*-- author :
      subroutine twob(ippp,nfl,avern)
c
c *** generation of momenta for elast. and quasi elast. 2 body react. ***
c *** nve 04-may-1988 cern geneva ***
c
c origin : h.fesefeldt 15-sep-1987
c
c the simple formula ds/d|t| = s0* exp(-b*|t|) is used.
c the b values are parametrizations from experimental data .
c not available values are taken from those of similar reactions
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
      dimension rndm(3)
c
c     data cb/3./
      data cb/0.01/
c
c --- statement functions ---
      bpp(x)=4.225+1.795*log(x)
c
c**
c**  for diffraction scattering on heavy nuclei use better routine
c**  "coscat"
c
      tarmas=rmass(14)
      if (nfl .eq. 2) tarmas=rmass(16)
      enp(8)=rmass(ippp)**2+tarmas**2+2.0*tarmas*enp(6)
      enp(9)=sqrt(enp(8))
      ek=enp(5)
      en=enp(6)
      p=enp(7)
      s=enp(8)
      rs=enp(9)
      cfa=0.025*((atno2-1.)/120.)*exp(-(atno2-1.)/120.)
      if(atno2.lt.1.5) goto 500
      ipa1=abs(ipa(1))
      ipa2=abs(ipa(2))
      rmc=rmass(ipa1)
      rmd=rmass(ipa2)
      rchc=rcharg(ipa1)
      rchd=rcharg(ipa2)
      if(abs(rmc-amas).gt.0.001) goto 500
      rmnve=rmass(14)
      if (nfl .eq. 2) rmnve=rmass(16)
      if(abs(rmd-rmnve).gt.0.001) goto 500
      if(abs(rchc-nch).gt.0.5) goto 500
      if(nfl.eq.1.and.rchd.lt.0.5) goto 500
      if(nfl.eq.2.and.abs(rchd).gt.0.5) goto 500
      if(enp(1).gt.0.0001.or.enp(3).gt.0.0001) goto 500
      call coscat
      go to 9999
c**
c**  set effective 4-momentum of initial particle
c**
  500 pv(1,199)=p*px
      pv(2,199)=p*py
      pv(3,199)=p*pz
      pv(4,199)=en
      pv(5,199)=amas
      pv(6,199)=nch
      pv(7,199)=tof
      pv(8,199)=ipart
      pv(9,199)=0.
      pv(10,199)=userw
      ier(47)=ier(47)+1
c      if (nprt(4)) write(newbcd,4001) (pv(j,199),j=1,10),ipa(1),ipa(2)
      do 2 j=1,6
    2 pv(j,1)=pv(j,199)
      pv(7,1)=1.
      if(pv(5,1).lt.0.) pv(7,1)=-1.
      pv(5,1)=abs(pv(5,1))
      nt=1
c**
c** two-body scattering possible?? if not, continue with original
c** particle, but spend the nuclear evaporation energy
c**
      if(p.lt.0.1) goto 200
      if(rs.lt.0.01) goto 200
c**
c** calculate slope b for elastic scattering on proton/neutron
c**
      b=bpp(p)
      if(b.lt.cb) b=cb
      if(abs(ipa(2)).gt.13) goto 9
      ipa(2)=14
      call grndm(rndm,1)
      if(rndm(1).lt.0.5) ipa(2)=16
c**
c** set masses and momenta for final state particles
c**
    9 rmc=rmass(abs(ipa(1)))
      rmd=rmass(abs(ipa(2)))
      pv(6,1)=rcharg(abs(ipa(1)))
      pv(6,2)=rcharg(abs(ipa(2)))
      pf=(s+rmd*rmd-rmc*rmc)**2 - 4*s*rmd*rmd
      if(nprt(4)) write(newbcd,4002) rmc,rmd,pv(6,1),pv(6,2),rs,s,pf
      if(pf.lt.0.001) go to 9999
      pf=sqrt(pf)/(2.*rs)
c**
c** set beam and target in cms
c**
      pv(1,3)=0.
      pv(2,3)=0.
      pv(3,3)=p
      pv(5,3)=abs(amas)
      pv(4,3)=sqrt(p*p+amas*amas)
      pv(1,4)=0.
      pv(2,4)=0.
      pv(3,4)=0.
      rmnve=rmass(14)
      if (nfl .eq. 2) rmnve=rmass(16)
      pv(4,4)=rmnve
      pv(5,4)=rmnve
c**
c** transform into cms.
c**
      call add(3,4,10)
      call lor(3,10,3)
      call lor(4,10,4)
c**
c** set final state masses and energies in cms
c**
      pv(5,1)=abs(rmc)
      pv(5,2)=abs(rmd)
      pv(7,1)=1.
      pv(7,2)=1.
      if(rmc.lt.0.) pv(7,1)=-1.
      if(rmd.lt.0.) pv(7,2)=-1.
      pv(4,1)=sqrt(pf*pf+pv(5,1)*pv(5,1))
      pv(4,2)=sqrt(pf*pf+pv(5,2)*pv(5,2))
c**
c** set |t| and |tmin|
c**
      call grndm(rndm,2)
      call lengtx(3,pin)
      btrang=b*4.*pin*pf
      tacmin=(pv(4,3)-pv(4,1))**2-(pin-pf)**2
c**
c** simply a protection against exponent overflow 1.e20 is big enough
c**
      exindt=-1.
      if(btrang.lt.46) exindt=exindt+exp(-btrang)
      tdn=log(1.+rndm(1)*exindt)/btrang
c**
c** caculate (sin(teta/2.)**2 and cos(teta), set azimuth angle phi
c**
      ctet=1.+2.*tdn
      if(abs(ctet).gt.1.) ctet=sign(1.,ctet)
      stet=sqrt((1.-ctet)*(1.+ctet))
      phi=rndm(2)*twpi
c**
c** calculate final state momenta in cms
c**
      pv(1,1)=pf*stet*sin(phi)
      pv(2,1)=pf*stet*cos(phi)
      pv(3,1)=pf*ctet
      pv(1,2)=-pv(1,1)
      pv(2,2)=-pv(2,1)
      pv(3,2)=-pv(3,1)
c**
c** transform into lab
c**
      do 11 i=1,2
      call lor(i,4,i)
      call defs1(i,199,i)
      if(atno2.lt.1.5) goto 11
      call lengtx(i,pp)
      if(pp.lt.0.001) goto 11
      ekin=pv(4,i)-abs(pv(5,i))
      call normal(ran)
      ekin=ekin-cfa*(1.+0.5*ran)
      if(ekin.lt.0.0001) ekin=0.0001
      pp1=sqrt(ekin*(ekin+2.*abs(pv(5,i))))
      pv(4,i)=ekin+abs(pv(5,i))
      pv(1,i)=pv(1,i)*pp1/pp
      pv(2,i)=pv(2,i)*pp1/pp
      pv(3,i)=pv(3,i)*pp1/pp
   11 continue
      nt=2
c**
c** add black track particles .
c** here the procedure is somewhat different as in 'twoclu' and 'genxpt'
c** the reason is, that we have to simulate also the nuclear reactions
c** at low energies like a(h,p)b, a(h,p p)b, a(h,n)b e.t.c.
c**
  200 if(enp(1).le.0.0001.and.enp(3).le.0.0001) goto 40
      spall=0.
      tex=enp(1)
      if(tex.lt.0.0001) goto 445
      black=tex/0.02
      call poisso(black,nbl)
      if(nbl.gt.atno2) nbl=atno2
      if(enp(1).gt.0.0001.and.nbl.le.0) nbl=1
c      if (nprt(4)) write(newbcd,3003) nbl,tex
      if(nt+nbl.gt.198) nbl=198-nt
      if(nbl.le.0) goto 445
      ekin=tex/nbl
      ekin2=0.
      call steep(xx)
      do 441 i=1,nbl
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
      if(ekin1.lt.0.) ekin1=0.0001
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
      pv(5,nt)=abs(rmass(ipa1))
      pv(6,nt)=rcharg(ipa1)
      pv(7,nt)=2.
      pv(4,nt)=ekin1+pv(5,nt)
      pp=sqrt(abs(pv(4,nt)**2-pv(5,nt)**2))
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
  445 tex=enp(3)
      if(tex.lt.0.0001) goto 40
      nbl=ifix(2.*log(atno2))
      if(nt+nbl.gt.198) nbl=198-nt
      if(nbl.le.0) goto 40
      ekin=tex/nbl
      ekin2=0.
      call steep(xx)
c      if (nprt(4)) write(newbcd,3004) nbl,tex
      do 442 i=1,nbl
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
      if(ekin1.lt.0.) ekin1=0.0001
      call grndm(rndm,3)
      cost=-1.+rndm(1)*2.
      sint=sqrt(abs(1.-cost*cost))
      phi=twpi*rndm(2)
      ran=rndm(3)
      ipa(nt+1)=-30
      if(ran.gt.0.60) ipa(nt+1)=-31
      if(ran.gt.0.90) ipa(nt+1)=-32
      inve=abs(ipa(nt+1))
      pv(5,nt+1)=rmass(inve)
      spall=spall+pv(5,nt+1)*1.066
      if(spall.gt.atno2) goto 40
      nt=nt+1
      pv(6,nt)=rcharg(inve)
      pv(7,nt)=2.
      pv(4,nt)=pv(5,nt)+ekin1
      pp=sqrt(abs(pv(4,nt)**2-pv(5,nt)**2))
      pv(1,nt)=pp*sint*sin(phi)
      pv(2,nt)=pp*sint*cos(phi)
      pv(3,nt)=pp*cost
  442 continue
c**
c** store on event common
c**
   40 ekin=pv(4,200)-abs(pv(5,200))
      ekin1=pv(4,199)-abs(pv(5,199))
      ekin2=0.
      call tdelay(tof1)
      call grndm(rndm,1)
      ran=rndm(1)
      tof=tof-tof1*log(ran)
      do 1 i=1,nt
      ekin2=ekin2+pv(4,i)-abs(pv(5,i))
      if(pv(7,i).lt.0.) pv(5,i)=-pv(5,i)
      pv(7,i)=tof
      pv(8,i)=abs(ipa(i))
      pv(9,i)=0.
    1 pv(10,i)=0.
c      if (nprt(4)) write(newbcd,1003) nt,ekin,ekin1,ekin2
      intct=intct+1.
      nmode=2
      if(spall.lt.0.5.and.atno2.gt.1.5) nmode=14
      call setcur(nt)
      ntk=ntk+1
      if(nt.eq.1) go to 9999
      do 50 ii=2,nt
      i=ii-1
      if(ntot.lt.nsize/12) goto 43
      return
   43 call settrk(i)
   50 continue
c
 1001 format(' *twob* elastic scattering, mass indices ',
     $ 2i3,' --> ',2i3,' lab momentum ',f8.4,
     $ ' masses ',2f8.4,' --> ',2f8.4/
     $ ' final state momentum ',f8.4,' t,tacmin,ctet ',3f8.4/
     $ ' center of mass system four vectors ')
 1002 format(' *twob* ',5f10.4,10x,5f10.4/1h ,7x,5f10.4,10x,5f10.4/
     $ ' lab system final state four vectors')
 1003 format(' *twob* comparison',2x,i5,2x,3f10.4)
 4001 format(' *twob* ',10f10.4,2x,2i3)
 4002 format(' *twob* ',7f10.4)
 3003 format(' *twob* ',i3,' black track particles produced',
     $ ' with total kinetic energy of ',f8.3,' gev')
 3004 format(' *twob* ',i5,' heavy fragments produced',
     $ ' with total energy of',f8.4,' gev')
c
 9999 continue
c
      return
      end
