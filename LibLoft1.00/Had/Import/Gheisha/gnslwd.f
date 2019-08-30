*cmz :  3.14/16 13/03/89  14.42.12  by  nick van eijndhoven (cern)
*-- author :
      subroutine gnslwd(nucflg,int,nfl,teklow,stop)
c
c *** neutron tracking routine for energies below the cut-off. ***
c *** take only elastic scattering, neutron capture            ***
c *** and nuclear fission.                                     ***
c *** nve 11-may-1988 cern geneva ***
c
c called by : gheish
c origin : h.fesefeldt (routine nsldow 20-oct-1987)
c
      integer nmec,lmec,namec,nstep ,maxnst,ignext,inwvol,istop,maxmec
     + ,igauto,iekbin,ilosl, imull,ingoto,nldown,nlevin,nlvsav,istory
      real  vect,getot,gekin,vout,destep,destel,safety,sleng ,step
     + ,snext,sfield,tofg  ,gekrat,upwght
      parameter (maxmec=30)
      common/gctrak/vect(7),getot,gekin,vout(7),nmec,lmec(maxmec)
     + ,namec(maxmec),nstep ,maxnst,destep,destel,safety,sleng
     + ,step  ,snext ,sfield,tofg  ,gekrat,upwght,ignext,inwvol
     + ,istop ,igauto,iekbin, ilosl, imull,ingoto,nldown,nlevin
     + ,nlvsav,istory
c
c --- gheisha commons ---
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
      common /vecuty/ pv(10,200)
c
      common/prntfl/inbcd,newbcd,inbin,newbin,npevt,nevtp,lprt,nprt(10)
                    logical lprt,nprt
c
      dimension rndm(2)
c
c --- flags to indicate the nucrec action ---
c nucflg = 0 ==> no action by nucrec
c          1 ==> action by nucrec ==> special treatment in gheish
      nopt=0
      nucflg=0
c
c --- in order to avoid troubles caused by arithmetic incertainties, ---
c --- recalculate some quantities. take kinetic energy ek as most ---
c --- relevant quantity. ---
c
c --- very low kinetic energy ==> neutron capture ---
      if (ek .lt. 1.e-9) go to 22
c
      en=ek+abs(amas)
      p=sqrt(abs(en*en-amas*amas))
      pu=sqrt(px**2+py**2+pz**2)
      if (pu .ge. 1.e-9) go to 7
c
      px=0.0
      py=0.0
      pz=0.0
      go to 22
c
 7    continue
      px=px/pu
      py=py/pu
      pz=pz/pu
c
c --- select process according to "int" ---
      go to (23,23,21,22), int
c
c *** nuclear fission ***
 21   continue
      stop=1
      tkin=fissio(ek)
      go to 9999
c
c *** neutron capture ***
 22   continue
      stop=1
      call captur(nopt)
      go to 9999
c
c *** elastic and inelastic scattering ***
 23   continue
      pv(1,200)=p*px
      pv(2,200)=p*py
      pv(3,200)=p*pz
      pv(4,200)=en
      pv(5,200)=amas
      pv(6,200)=nch
      pv(7,200)=tof
      pv(8,200)=ipart
      pv(9,200)=0.0
      pv(10,200)=userw
c
c --- special treatment for inelastic scattering in heavy media ---
      if ((int .eq. 2) .and. (atno2 .ge. 1.5)) go to 29
c
c *** elastic scattering ***
 30   continue
c
c      if (nprt(9)) print 1000
 1000 format(' *gnslwd* elastic scattering')
c
      do 24 j=4,9
      pv(j,1)=pv(j,200)
 24   continue
      pv(10,1)=0.0
c
c --- very simple simulation of scattering angle and energy ---
c --- nonrelativistic approximation with isotropic angular ---
c --- distribution in the cms system ---
      call grndm(rndm,2)
      ran=rndm(1)
      cost1=-1.0+2.0*ran
      eka=1.0+2.0*cost1*atno2+atno2**2
      cost=(atno2*cost1+1.0)/sqrt(eka)
      if (cost .lt. -1.0) cost=-1.0
      if (cost .gt. 1.0) cost=1.0
      eka=eka/(1.0+atno2)**2
      ek=ek*eka
      en=ek+abs(amas)
      p=sqrt(abs(en*en-amas*amas))
      sint=sqrt(abs(1.0-cost*cost))
      phi=rndm(2)*twpi
      pv(1,2)=sint*sin(phi)
      pv(2,2)=sint*cos(phi)
      pv(3,2)=cost
      call defs1(2,200,2)
      pu=sqrt(pv(1,2)**2+pv(2,2)**2+pv(3,2)**2)
      px=pv(1,2)/pu
      py=pv(2,2)/pu
      pz=pv(3,2)/pu
      pv(1,1)=px*p
      pv(2,1)=py*p
      pv(3,1)=pz*p
      pv(4,1)=en
c
c --- store backscattered particle for atno < 4.5 ---
      if (atno2 .gt. 4.5) go to 27
c
c      if (nprt(9)) print 1001,atno2
 1001 format(' *gnslwd* backscattered particle stored for atno ',g12.5)
c
      pv(1,2)=pv(1,200)-pv(1,1)
      pv(2,2)=pv(2,200)-pv(2,1)
      pv(3,2)=pv(3,200)-pv(3,1)
      call lengtx(2,pp)
      pv(9,2)=0.0
      pv(10,2)=0.0
      pv(7,2)=tof
c
      if (atno2 .gt. 3.5) go to 274
      if (atno2 .gt. 2.5) go to 273
      if (atno2 .gt. 1.5) go to 272
c
 271  continue
      pv(5,2)=rmass(14)
      pv(4,2)=sqrt(pp*pp+pv(5,2)*pv(5,2))
      pv(6,2)=rcharg(14)
      pv(8,2)=14.0
      go to 275
c
 272  continue
      pv(5,2)=rmass(30)
      pv(4,2)=sqrt(pp*pp+pv(5,2)*pv(5,2))
      pv(6,2)=rcharg(30)
      pv(8,2)=30.0
      go to 275
c
 273  continue
      pv(5,2)=rmass(31)
      pv(4,2)=sqrt(pp*pp+pv(5,2)*pv(5,2))
      pv(6,2)=rcharg(31)
      pv(8,2)=31.0
      go to 275
c
 274  continue
      pv(5,2)=rmass(32)
      pv(4,2)=sqrt(pp*pp+pv(5,2)*pv(5,2))
      pv(6,2)=rcharg(32)
      pv(8,2)=32.0
c
 275  continue
      intct=intct+1.0
      call setcur(1)
      ntk=ntk+1
      call settrk(2)
      go to 9999
c
c --- put quantities in common /result/ ---
 27   continue
      if (pv(10,1) .ne. 0.0) userw=pv(10,1)
      sinl=pz
      cosl=sqrt(abs(1.0-sinl*sinl))
      if (abs(cosl) .lt. 1.e-10) go to 28
c
      sinp=py/cosl
      cosp=px/cosl
      go to 9999
c
 28   continue
      call grndm(rndm,1)
      phi=rndm(1)*twpi
      sinp=sin(phi)
      cosp=cos(phi)
      go to 9999
c
c *** inelastic scattering on heavy nuclei ***
 29   continue
c
c      if (nprt(9)) print 1002
 1002 format(' *gnslwd* inelastic scattering on heavy nucleus')
c
c --- decide between spallation or simple nuclear reaction ---
      call grndm(rndm,1)
      test1=rndm(1)
      test2=4.5*(ek-0.01)
      if (test1 .gt. test2) go to 40
c
c *** spallation ***
c
c      if (nprt(9)) print 1003
 1003 format(' *gnslwd* spallation')
c
      pv(1,200)=p*px
      pv(2,200)=p*py
      pv(3,200)=p*pz
      pv(4,200)=en
      pv(5,200)=amas
      pv(6,200)=nch
      pv(7,200)=tof
      pv(8,200)=ipart
      pv(9,200)=0.0
      pv(10,200)=userw
c
c --- fermi-motion and evaporation ---
      tkin=cinema(ek)
      enp(5)=ek+tkin
c --- check for lowerbound of ekin in cross-section tables ---
      if (enp(5) .le. teklow) enp(5)=teklow
      enp(6)=enp(5)+abs(amas)
      enp(7)=enp(6)*enp(6)-amasq
      enp(7)=sqrt(enp(7))
      tkin=fermi(enp(5))
      enp(5)=enp(5)+tkin
c --- check for lowerbound of ekin in cross-section tables ---
      if (enp(5) .le. teklow) enp(5)=teklow
      enp(6)=enp(5)+abs(amas)
      enp(7)=enp(6)*enp(6)-amasq
      enp(7)=sqrt(enp(7))
      tkin=exnu(enp(5))
      enp(5)=enp(5)-tkin
c --- check for lowerbound of ekin in cross-section tables ---
      if (enp(5) .le. teklow) enp(5)=teklow
      enp(6)=enp(5)+abs(amas)
      enp(7)=enp(6)*enp(6)-amasq
      enp(7)=sqrt(enp(7))
c
c --- neutron cascade ---
      k=2
      call vzero(ipa(1),100)
      call casn(k,int,nfl)
      go to 9999
c
 40   continue
c      if (nprt(9)) print 1004
 1004 format(' *gnslwd* nuclear reaction')
      call nucrec(nopt,1)
      if (nopt .ne. 0) nucflg=1
      if (nopt .eq. 0) go to 30
c
 9999 continue
      return
      end
