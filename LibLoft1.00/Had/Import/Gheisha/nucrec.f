*cmz :  3.14/16 29/10/90  15.27.30  by  rene brun
*-- author :
      subroutine nucrec(nopt,irec)
c
c *** nuclear reaction kinematics at low energies ***
c *** nve 18-may-1988 cern geneva ***
c
c called by : gheish, gnslwd
c origin    : h.fesefeldt (12-feb-1987)
c
c nopt=1   n m(a,z) --> g (g) m(a+1,z  )    neutron capture
c nopt=2   n m(a,z) --> n (g) m(a  ,z  )    inelastic neutron scatt.
c nopt=3   n m(a,z) --> p (g) m(a  ,z-1)
c nopt=4   n m(a,z) --> d (g) m(a-1,z-1)
c nopt=5   n m(a,z) --> t (g) m(a-2,z-1)
c nopt=6   n m(a,z) --> alp.  m(a-3,z-2)
c nopt=7   n m(a,z) --> n n   m(a-1,z  )
c nopt=8   n m(a,z) --> n p   m(a-1,z-1)
c nopt=9   n m(a,z) --> p p   m(a-1,z-2)
c nopt=11  p m(a,z) --> g (g) m(a+1,z+1)    proton capture
c nopt=12  p m(a,z) --> n (g) m(a  ,z  )    inelastic proton scatt.
c nopt=13  p m(a,z) --> p (g) m(a  ,z+1)
c nopt=14  p m(a,z) --> d (g) m(a-1,z  )
c nopt=15  p m(a,z) --> t (g) m(a-2,z  )
c nopt=16  p m(a,z) --> alp.  m(a-3,z-1)
c nopt=17  p m(a,z) --> n n   m(a-1,z+1)
c nopt=18  p m(a,z) --> n p   m(a-1,z  )
c nopt=19  p m(a,z) --> p p   m(a-1,z-1)
c similar for d,t,alpha scattering on nuclei
c
c note : double precision calculations are vital for these low
c        energy processes
c        therefore the vars of /genio/ are declared double precision
c        also a double precision version of the phase space package
c        "phpnuc" has been introduced
c *** hmf 29-aug-1989 rwth aachen ***
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
      common/nucin /tecm,amass(18),npg,kgenev
      common/nucout/pcm(5,18),wgt
      double precision tecm,amass,pcm,wgt
c
c
c
      dimension qval(10),tch(10)
      dimension rndm(2)
c
c** program returns with nopt=0, if inelastic scattering energetically
c** not possible, or if wrong particles enter this routine: only for
c** protons,neutrons, deuterium, tritium and alphas.
c** if ek > 100 mev, this routine is certainly not adequate.
c
      nopt=0
      if (irec .eq. 0) go to 9999
c
c      if (nprt(9) .and. (ek .gt. 0.1)) print 9000,ek,ipart
 9000 format(' *nucrec* energy too high ek = ',g12.5,' gev ',
     $ ' kpart = ',i3)
      if (ek .gt. 0.1) go to 9999
c
c%%%      if(ipart.eq.16) goto 2
c%%%      if(ipart.eq.14) goto 3
c%%%      if(ipart.eq.30) goto 4
c%%%      if(ipart.eq.31) goto 5
c%%%      if(ipart.eq.32) goto 6
c%%%      go to 9999
c%%%    2 amas = atomas(1.,0.)
c%%%      goto 8
c%%%    3 amas = atomas(1.,1.)
c%%%      goto 8
c%%%    4 amas = atomas(2.,1.)
c%%%      goto 8
c%%%    5 amas = atomas(3.,1.)
c%%%      goto 8
c%%%    6 amas = atomas(4.,2.)
c
      if (ipart .eq. 16) go to 8
      if (ipart .eq. 14) go to 8
      if (ipart .eq. 30) go to 8
      if (ipart .eq. 31) go to 8
      if( ipart .eq. 32) go to 8
      go to 9999
c** set beam particle, take ek as fundamental quantity
c** due to the difficult kinematic, all masses have to be assigned
c** the best measured values.
 8    continue
      call vzero(qval,10)
      call vzero(tch ,10)
c --- get mass which matches geant ---
      amas=rmass(ipart)
      en=ek+abs(amas)
      p =sqrt(abs(en*en-amas*amas))
      pp=sqrt(px*px+py*py+pz*pz)
      if (pp .gt. 1.0e-6) go to 8000
      call grndm(rndm,2)
      phinve=twpi*rndm(1)
      cost=-1.+2.*rndm(2)
      if (cost .le. -1.) cost=-1.
      if (cost .ge.  1.) cost= 1.
      rthnve=acos(cost)
      px=sin(rthnve)*cos(phinve)
      py=sin(rthnve)*sin(phinve)
      pz=cos(rthnve)
      pp=1.
 8000 continue
      px=px/pp
      py=py/pp
      pz=pz/pp
      call vzero(pv,100)
      pv(1,1) =px*p
      pv(2,1) =py*p
      pv(3,1) =pz*p
      pv(4,1) =en
      pv(5,1) =amas
      pv(6,1) =nch
      pv(7,1) =tof
      pv(8,1) =ipart
      pv(9,1) =0.
      pv(10,1)=userw
      pv(1,2) =0.
      pv(2,2) =0.
      pv(3,2) =0.
      pv(4,2) =0.
      pv(5,2) =atomas(atno2,zno2)
      pv(6,2) =zno2
      pv(7,2) =tof
      pv(8,2) =0.
      pv(9,2) =0.
      pv(10,2)=0.
c** calculate q-value of reactions
      if(ipart.eq.16) goto 20
      if(ipart.eq.14) goto 30
      if(ipart.eq.30) goto 40
      if(ipart.eq.31) goto 50
      if(ipart.eq.32) goto 60
   20 pv(5,11)=atomas(atno2+1.,zno2   )
      pv(6,11)=zno2
      pv(5,21)=0.
      pv(6,21)=0.
      pv(8,21)=1.
      pv(5,31)=0.
      pv(6,31)=0.
      pv(8,31)=1.
c
      pv(5,12)=pv(5,2)
      pv(6,12)=pv(6,2)
      pv(5,22)=rmass(16)
      pv(6,22)=0.
      pv(8,22)=16.
      pv(5,32)=0.
      pv(6,32)=0.
      pv(8,32)=1.
c
      pv(5,13)=atomas(atno2   ,zno2-1.)
      pv(6,13)=zno2-1.
      pv(5,23)=rmass(14)
      pv(6,23)=1.
      pv(8,23)=14.
      pv(5,33)=0.
      pv(6,33)=0.
      pv(8,33)=1.
c
      pv(5,14)=atomas(atno2-1.,zno2-1.)
      pv(6,14)=zno2-1.
      pv(5,24)=rmass(30)
      pv(6,24)=1.
      pv(8,24)=30.
      pv(5,34)=0.
      pv(6,34)=0.
      pv(8,34)=1.
c
      pv(5,15)=atomas(atno2-2.,zno2-1.)
      pv(6,15)=zno2-1.
      pv(5,25)=rmass(31)
      pv(6,25)=1.
      pv(8,25)=31.
      pv(5,35)=0.
      pv(6,35)=0.
      pv(8,35)=1.
c
      pv(5,16)=atomas(atno2-3.,zno2-2.)
      pv(6,16)=zno2-2.
      pv(5,26)=rmass(32)
      pv(6,26)=2.
      pv(8,26)=32.
      pv(5,36)=0.
      pv(6,36)=0.
      pv(8,36)=1.
c
      pv(5,17)=atomas(atno2-1.,zno2   )
      pv(6,17)=zno2
      pv(5,27)=pv(5,22)
      pv(6,27)=0.
      pv(8,27)=16.
      pv(5,37)=pv(5,22)
      pv(6,37)=0.
      pv(8,37)=16.
c
      pv(5,18)=pv(5,14)
      pv(6,18)=pv(6,14)
      pv(5,28)=pv(5,22)
      pv(6,28)=0.
      pv(8,28)=16.
      pv(5,38)=pv(5,23)
      pv(6,38)=1.
      pv(8,38)=14.
c
      pv(5,19)=atomas(atno2-1.,zno2-2.)
      pv(6,19)=zno2-2.
      pv(5,29)=pv(5,23)
      pv(6,29)=1.
      pv(8,29)=14.
      pv(5,39)=pv(5,23)
      pv(6,39)=1.
      pv(8,39)=14.
c
      goto 70
   30 pv(5,11)=atomas(atno2+1.,zno2+1.)
      pv(6,11)=zno2+1.
      pv(5,21)=0.
      pv(6,21)=0.
      pv(8,21)=1.
      pv(5,31)=0.
      pv(6,31)=0.
      pv(8,31)=1.
c
      pv(5,12)=atomas(atno2   ,zno2+1.)
      pv(6,12)=zno2+1.
      pv(5,22)=rmass(16)
      pv(6,22)=0.
      pv(8,22)=16.
      pv(5,32)=0.
      pv(6,32)=0.
      pv(8,32)=1.
c
      pv(5,13)=pv(5,2)
      pv(6,13)=pv(6,2)
      pv(5,23)=rmass(14)
      pv(6,23)=1.
      pv(8,23)=14.
      pv(5,33)=0.
      pv(6,33)=0.
      pv(8,33)=1.
c
      pv(5,14)=atomas(atno2-1.,zno2   )
      pv(6,14)=zno2
      pv(5,24)=rmass(30)
      pv(6,24)=1.
      pv(8,24)=30.
      pv(5,34)=0.
      pv(6,34)=0.
      pv(8,34)=1.
c
      pv(5,15)=atomas(atno2-2.,zno2   )
      pv(6,15)=zno2
      pv(5,25)=rmass(31)
      pv(6,25)=1.
      pv(8,25)=31.
      pv(5,35)=0.
      pv(6,35)=0.
      pv(8,35)=1.
c
      pv(5,16)=atomas(atno2-3.,zno2-1.)
      pv(6,16)=zno2-1.
      pv(5,26)=rmass(32)
      pv(6,26)=2.
      pv(8,26)=32.
      pv(5,36)=0.
      pv(6,36)=0.
      pv(8,36)=1.
c
      pv(5,17)=atomas(atno2-1.,zno2+1.)
      pv(6,17)=zno2+1.
      pv(5,27)=pv(5,22)
      pv(6,27)=0.
      pv(8,27)=16.
      pv(5,37)=pv(5,22)
      pv(6,37)=0.
      pv(8,37)=16.
c
      pv(5,18)=pv(5,14)
      pv(6,18)=pv(6,14)
      pv(5,28)=pv(5,22)
      pv(6,28)=0.
      pv(8,28)=16.
      pv(5,38)=pv(5,23)
      pv(6,38)=1.
      pv(8,38)=14.
c
      pv(5,19)=atomas(atno2-1.,zno2-1.)
      pv(6,19)=zno2-1.
      pv(5,29)=pv(5,23)
      pv(6,29)=1.
      pv(8,29)=14.
      pv(5,39)=pv(5,23)
      pv(6,39)=1.
      pv(8,39)=14.
c
      nopt=10
      goto 70
   40 pv(5,11)=atomas(atno2+2.,zno2+1.)
      pv(6,11)=zno2+1.
      pv(5,21)=0.
      pv(6,21)=0.
      pv(8,21)=1.
      pv(5,31)=0.
      pv(6,31)=0.
      pv(8,31)=1.
c
      pv(5,12)=atomas(atno2+1.,zno2+1.)
      pv(6,12)=zno2+1.
      pv(5,22)=rmass(16)
      pv(6,22)=0.
      pv(8,22)=16.
      pv(5,32)=0.
      pv(6,32)=0.
      pv(8,32)=1.
c
      pv(5,13)=atomas(atno2+1.,zno2   )
      pv(6,13)=zno2
      pv(5,23)=rmass(14)
      pv(6,23)=1.
      pv(8,23)=14.
      pv(5,33)=0.
      pv(6,33)=0.
      pv(8,33)=1.
c
      pv(5,14)=pv(5,2)
      pv(6,14)=pv(6,2)
      pv(5,24)=rmass(30)
      pv(6,24)=1.
      pv(8,24)=30.
      pv(5,34)=0.
      pv(6,34)=0.
      pv(8,34)=1.
c
      pv(5,15)=atomas(atno2-1.,zno2   )
      pv(6,15)=zno2
      pv(5,25)=rmass(31)
      pv(6,25)=1.
      pv(8,25)=31.
      pv(5,35)=0.
      pv(6,35)=0.
      pv(8,35)=1.
c
      pv(5,16)=atomas(atno2-2.,zno2-1.)
      pv(6,16)=zno2-1.
      pv(5,26)=rmass(32)
      pv(6,26)=2.
      pv(8,26)=32.
      pv(5,36)=0.
      pv(6,36)=0.
      pv(8,36)=1.
c
      pv(5,17)=atomas(atno2   ,zno2+1.)
      pv(6,17)=zno2+1.
      pv(5,27)=pv(5,22)
      pv(6,27)=0.
      pv(8,27)=16.
      pv(5,37)=pv(5,22)
      pv(6,37)=0.
      pv(8,37)=16.
c
      pv(5,18)=pv(5,14)
      pv(6,18)=pv(6,14)
      pv(5,28)=pv(5,22)
      pv(6,28)=0.
      pv(8,28)=16.
      pv(5,38)=pv(5,23)
      pv(6,38)=1.
      pv(8,38)=14.
c
      pv(5,19)=atomas(atno2   ,zno2-1.)
      pv(6,19)=zno2-1.
      pv(5,29)=pv(5,23)
      pv(6,29)=1.
      pv(8,29)=14.
      pv(5,39)=pv(5,23)
      pv(6,39)=1.
      pv(8,39)=14.
c
      nopt=20
      goto 70
   50 pv(5,11)=atomas(atno2+3.,zno2+1.)
      pv(6,11)=zno2+1.
      pv(5,21)=0.
      pv(6,21)=0.
      pv(8,21)=1.
      pv(5,31)=0.
      pv(6,31)=0.
      pv(8,31)=1.
c
      pv(5,12)=atomas(atno2+2.,zno2+1.)
      pv(6,12)=zno2+1.
      pv(5,22)=rmass(16)
      pv(6,22)=0.
      pv(8,22)=16.
      pv(5,32)=0.
      pv(6,32)=0.
      pv(8,32)=1.
c
      pv(5,13)=atomas(atno2+2.,zno2   )
      pv(6,13)=zno2
      pv(5,23)=rmass(14)
      pv(6,23)=1.
      pv(8,23)=14.
      pv(5,33)=0.
      pv(6,33)=0.
      pv(8,33)=1.
c
      pv(5,14)=atomas(atno2+1.,zno2   )
      pv(6,14)=zno2
      pv(5,24)=rmass(30)
      pv(6,24)=1.
      pv(8,24)=30.
      pv(5,34)=0.
      pv(6,34)=0.
      pv(8,34)=1.
c
      pv(5,15)=pv(5,2)
      pv(6,15)=pv(6,2)
      pv(5,25)=rmass(31)
      pv(6,25)=1.
      pv(8,25)=31.
      pv(5,35)=0.
      pv(6,35)=0.
      pv(8,35)=1.
c
      pv(5,16)=atomas(atno2-1.,zno2-1.)
      pv(6,16)=zno2-1.
      pv(5,26)=rmass(32)
      pv(6,26)=2.
      pv(8,26)=32.
      pv(5,36)=0.
      pv(6,36)=0.
      pv(8,36)=1.
c
      pv(5,17)=atomas(atno2+1.,zno2+1.)
      pv(6,17)=zno2+1.
      pv(5,27)=pv(5,22)
      pv(6,27)=0.
      pv(8,27)=16.
      pv(5,37)=pv(5,22)
      pv(6,37)=0.
      pv(8,37)=16.
c
      pv(5,18)=pv(5,14)
      pv(6,18)=pv(6,14)
      pv(5,28)=pv(5,22)
      pv(6,28)=0.
      pv(8,28)=16.
      pv(5,38)=pv(5,23)
      pv(6,38)=1.
      pv(8,38)=14.
c
      pv(5,19)=atomas(atno2+1.,zno2-1.)
      pv(6,19)=zno2-1.
      pv(5,29)=pv(5,23)
      pv(6,29)=1.
      pv(8,29)=14.
      pv(5,39)=pv(5,23)
      pv(6,39)=1.
      pv(8,39)=14.
c
      nopt=30
      goto 70
   60 pv(5,11)=atomas(atno2+4.,zno2+2.)
      pv(6,11)=zno2+2.
      pv(5,21)=0.
      pv(6,21)=0.
      pv(8,21)=1.
      pv(5,31)=0.
      pv(6,31)=0.
      pv(8,31)=1.
c
      pv(5,12)=atomas(atno2+3.,zno2+2.)
      pv(6,12)=zno2+2.
      pv(5,22)=rmass(16)
      pv(6,22)=0.
      pv(8,22)=16.
      pv(5,32)=0.
      pv(6,32)=0.
      pv(8,32)=1.
c
      pv(5,13)=atomas(atno2+3.,zno2+1.)
      pv(6,13)=zno2+1.
      pv(5,23)=rmass(14)
      pv(6,23)=1.
      pv(8,23)=14.
      pv(5,33)=0.
      pv(6,33)=0.
      pv(8,33)=1.
c
      pv(5,14)=atomas(atno2+2.,zno2+1.)
      pv(6,14)=zno2+1.
      pv(5,24)=rmass(30)
      pv(6,24)=1.
      pv(8,24)=30.
      pv(5,34)=0.
      pv(6,34)=0.
      pv(8,34)=1.
c
      pv(5,15)=atomas(atno2+1.,zno2+1.)
      pv(6,15)=zno2+1.
      pv(5,25)=rmass(31)
      pv(6,25)=1.
      pv(8,25)=31.
      pv(5,35)=0.
      pv(6,35)=0.
      pv(8,35)=1.
c
      pv(5,16)=pv(5,2)
      pv(6,16)=pv(6,2)
      pv(5,26)=rmass(32)
      pv(6,26)=2.
      pv(8,26)=32.
      pv(5,36)=0.
      pv(6,36)=0.
      pv(8,36)=1.
c
      pv(5,17)=atomas(atno2+2.,zno2+2.)
      pv(6,17)=zno2+2.
      pv(5,27)=pv(5,22)
      pv(6,27)=0.
      pv(8,27)=16.
      pv(5,37)=pv(5,22)
      pv(6,37)=0.
      pv(8,37)=16.
c
      pv(5,18)=pv(5,14)
      pv(6,18)=pv(6,14)
      pv(5,28)=pv(5,22)
      pv(6,28)=0.
      pv(8,28)=16.
      pv(5,38)=pv(5,23)
      pv(6,38)=1.
      pv(8,38)=14.
c
      pv(5,19)=atomas(atno2+2.,zno2   )
      pv(6,19)=zno2
      pv(5,29)=pv(5,23)
      pv(6,29)=1.
      pv(8,29)=14.
      pv(5,39)=pv(5,23)
      pv(6,39)=1.
      pv(8,39)=14.
c
      nopt=40
   70 qv     =ek+pv(5,2)+pv(5,1)
      tc     =   pv(6,2)+pv(6,1)
      qval(1)=qv - pv(5,11)
      tch (1)=tc - pv(6,11)
      qval(2)=qv - pv(5,12) - pv(5,22)
      tch (2)=tc - pv(6,12) - pv(6,22)
      qval(3)=qv - pv(5,13) - pv(5,23)
      tch (3)=tc - pv(6,13) - pv(6,23)
      qval(4)=qv - pv(5,14) - pv(5,24)
      tch (4)=tc - pv(6,14) - pv(6,24)
      qval(5)=qv - pv(5,15) - pv(5,25)
      tch (5)=tc - pv(6,15) - pv(6,25)
      qval(6)=qv - pv(5,16) - pv(5,26)
      tch (6)=tc - pv(6,16) - pv(6,26)
      qval(7)=qv - pv(5,17) - pv(5,27) - pv(5,37)
      tch (7)=tc - pv(6,17) - pv(6,27) - pv(6,37)
      qval(8)=qv - pv(5,18) - pv(5,28) - pv(5,38)
      tch (8)=tc - pv(6,18) - pv(6,28) - pv(6,38)
      qval(9)=qv - pv(5,19) - pv(5,29) - pv(5,39)
      tch (9)=tc - pv(6,19) - pv(6,29) - pv(6,39)
   74 qv = 0
      if(irec.eq.2) qval(1)=0.
      if(ipart.ne.16) goto 75
      call grndm(rndm,2)
      if(rndm(1).gt.((atno2-1.)/230.)**2) qval(1)=0.
      eka=7.9254/atno2
      if(rndm(2).lt.ek/eka) goto 75
      qval(3)=0.
      qval(4)=0.
      qval(5)=0.
      qval(6)=0.
      qval(9)=0.
   75 do 71 i=1,9
      if(pv(5,10+i).lt.0.5) qval(i)=0.
      if(qval(i).lt.0.    ) qval(i)=0.
      if(abs(tch(i)-0.1).gt.0.5 ) qval(i)=0.
      qv=qv+qval(i)
   71 continue
      call grndm(rndm,1)
      ran=rndm(1)
      qv1=0.
      do 72 i=1,9
      if(qval(i).eq.0.) goto 72
      qv1=qv1+qval(i)/qv
      if(ran.le.qv1) goto 73
   72 continue
c** reaction kinematically not possible
      nopt=0
      go to 9999
   73 nopt=nopt+i
      pv(5,3)=pv(5,10+i)
      pv(6,3)=pv(6,10+i)
      pv(8,3)=0.
      pv(5,4)=pv(5,20+i)
      pv(6,4)=pv(6,20+i)
      pv(8,4)=pv(8,20+i)
      pv(5,5)=pv(5,30+i)
      pv(6,5)=pv(6,30+i)
      pv(8,5)=pv(8,30+i)
      nt=2
      ran=ek*10.
      if(ran.gt.0.5) ran=0.5
      call grndm(rndm,1)
      if(rndm(1).lt.ran) nt=3
      if(mod(nopt,10).ge.7) nt=3
c** calculate cms energy
   80 pv(4,2)=pv(5,2)
      call add(1,2,200)
      pv(1,200)=-pv(1,200)
      pv(2,200)=-pv(2,200)
      pv(3,200)=-pv(3,200)
c** set quantities for phase space routine in cms
      tecm=pv(5,200)
      npg=nt
      kgenev=1
      do 81 i=1,npg
   81 amass(i)=pv(5,2+i)
c --- invoke double precision version of the phase space package ---
      call phpnuc
      do 83 i=1,npg
      do 82 j=1,4
   82 pv(j,2+i)=pcm(j,i)
c** transform into lab.system
      call lor(2+i,200,2+i)
      pv(7,2+i)=tof
   83 continue
c** set charges and particle index for low mass fragments
      if (abs(pv(5,3)-rmass(14)) .lt. 0.0001) go to 84
      if (abs(pv(5,3)-rmass(16)) .lt. 0.0001) go to 85
      if (abs(pv(5,3)-rmass(30)) .lt. 0.0001) go to 86
      if (abs(pv(5,3)-rmass(31)) .lt. 0.0001) go to 87
      if (abs(pv(5,3)-rmass(32)) .lt. 0.0001) go to 88
      goto 89
   84 pv(6,3)=1.
      pv(8,3)=14.
      goto 89
   85 pv(6,3)=0.
      pv(8,3)=16.
      goto 89
   86 pv(6,3)=1.
      pv(8,3)=30.
      goto 89
   87 pv(6,3)=1.
      pv(8,3)=31.
      goto 89
   88 pv(6,3)=2.
      pv(8,3)=32.
   89 ntt=2+nt
      do 90 i=1,ntt
      ipp=ifix(pv(8,i)+0.01)
      if(ipp.eq.0) goto 90
      ek=pv(4,i)-pv(5,i)
      if(i.lt.3) goto 92
      if(ipp.lt.30) goto 92
      call grndm(rndm,1)
      ek=ek*0.5*rndm(1)
   92 if(ek.lt.1.e-6) ek=1.e-6
      pv(5,i)=rmass(ipp)
      pv(4,i)=ek+pv(5,i)
      p=sqrt(abs(pv(4,i)**2-pv(5,i)**2))
      pp=sqrt(pv(1,i)**2+pv(2,i)**2+pv(3,i)**2)
      if(pp.gt.1.e-6) goto 91
      call grndm(rndm,2)
      phinve=twpi*rndm(1)
      cost=-1.+2.*rndm(2)
      if (cost .le. -1.) cost=-1.
      if (cost .ge.  1.) cost= 1.
      rthnve=acos(cost)
      pv(1,i)=sin(rthnve)*cos(phinve)
      pv(2,i)=sin(rthnve)*sin(phinve)
      pv(3,i)=cos(rthnve)
      pp=1.
   91 pv(1,i)=pv(1,i)*p/pp
      pv(2,i)=pv(2,i)*p/pp
      pv(3,i)=pv(3,i)*p/pp
   90 continue
      if(.not.nprt(4)) goto 100
      write(newbcd,1000) xend,yend,zend,ind,nopt
 1000 format(1h ,'nuclear reaction at (x,y,z) ',3(g12.5,1x)
     +,' material ',i5,' nopt ',i5)
      do 95 i=1,ntt
         write(newbcd,1001) i,(pv(j,i),j=1,10)
   95 continue
 1001 format(1h ,i3,1x,10(g10.3,1x))
  100 intct=intct+1.
c** set interaction mode according to gheisha-convention
      nmode=14
c** n-capture : nmode=10
      if(nopt.eq.1) nmode=10
      if(pv(8,3).gt.0.) goto 110
      call setcur(4)
      ntk=ntk+1
      if(nt.eq.3) call settrk(5)
      go to 9999
 110  continue
      call setcur(4)
      ntk=ntk+1
      call settrk(3)
      if(nt.eq.3) call settrk(5)
      call settrk(3)
c
 9999 continue
      return
      end
