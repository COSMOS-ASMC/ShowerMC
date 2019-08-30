*cmz :  3.14/16 14/09/89  11.24.42  by  nick van eijndhoven (cern)
*-- author :
      subroutine settrk(ntr)
c
c *** fill the stack via common /event/ ***
c *** instead of the userword, the particle index is stored ***
c *** nve 01-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (10-nov-1983)
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
c
c --- check pv array boundary ---
      if (ntr .le. 200) go to 10
      print 1000,ntr
 1000 format(' *settrk* ntr = ',i3,' would adress outside pv array'/
     $ ' ===> track will not be put on stack and will be lost')
      go to 9999
c
c --- check total number of produced particles ---
 10   continue
      nvedum=ntot+1
      if (nvedum .le. 100) go to 20
      if (nvedum .eq. 101) print 1001,nvedum
 1001 format(' *settrk* storage of particle no. ',i3, 'not allowed'/
     $ ' maximum number of generated particles is 100'/
     $ ' ===> from now on all generated particles will be discarded')
      go to 9999
c
c --- store generated particle on the stack ---
 20   continue
      eve(next   )=xend
      eve(next+ 1)=yend
      eve(next+ 2)=zend
      eve(next+ 3)=rca
      eve(next+ 4)=rce
      eve(next+ 5)=pv(5,ntr)
      eve(next+ 6)=pv(6,ntr)
      eve(next+ 7)=pv(7,ntr)
      eve(next+ 8)=pv(1,ntr)
      eve(next+ 9)=pv(2,ntr)
      eve(next+10)=pv(3,ntr)
      eve(next+11)=pv(8,ntr)
      next=next+12
      ntot=ntot+1
      next1=next-12
      next2=next-1
      ntot1=ntot-1
c
 9999 continue
      return
      end
