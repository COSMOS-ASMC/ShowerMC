*cmz :  3.14/13 28/05/90  10.56.04  by  rene brun
*-- author :
      subroutine chanwx
c ped
c ped ---------------------------------------------------------------
c ped weights for the sampling procedure (added one to each other in
c ped corresp. channels) for nucrin and hadrin (chanwh+chanwn)
c ped calculation of threshold energy of the reaction channels
c ped                    p.pedroni                02/90
c --------------------------------------------------------------------
*      integer*2 ich,ibar,k1,k2
      common/abltis/am(110),ga(110),tau(110),ich(110),
     +              ibar(110),k1(110),k2(110)
*      integer*2 nzk
      common/split/ nzk(460,3),wt(460)
*
      common /redver/ irii(17),ikii(17),ieii(17),
     +        thresh(268)
      common /reach/umo(296),plabf(296),siinh(296),wkh(5184),
     +              nrk(2,268),nure(30,2)
      common/reacn/siin(296),wkn(5184)
      dimension hwt(460)
      dimension hwkh(40),hwkn(40)
      dimension si(5184),sih(5184)
      equivalence (wkn(1),si(1))
      equivalence (wkh(1),sih(1))
      do 10 i=1,5184
         wkn(i)=wkh(i)
   10 continue
c
      ireg=16
      do 100 ire=1,ireg
         iwko=irii(ire)
         iee=ieii(ire+1)-ieii(ire)
         ike=ikii(ire+1)-ikii(ire)
         ieo=ieii(ire)+1
         iika=ikii(ire)
         do 90  ie=1,iee
            sisn=1.e-14
            sish=1.e-14
            sinorc=0.1
            do ik=1,ike
            iwk=iwko+iee*(ik-1)+ie
            if (nrk(2,iika+ik).eq.0) sinorc=1.
            if (sinorc.lt.1..and.ie.gt.3) si(iwk)=si(iwk)*0.1**(ie-3)
            if (sinorc.lt.1..and.ie.gt.8) si(iwk)=1.e-5
            sish=sish+sih(iwk)
            sisn=sisn+si(iwk)*sinorc
            enddo
            siinh(ieo+ie-1)=sish
            siin(ieo+ie-1)=sisn
            sioh=0.
            sion=0.
            if (sish.ge.1.e-12) go to 20
            sish=1.
            sioh=1.
   20       continue
            if (sisn.ge.1.e-12) go to 30
            sisn=1.
            sion=1.
   30       continue
            sinorc=0.1
            do 40 ik=1,ike
               if(nrk(2,iika+ik).eq.0) sinorc=1.
               iwk=iwko+iee*(ik-1)+ie
               sioh=sioh+sih(iwk)/sish
               sion=sion+si(iwk)/sisn*sinorc
               hwkn(ik)=sion
   40       hwkh(ik)=sioh
            do 50 ik=1,ike
               iwk=iwko+iee*(ik-1)+ie
               wkh(iwk)=hwkh(ik)
   50       wkn(iwk)=hwkn(ik)
            iiki=ikii(ire)
            do 80 ik=1,ike
               am111=0.
               inrk1=nrk(1,iiki+ik)
               if (inrk1.gt.0) am111=am(inrk1)
               am222=0.
               inrk2=nrk(2,iiki+ik)
               if (inrk2.gt.0) am222=am(inrk2)
               thresh(iiki+ik)=am111 +am222
               if (inrk2-1.ge.0) go to 70
               inrkk=k1(inrk1)
               amss=5.
               inrko=k2(inrk1)
               do 60 inrk1=inrkk,inrko
                  inzk1=nzk(inrk1,1)
                  inzk2=nzk(inrk1,2)
                  inzk3=nzk(inrk1,3)
                  ams=am(inzk1)+am(inzk2)
                  if (inzk3-1.ge.0) ams=ams+am(inzk3)
                  if (amss.gt.ams) amss=ams
   60          continue
               ams=amss
               if (ams.lt.umo(ieo)) ams=umo(ieo)
               thresh(iiki+ik)=ams
   70          continue
   80       continue
   90    continue
  100 continue
      do 110 j=1,460
  110 hwt(j)=0.
      do 130 i=1,110
         ik1=k1(i)
         ik2=k2(i)
         hv=0.
         do 120 j=ik1,ik2
            hv=hv+wt(j)
            hwt(j)=hv
            ji=j
  120    continue
         if (abs(hv-1.).gt.1.e-4)write(6,10000)
10000 format(' error in hwt because of false use of chanwx')
  130 continue
      do 140 j=1,460
  140 wt(j)=hwt(j)
c
      end
