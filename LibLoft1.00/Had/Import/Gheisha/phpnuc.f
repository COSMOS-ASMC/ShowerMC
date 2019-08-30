*cmz :  3.14/16 29/10/90  14.42.30  by  rene brun
*-- author :
      subroutine phpnuc
c
c *** double precision version of the phase space routine "phasp"
c *** this routine must be called by the nuclear interaction routine
c *** "nucrec" (see also comments therein). the reason is simply that
c *** energy-momentum calculations are not possible within only
c *** 6 digits of accuracy for total energies
c *** in the order of hundreds of gev (uranium nucleus), compared with
c *** kinetic energies in the order of mev (neutrons, protons and
c *** photons in the reactions a(x,y(gamma,gamma))a'). in the original
c *** gheisha8 code all these calculations are done in double precision
c *** hmf 29-aug-1989 rwth aachen
c
c called by : nucrec
c origin    : h.fesefeldt
c
      implicit double precision (a-h,o-z)
      real*4 rndm(1)
c
      common/prntfl/inbcd,newbcd,inbin,newbin,npevt,nevtp,lprt,nprt(10)
                    logical lprt,nprt
c
      common/nucin /tecm,amass(18),npg,kgenev
      common/nucout/pcm(5,18),wgt
      double precision tecm,amass,pcm,wgt
c
c
c
      save  knt, twopi, ffq
      dimension emm(18)
      dimension rno(50)
      dimension em(18),pd(18),ems(18),sm(18),ffq(18),pcm1(90)
      equivalence (nt,npg),(amass(1),em(1)),(pcm1(1),pcm(1,1))
      data  ffq/0.,3.141592, 19.73921, 62.01255, 129.8788, 204.0131,
     2                       256.3704, 268.4705, 240.9780, 189.2637,
     3                       132.1308,  83.0202,  47.4210,  24.8295,
     4                        12.0006,   5.3858,   2.2560,   0.8859/
      data  knt , twopi /  1 , 6.2831853073 /
c
c --- initialise local arrays and the result array pcm ---
      call vzero(pcm,90)
      call vzero(emm,18)
      call vzero(pd,18)
      call vzero(ems,18)
      call vzero(sm,18)
c
      knt = knt + 1
      if (.not.nprt(3).and..not.nprt(4)) goto 100
      write(newbcd,1200) npg,tecm,(amass(jk),jk=1,npg)
  100 continue
  150 if (nt .lt. 2)  go to 1001
      if (nt .gt. 18)  go to 1002
      ntm1=nt-1
      ntm2=nt-2
      ntp1=nt+1
      ntnm4 = 3*nt - 4
      emm(1)=em(1)
      tm=0.0
      do 200 i=1,nt
      ems(i)=em(i)**2
      tm=tm+em(i)
 200  sm(i)=tm
      wgt=1.
 210  tecmtm=tecm-tm
      if (tecmtm .le. 0.0)  go to 1000
      emm(nt)=tecm
      if (kgenev.gt.1) go to 400
      emmax=tecmtm+em(1)
      emmin=0.0
      wtmax=1.0
      do 350 i=2,nt
      emmin=emmin+em(i-1)
      emmax=emmax+em(i)
  350 wtmax=wtmax*dpdnuc(emmax,emmin,em(i))
      wtmaxq=1.0/wtmax
      go to 455
  400 wtmaxq=tecmtm**ntm2*ffq(nt) / tecm
  455 continue
      do 457 i= 1, ntnm4
      call grndm(rndm,1)
  457 rno(i) = dble(rndm(1))
      if(ntm2) 900,509,460
  460 continue
      call dlpnuc(rno,ntm2)
      do 508 j=2,ntm1
  508 emm(j)=rno(j-1)*(tecmtm)+sm(j)
  509 wgt=wtmaxq
      ir=ntm2
      do 530 i=1,ntm1
      pd(i)=dpdnuc(emm(i+1),emm(i),em(i+1))
  530 wgt=wgt*pd(i)
      pcm(1,1)=0.0
      pcm(2,1)=pd(1)
      pcm(3,1)=0.0
      do 570 i=2,nt
      pcm(1,i)=0.0
      pcm(2,i) = -pd(i-1)
      pcm(3,i)=0.0
      ir=ir+1
      bang=twopi*rno(ir)
      cb=cos(bang)
      sb=sin(bang)
      ir=ir+1
      c=2.0*rno(ir)-1.0
      s=sqrt(1.0-c*c)
      if(i.eq.nt) go to 1567
      esys=sqrt(pd(i)**2+emm(i)**2)
      beta=pd(i)/esys
      gama=esys/emm(i)
      do 568 j=1,i
      ndx = 5*j - 5
      aa= pcm1(ndx+1)**2 + pcm1(ndx+2)**2 + pcm1(ndx+3)**2
      pcm1(ndx+5) = sqrt(aa)
      pcm1(ndx+4) = sqrt(aa+ems(j))
      call dotnuc(c,s,cb,sb,pcm,j)
      psave = gama*(pcm(2,j)+beta*pcm(4,j))
  568 pcm(2,j)=psave
      go to 570
 1567 do 1568 j=1,i
      aa=pcm(1,j)**2 + pcm(2,j)**2 + pcm(3,j)**2
      pcm(5,j)=sqrt(aa)
      pcm(4,j)=sqrt(aa+ems(j))
      call dotnuc(c,s,cb,sb,pcm,j)
 1568 continue
  570 continue
  900 continue
      return
 1000 do 212 i=1,npg
      pcm(1,i)=0.
      pcm(2,i)=0.
      pcm(3,i)=0.
      pcm(4,i)=amass(i)
  212 pcm(5,i)=amass(i)
      wgt=0.
      return
 1001 if(nprt(3).or.nprt(4)) write(newbcd,1101)
      go to 1050
 1002 write(newbcd,1102)
 1050 write(newbcd,1150) knt
      write(newbcd,1200) npg,tecm,(amass(jk),jk=1,npg)
      stop
 1100 format ('0*phpnuc* available energy negative')
 1101 format ('0*phpnuc* less than 2 outgoing particles')
 1102 format ('0*phpnuc* more than 18 outgoing particles')
 1150 format ('0*phpnuc* above error detected in phasp',
     $ ' at call number ',i7)
 1200 format ('0*phpnuc* input data to phpnuc. npg = ',i6/
     $ ' tecm = ',e15.7,' particle masses = ',5e15.7/(42x,5e15.7))
      end
