************************************************************************
*                                                                      *
      subroutine isbert(finput,npart,nomp,ityp,eout,alph,beta,gam)
*                                                                      *
*       main routine of isobert model                                  *
*       last modified by K.Niita on 01/06/2001                         *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       finput  : input data                                           *
*        finput(1)  : target mass                                      *
*        finput(2)  : target charge                                    *
*        finput(3)  : projectile energy (MeV)                          *
*        finput(4)  : 0.0 fixed                                        *
*        finput(5)  : 1.0 fixed                                        *
*        finput(6)  : andit ,option of angular distribution of delta   *
*        finput(7)  : projectile id, ityp - 1                          *
*                                                                      *
*       nomp    : dimension size                                       *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       npart   : number of out going particles                        *
*       ityp(i) : particle type of i-th particle                       *
*       eout(i) : energy of i-th particle (MeV)                        *
*       alph(i),beta(i),gam(i) : unit momentum vector of i-th particle *
*                                                                      *
*      npcle,nhole : particle and hole number for pre-equ.             *
*      efermi      : fermi energy for pre-equ.                         *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /preeq/  npcles, nholes, efermi, atar, ztar

*-----------------------------------------------------------------------

      common/mecc/eray(60),aray(60),bray(60),gray(60),nopart,kind(60)
      common/matis/b(7,100),mont
      common/statis/bm(4),sb(4),sbs(4),mob
      common/hiro/ st(8000)
      common/palcom/npl(2),monin
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1 denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2 amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3 wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
      common/tapes/jkey,kont,b21,mprin
      common/aedit/data1(5),data4(800),ad4(20),sd4(20),aa(180),
     &             transp,nd4(20),n(5),kid(5000),ket(250),kct(100)
      common/switch/reat,vpion,isonsw,idistr,icb,nzman
      common/aden1/x(10),dx(8),ro(8),xr(8),eva(14),pt
      common /qmri/ npcle, nhole
      common /kstep/nst(30), ncol

*-----------------------------------------------------------------------

      integer   kid,data4
      dimension rmass(5)

      dimension finput(7)
      dimension ityp(nomp)
      dimension eout(nomp),alph(nomp),beta(nomp),gam(nomp)

      data rmass/ 938.272, 939.566, 139.568, 134.973, 139.568/

c ----------------------------------------------------------------------
c
8888  mprin=0
      jkey=1
      kont=0
      mont=0
      idistr=0
      vpion=0.0
      isonsw=0
      pt=1.07
      icb=1
      nzman=10000
      dx(1)=-2.5d0
      dx(2)=1.  d0
      dx(3)=0.6 d0
      dx(4)=0.6 d0
      dx(5)=0.6 d0
      dx(6)=0.6 d0
      dx(7)=0.6 d0
      dx(8)=1.  d0
      poe1=0.0
      b21=0.
   50 format(36h    i++       i+        i0        i-)
   51 format(4f10.0)
  152 format(1h1,46x, 'ang. momentum       ',2x, 'lin. momentum'/
     1' casc #   z   a     ex    p   n  p+  p0  p-',
     2 t48,'!i!   iz/!i!          !p!   pz/!p!',9x
     4,'      rand #'/)
  100 format(3f10.5)
  101 format(i10)
  200 do 1002 i=1,8000
 1002 st(i)=0.0
      do 1009 i=1,4
      bm(i)=0.
      sb(i)=0.
      sbs(i)=0.
 1009 continue
      mob = 0.
      do 115 i=1,2
  115 npl(i)=0
      monin=0
      do 2003 i=1,5
 2003 data1(i)=0.0
      do 2004 i=1,5
 2004 n(i)=0
      do 2005 i=1,5000
 2005 kid(i)=0
      do 2006 i=1,250
 2006 ket(i)=0
      do 2007 i=1,100
 2007 kct(i)=0
      do 2008 i=1,800
 2008 data4(i)=0
      do 2009 i=1,20
 2009 nd4(i)=0
      do 2010 i=1,20
 2010 ad4(i)=0.
      do 2011 i=1,20
 2011 sd4(i)=0.
      do 2012 i=1,180
 2012 aa(i)=0.
      transp=0.
      do 1011 i=1,60
      ityp(i) =0
      kind(i) =0
      eout(i) =0.
      eray(i) =0.
      alph(i) =0.
      aray(i) =0.
      beta(i) =0.
      bray(i) =0.
      gam(i)  =0.
      gray(i) =0.
 1011 continue
      nopart= -1
      ncol=0
      do 3003 i=1,14
      csa(i) = 0.
 3003 continue
      do 3004 i=1,16
 3004 icsb(i) = 0
      do 3005 i=1,30
      nst(i) =0
      ws(i)  =0.
      wsp(i) =0.
      wws(i) =0.
 3005 wwsp(i)=0.
      do 3006 i=1,9
      den(i) =0.
      denn(i)=0.
      denp(i)=0.
      cutofa(i)=0.
         do 3007 j=1,2
         ef(j,i) =0.
         pf(j,i) =0.
 3007 continue
 3006 continue
      do 3008 i=1,200
 3008 stp(i) = 0.
      do 3009 i=1,40
 3009 ind(i) = 0
      do 3010 i=1,5
 3010 amass(i)=0.
      do 3011 i=1,15
 3011 gws(i)=0.
      do 3012 i=1,25
 3012  cm(i)=0.
      do 3013 i=1,20
 3013 ela(i)=0.
      do 3014 i=1,21
 3014 cm1(i)=0.
      do 3015 i=1,14
 3015 cm2(i)=0.
      do 3016 i=1,3
      aip(i)=0.
 3016 ain(i)=0.
      aiscat =0.
      aidec  =0.
      aicap  =0.
      aitot  =0.
      trick  =0.
      fscat  =0.
      poe1   =0.
c
      csa(6)=1.0e0
      csa(1)=finput(7)+ 1.0
c     csa(1)=sngl(finput(7))+ 1.0
      inid = int(csa(1) + 0.0001 )
      restm = rmass(inid)
c     csa(2)=sngl(finput(3))+restm
c     csa(3)=sngl(finput(1))
c     csa(4)=sngl(finput(2))
      csa(2)=finput(3)+restm
      csa(3)=finput(1)
      csa(4)=finput(2)
      amasno=csa(3)
      zee   =csa(4)
c     if (zee.eq.6) pt=1.30
      aneu  =amasno -zee
      bet0=zee*938.256 + aneu*939.550 - amasno*931.478
      betz=(zee-1.)*938.256 + aneu*939.550 -(amasno-1.)*931.478
      betn= zee*938.256 + (aneu-1.)*939.550 -(amasno-1.)*931.478
      csa(5)=0.5*(2.*bet0-betz-betn)
c
      cutofa(1)=0.
      cutofa(2)=11.09
      cutofa(3)=17.21
      cutofa(4)=22.18
      vpion = 0.
      model = 72
      isonsw=  0
      mprin =  0
      jkey  =  0
      kont  =  0
      ind(21) = model-70
      npcle = 1
      nhole = 0
c
    3 call den1i
      if(jkey.eq.1) call begin(model)
      if(mprin.ne.0) print 152
      call amfp
      csa(14)=ws(13)
    5 call poe
      if(ind(23)-1)6,13,6
    6 call selpt(gdc9,gdc10)
      if(ind(22)-1) 20,8,20
   20 if(ws(14)-csa(9)) 8,7,7
    7 call rlar
      if(ind(22)-1) 23,6,23
   23 ind21=ind(21)
      go to (21,21,22,21,21,22),ind21
   22 if(ind(20)-1) 21,10,21
   21 if(ind(17)-1)11,12,11
   11 call deltas
   12 call ajtest
      if(icsb(1))6,13,6
   13 ind(7)=1
   26 call edit(ws)
      icsb(14)=icsb(14)-1
      if(icsb(14).ne.0) go to 5
   31 icsb(14)=0
      call edit(ws)
      if (nopart.eq.-1 .and. ind(3).ne.1) then
          nopart=0
          npcle=2
          nhole=1
      endif
      go to 201
      print 50
      print 51,(aip(i),i=1,2),(ain(i),i=1,2)
chiro go to 200
    8 call psel(gdc9,gdc10)
      if(ind(16)-1) 18,19,18
   18 if(ws(3)-10.0)17,9,9
   17 if(ws(3)-3.0)9,16,16
    9 call ainter(gdc9,gdc10)
      if(ind(22)-1) 16,10,16
   16 if(ind(4)-1) 8887,11,10
 8887 continue
         if(ind(4).eq.-2) then
         write(6,*) 'error in ainter -> retry this cascade '
         nopart=-1
         return
         else
         go to 10
         endif
   10 call gamma
      ind(26)=0
24    if(ind(19).eq.1)go to 5
      if(ind(22).eq.1)go to 6
   19 if(ind(20)-1)12,15,12
   15 ind(20)=0
      go to 6
chiro
  201 continue

      do 1012 i=1,60

      ityp(i)=kind(i)
      eout(i)=eray(i)
      alph(i)=aray(i)
      beta(i)=bray(i)
       gam(i)=gray(i)

 1012 continue

      npart = nopart

      efermi = 0.5 * 0.511976 * (ef(1,1) + ef(2,1))
      npcles = npcle
      nholes = nhole

      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine den1i
c ----------------------------------------------------------------------
c     subroutine den1i calculates the boundaries, fermi energy, and
c     density of each region.  there are 8 regions defined by x(i).
c     in the isobar program efp and efn are the fermi kinetic energy of
c     the protons and neutrons.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/hiro/st(8000)
      common/switch/reat,vpion,isonsw,idistr,icb,nzman
      common/aden1/x(10),dx(8),ro(8),xr(8),eva(14),pt
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat

      dimension efp1(8),efn1(8)

c ----------------------------------------------------------------------
c
c  ef(1,i) and ef(2,i) are the kinetic energy of protons and neutrons.
c  denp(i) and denn(i),and den(i),are the density of protons, neutrons,
c  and total density. pf(1,i) and pf(2,i) are fermi momentum for protons
c  and neutrons.
   11 format(4e13.6)
      csa(8)=csa(4)/csa(3)
      amass(1)=1836.14
      amass(2)=1838.67
      amass(3)=273.18
      amass(4)=264.20
      amass(5)=273.18
      icsb(14)=idint(csa(6))
      i=idint(csa(1))
      ws(4)=amass(i)
      c=pt*csa(3)**(1./3.)
c
      x(1)=c+2.5
c-----reat is a rearrangement time
c-----based on fermi velocity vf/c=.258 corresponding to pf=250 mev/c
      reat=.1*.5/.258
c-----
      do 5 k=2,8
      x(k)=x(k-1)-dx(10-k)
c     x(k)=x(k-1)-sngl(dx(10-k))
    5 if(x(k).lt.0.) x(k)=0.
      x(9)=0.0
      enr=0.0
      epr=1.44*csa(4)/x(1)
      if(cutofa(1).gt.0.0) go to 20
      cutofa(1)=csa(5)+epr
      cutofa(2)=csa(5)
      cutofa(4)=2.*csa(5)
      cutofa(3)=dmax1(cutofa(1),cutofa(4))
   20 cutofa(5)=cutofa(1)+cutofa(1)-csa(5)
      cutofa(6)=cutofa(2)
      eva(1)=pt
      eva(2)=csa(2)
      eva(3)=csa(4)
      eva(4)=csa(3)
      eva(5)=csa(5)
      eva(6)=csa(5)
      eva(7)=cutofa(3)
      eva(8)=cutofa(4)
      csa(7)=csa(3)-csa(4)
      eva(9)=enr
      eva(10)=epr
      csa(9)=x(1)
      csa(10)=csa(9)*csa(9)
      csa(11)=1.44*csa(4)/csa(9)
c
      den(8)=3.0*csa(3)/4.0/3.14159/(0.2*(x(3)*x(3)*x(3)+x(4)*x(4)*x(4)+
     1x(5)*x(5)*x(5)+x(6)*x(6)*x(6))+0.075*(x(7)*x(7)*x(7)+x(2)*x(2)*x(2
     2))+0.025*(x(8)*x(8)*x(8)+x(1)*x(1)*x(1)))
      den(7)=0.975*den(8)
      den(6)=0.9*den(8)
      den(5)=0.7*den(8)
      den(4)=0.5*den(8)
      den(3)=0.3*den(8)
      den(2)=0.1*den(8)
      den(1)=0.025*den(8)
      csa(13)=den(8)
      za=csa(4)/csa(3)
      ana=csa(7)/csa(3)
      do2i=1,8
      denp(i)=za*den(i)
      denn(i)=ana*den(i)
    2 continue
      cnst=1.0/8.3776
      do4i=1,8
      pf(1,i)=2426.25*(cubert(cnst*denp(i)))
      pf(2,i)=2426.25*(cubert(cnst*denn(i)))
    4 continue
      do 6 i=1,8
      do 6 ll=1,2
    6 ef(ll,i)=dsqrt(pf(ll,i)*pf(ll,i)+amass(ll)*amass(ll))-amass(ll)
      ef(1,9)=0.0
      ef(2,9)=0.0
      do 88 j=1,8
      jj=9-j
      xr(j)=x(jj)
      ro(j)=den(jj)
      efp1(j)=ef(1,jj)*0.510976
   88 efn1(j)=ef(2,jj)*0.510976
c
c     call backup
c
      vpion=vpion/.510976
      do 10i=1,6
      cutofa(i)=cutofa(i)/0.510976
   10 continue
      csa(5)=csa(5)/0.510976
      csa(2)=csa(2)/0.510976
      csa(11)=csa(11)/0.510976
      ws(3)=csa(1)
      ws(5)=csa(2)
      ws(22)=0.
      ws(23)=0.
c
      go to 555
c
  100 print 18,(xr(i),ro(i),efn1(i),efp1(i),i=1,8)
   18 format(1h0,7x,4hr(i),10x,5hro(i),10x,6hefn(i),10x,6hefp(i)// (2x,f
     111.4,3x,f11.4,5x,f11.4,5x,f11.4))
      print 16
   16 format(1h0,20x,13hcutoff energy,5x,14hbinding energy,5x,17hpotenti
     1al barrier)
      pip=3.1415927*x(1)**2*10.
      print 9876,eva(7),eva(6),eva(10),eva(8),eva(5),eva(9)
 9876 format(1h0,7x,13hfor protons  ,f11.4,7x,f11.4,8x,f14.4 // 8x,13hfo
     1r neutrons ,f11.4,7x,f11.4,8x,f14.4)
      print 8765,pip
 8765 format(24h0geometric cross section,f11.4,4h mb.)
      print 2289,pt
 2289 format (1h0,' radius parameter =',f6.3,'  fermi')
 555  return
      end subroutine


************************************************************************
*                                                                      *
      subroutine begin(model)
c ----------------------------------------------------------------------
c     subroutine begin writes input card information on tape.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/matis/b(7,100),mont
      common/aden1/x(10),dx(8),ro(8),xr(8),eva(14),pt
      common/tapes/jkey,kont,b21,mprin
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat

c ----------------------------------------------------------------------
c
c**** jkey = 1  write output tape for eva (tape no. 1)
c            0  do not write output tape
c**** kont = 1  continue the existing file
c            0  start a new file on the tape
c
      if(kont .eq. 1) go to 3
      mont=0
      nsiso=8
c
c****  write output tape for e v a - tape no.1 (changed 6.5.74)
c****  tape no. 10 is a backup tape of the storage-st
      ntype=idint(csa(1))
      ng=idint(csa(6))
      return
c
    3 read (1) nsiso,ntype
   30 read (1) mont,((b(i,j),i=1,7),j=1,mont)
      if(mont.eq.1.and.b(1,1).eq.-10000.) go to 4
      if(mont.ne.1000) go to 32
      go to 30
    4 mont=0
      b21=b(2,1)
   31 backspace 1
    2 return
   32 read (1) m,b11,b21
      backspace 1
      go to 31
      end subroutine


************************************************************************
*                                                                      *
      subroutine amfp
c ----------------------------------------------------------------------
c     subroutine amfp calculates the mean free path of the particle.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/hiro/ st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      amfpa=20.0*den(8)
c e is the kinetic energy of the particle in mev.
      e=(ws(5)-ws(4))*0.510976
      if(ws(3)-10.0)1,4,4
    1 if(ws(3)-2.0)3,3,5
c statement 3 calculates the cross section for nucleon interaction.
3     continue
      if(ws(3).eq.1.)t=(cii(e)*csa(4)+cij(e)*csa(7))/csa(3)
      if(ws(3).eq.2.)t=(cii(e)*csa(7)+cij(e)*csa(4))/csa(3)
c         call ciic(cx)
c         cpp = cx
c         call cijc(cx)
c         cpn = cx
      if(ws(3).eq.1.)t=(ciicug(e)*csa(4)+cijcug(e)*csa(7))/csa(3)
      if(ws(3).eq.2.)t=(ciicug(e)*csa(7)+cijcug(e)*csa(4))/csa(3)
cii   call look(ws(14),nstep)
cii   if (nstep.ge.1 .and. nstep.le.7) then
cii   if(ws(3).eq.1.)
cii  &   t=(ciim(e,nstep)*csa(4)+cijm(e,nstep)*csa(7))/csa(3)
cii   if(ws(3).eq.2.)
cii  &   t=(ciim(e,nstep)*csa(7)+cijm(e,nstep)*csa(4))/csa(3)
cii   endif
c cross section are needed to calculate the characteristic time of a
c  particle.
   10 amfpa=1.0/(t*amfpa)
      ws(13)=amfpa*ws(5)/dsqrt(ws(5)*ws(5)-ws(4)*ws(4))
    2 return
c statement 5 to not including 4 calculates the cross section for pion
c  interaction.
c the following change was made 10/16/70 because ptotal is tabulated in
c terms of the k.e.in the c.m.not the k.e.in the laboratory
c to calculate the k.e.in the c.m.the average of a proton and a neutron
c mass was used for the nucleon mass.
5     continue
      e=dsqrt(ws(4)*ws(4)+1837.41*(1837.41+2.*ws(5)))-ws(4)-1837.41
      e=e*0.510976
      if(ws(3)-4.0)7,6,7
    7 ind(9)=1
      call ptotal(e,ind(9),ind(10),ind(11),scat,absor)
      ind(9)=0
      ind(10)=1
      call ptotal(e,ind(9),ind(10),ind(11),scatt,absor)
      ind(10)=0
      amfpc2=1.0-csa(8)
      if(ws(3)-5.0)9,8,9
    9 t=scatt*amfpc2+scat*csa(8)
      go to 10
    8 t=scat*amfpc2+scatt*csa(8)
      go to 10
    6 ind(11)=1
      call ptotal(e,ind(9),ind(10),ind(11),t,absor)
      ind(11)=0
      go to 10
c statement 4 to the end of subroutine calculates the cross section for
c  isobar interaction.
    4 if(ws(3)-13.0)12,11,12
   12 if(ws(3)-25.0)13,11,13
   11 ind(12)=1
   13 aiso=dsqrt(ws(5)*ws(5)-ws(4)*ws(4))/ws(5)
      call look(gws(1),istep)
      aiso1=den(istep)
      wsp(3)=1.0
      call calcu(amass(1),ws(5),ws(4),u,v)
      call aitota(e,aiso,aiso1,u,v)
      amfpi1=csa(8)*aitot
      wsp(3)=2.0
      call calcu(amass(2),ws(5),ws(4),u,v)
      call aitota(e,aiso,aiso1,u,v)
      t=(1.0-csa(8))*aitot+amfpi1
      ind(12)=0
      go to 10
      end subroutine


************************************************************************
*                                                                      *
      subroutine poe
c ----------------------------------------------------------------------
c     subroutine poe calculates the point of entry into the nucleus and
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/hiro/st(8000)
      common/angcom/pxut,pyut,pzut,axut,ayut,azut
      common/switch/reat,vpion,isonsw,idistr,icb,nzman
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),denn(
     19),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40
     2),amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),wwsp(30)
     3,aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      ws(3)=csa(1)
c csa(12) contains the running sum of cascades that are calculated.
c     write(6,*) 'running sub of cascades is  ',csa(12)
      csa(12)=csa(12)+1.0
c  initiates the counters, icsb(i),and zero most of the indicators, ind
c  (i).
      do 1 i=1,20
      ind(i)=0
    1 continue
      ind(23)=0
      ind(3)=1
      if(ws(3)-2.0)11,11,14
   14 ws(24)=0.0
      ind(2)=1
      go to 4
   11 ind(2)=1
      if (ws(3)-1.)2,3,2
    3 ws(24)=ef(1,1)
      go to 4
    2 ws(24)=ef(2,1)
    4 ws(21)=0.
      icsb(1)=1
      icsb(2)=1
      icsb(3)=0
      icsb(4)=0
      icsb(5)=80
      icsb(7)=0
      icsb(9)=21
      icsb(11)=0
      icsb(16)=0
      icsb(15)=0
      icsb(6)=icsb(7)*icsb(5)
      icsb(8)=icsb(15)*icsb(5)
      icsb(10)=icsb(11)*icsb(9)
      icsb(12)=icsb(16)*icsb(9)
      icsb(13)=icsb(6)
      ws(13)=csa(14)
      ws(16)=csa(14)
      ws(17)=csa(14)
      ws(18)=csa(14)
      ws(19)=csa(14)
      ws(15)=0.
      i=idint(csa(1))
      ws(4)=amass(i)
      ws(5)=csa(2)
      ws(14)=csa(9)
      ws(25)=csa(10)
      pp=dsqrt(csa(2)*csa(2)-amass(i)*amass(i))
c rand is between zero and one.
c   5 irand=irand*513+300000001
    5 random=unirn(dummy)
c   5 randh =ransu(0)
      randh=random*1.d0
c     write(6,*) randh, random
c first choose x (ws(10)),then y (ws(11)). since the incident particle
c  enters the nucleus along the z-axis, r=sqrt(x*x+y*y). thus if r is
c  less than the radius of the nucleus, csa(9), the incident particle
c  has hit the nucleus.  the square of the nuclear radius is in csa(10).
      ws(10)=csa(9)*(randh+randh-1.0)
c     irand=irand*513+300000001
      random=unirn(dummy)
c     randh =ransu(0)
      randh=random*1.d0
c     write(6,*) 'random chk ',randh
      ws(11)=csa(9)*(randh+randh-1.0)
      poea=ws(10)*ws(10)+ws(11)*ws(11)-csa(10)
      if(poea)6,5,5
    6 ws(12)=-dsqrt(-poea)
      ws(6)=0.
      ws(7)=0.
      ws(8)=1.
      pxut=0.
      pyut=0.
      pzut=pp
      axut=ws(11)*pp
      ayut=-ws(10)*pp
      azut=0.
      ws(9)=0.0
      ind21=ind(21)
c for model without refraction skip refraction at the surface
c  of the nucleus.
      go to (21,21,22,21,21,22),ind21
c following changed to let pi potential=efn 3/2/71
22    if(ws(3)-2.0)19,20,30
30    continue
      if(ind21.eq.6)ws(5)=ws(5)+vpion
c change made 3/23/71 to allow for attractive coulomb field on pi- i-
      if(ws(3).eq.5.)ws(5)=ws(5)+icb*csa(11)
      if(ws(3).eq.3.)ws(5)=ws(5)-icb*csa(11)
      go to 13
c efn and efp are fermi kinetic energy not including mass.
c for models without refraction or reflection, the incident nucleon just
c  gains energy as it enters the nucleus.
   19 ws(5)=ef(1,1)+csa(5)+ws(5)
      go to 13
   20 ws(5)=ef(2,1)+csa(5)+ws(5)
      go to 13
   21 if(ws(3)-2.0)8,7,12
c for models with refraction pions only feel coulomb potential negative
c  pions could suffer negative coulomb scattering .  then ind(23) is set
c  to one in subroutine refra.
c following changed to let pi potential=efn 3/2/71
12    if(ws(3)-4.0)15,35,16
35    if(ind21.lt.4)go to 13
      refra6=-vpion
      go to 9
c refra6 contains the change of potential when the particle comes in.
   15 refra6=-csa(11)
      ws(5)=ws(5)-csa(11)
      if(ind21.lt.4)go to 9
      refra6=refra6-vpion
      go to 9
c change made 3/23/71 to allow for attractive coulomb field on pi- i-
16    refra6=0.0
      ws(5)=ws(5)+csa(11)
      if(ind21.lt.4)go to 9
      refra6=refra6-vpion
      go to 9
    7 refra6=-ef(2,1)-csa(5)
      go to 9
    8 refra6=-ef(1,1)-csa(5)-csa(11)
      ws(5)=ws(5)-csa(11)
    9 call refra(refra6,refra3,refra2)
c     write(6,*) 'now program has called sub.refra in sub.poe.'
   13 poe1=poe1+1.0
      ws(1)=poe1
      ws(2)=5.0
      ws(14)=csa(9)
      do10i=1,14
      st(i)=ws(i)
   10 continue
c for the incident particle st(17),st(21) are zero and st(18),st(19),
c  st(20) are the x,y,z, of the point of entry.
      st(16)=ws(24)
      st(18)=ws(10)
      st(19)=ws(11)
      st(20)=ws(12)
      st(17)=0.
      st(21)=0.0
c when ind(23)=1, that cascade is counted as a transparency.
      if (ind(23)-1) 17,18,17
   18 ws(1)=-ws(1)
      ws(2)=-ws(2)
      st(1)=ws(1)
      st(2)=ws(2)
c energy is changed to outside energy for storing on tape.
c change made 3/23/71 to allow for attractive coulomb barrier
23    st(5)=ws(5)-refra6-csa(11)
      call edit(ws)
   17 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine selpt(gdc9,gdc10)
c ----------------------------------------------------------------------
c     subroutine selpt selects partner of a collision and also advances
c     the particle to a new position.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      do6i=8,18
      ind(i)=0
    6 continue
      ind(20)=0
c***added for high energy
      ind(27)=0
      ind(28)=0
      ind(29)=0
      ind21=ind(21)
c     irand=irand*513+300000001
      random=unirn(dummy)
c     randh =ransu(0)
      randh=random*1.d0
      if(randh-csa(8))1,2,2
    1 wsp(3)=1.
      wsp(4)=amass(1)
      go to 3
    2 wsp(3)=2.
      wsp(4)=amass(2)
    3 wsp(19)=wsp(4)*wsp(4)
c ind(22)=1 whenever isobar is forced to do something on the spot.
c  another partner may be chosen but should not advance the particle.
c otherwise one would be outside the nucleus if isobar is trying to
c  leave but cannot decay because there is not enough momentum to be
c  divided between the two decay particles.
    8 do4i=1,8
      gws(i+4)=ws(i+4)
    4 continue
      gws(1)=ws(14)
      gws(2)=ws(24)
      gws(3)=ws(25)
c     irand=irand*513+300000001
      random=unirn(dummy)
c     rando =ransu(0)
      gdc10=random
c   7 irand=irand*513+300000001
    7 random=unirn(dummy)
c   7 randh =ransu(0)
      randh=random*1.d0
      wsp(15)=randh+randh-1.
c change made 10/16/70 so that isobars below cut-off are not advanced
      if(ind(22).eq.1)go to 9
      ws(20)=(dsqrt(ws(5)*ws(5)-ws(4)*ws(4))/ws(5))*ws(16)
      gdc9=10.*ws(20)
      ws(25)=0.
      do 5i=1,3
      ws(i+9)=gdc9*ws(i+5)+ws(i+9)
      ws(25)=ws(25)+ws(i+9)*ws(i+9)
    5 continue
      ws(14)=dsqrt(ws(25))
    9 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine rlar
c ----------------------------------------------------------------------
c     subroutine rlar takes care of the particle that has advanced out
c     of the nucleus.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/hiro/st(8000)
      common/switch/reat,vpion,isonsw,idistr,icb,nzman
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),denn(
     19),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40
     2),amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),wwsp(30)
     3,aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
      common/kstep/nst(30),ncol
      data qvpn, qvnp / 18.1203, 12.5867/
      data qalpn, qalnp /  5.8023,  1.8960/
      data qfepn, qfenp /  5.4447,  2.9662/
      data qpbpn, qpbnp /  3.6781,  4.2321/
      data q90pn, q90np /  6.9696,  1.5515/
      data q63pn, q63np /  4.1500,  0.0000/
      data q65pn, q65np /  2.1333,  1.3750/
c ----------------------------------------------------------------------
c
      csa11=cutofa(1)-csa(5)
c ind(4)=1 when a particle has advanced out of the nucleus.
      ind(4)=1
c part of change made 10/19/70 to force pi+ below cut-off to be
c reflected
      ipisw=0
c refra6 is for checking the  energy of the particle to see if it can
c  get out.  it is the sum of binding energy plus the fermi potential
c  for neutron.  it is binding energy, fermi potential, and coulomb
c  barrier for proton.
      if(ws(3)-13.0) 20,21,19
c i++ has cutoff of binding energy + 2 e coulomb
c i+ has binding energy + e coulomb
   21 j=5
   23 refra6=cutofa(j)+ws(24)
      go to 26
c change made 3/23/71 to allow for attractive coulomb field on pi- i-
19    if(ws(3)-24.0)42,41,42
c i0 and i- has the same cutoff as neutrons.
   41 j=2
      go to 23
   42 j=1
      go to 23
   20 if(ws(3)-3.0) 27,9,9
c pion0 has no restriction when leaving the nucleus since it feels
c  neither nuclear nor coulomb potential.
9     if(ws(3)-4.0)22,60,28
c changes made 3/2/71 so pions feel a neutron fermi potential
60    if(ind(21).lt.4)go to 18
      refra6=vpion
      go to 26
c csa11 can be different from the coulomb barrier if in the input cards
c  cutofa(1) uses an effective coulomb barrier.
   22 refra6=csa11
      if(ind(21).gt.3)refra6=refra6+vpion
      go to 26
c change made 3/23/71 to allow for attractive coulomb field on pi- i-
28    refra6=csa11
      if(ind(21).gt.3)refra6=refra6+vpion
      go to 26
   27 i=idint(ws(3))
      refra6=cutofa(i)+ws(24)
c 
c q-value correction for c, al
      qv=0.
      IZ=CSA(4)
      IA=CSA(3)
      IF(CSA(1).EQ.1..AND.WS(3).EQ.2.) CALL TEPN(IZ,IA,QV)
      IF(CSA(1).EQ.2..AND.WS(3).EQ.1.) CALL TENP(IZ,IA,QV)
      if( qv .lt. 1.d-10 ) qv = 0.0d0
      if (ncol.eq.1)  refra6=refra6 + qv/0.510976
c 
c                     refra6=refra6 + qv/0.510976
c the following statements bring the particle back to the surface of
c  the nucleus because in advancing the particle it usually is moved to
c  a radial position greater than the radius of the nucleus.
   26 refra1=0.0
      do1i=1,3
      refra1=refra1+ws(i+5)*ws(i+9)
    1 continue
      refra1=refra1+refra1
      refra2=0.
      do2i=1,3
      refra2=refra2+ws(i+9)*ws(i+9)
    2 continue
      refra1=(refra2-csa(10))/refra1
      ws(25)=0.
      do3i=1,3
      ws(i+9)=ws(i+9)-refra1*ws(i+5)
      ws(25)=ws(25)+ws(i+9)*ws(i+9)
    3 continue
      ws(14)=dsqrt(ws(25))
      if(ws(3).gt.10.)refra6=refra6+vpion
      ind21=ind(21)
      go to (30,30,35,30,30,35) ,ind21
c statement 35 up to not including 30 is for models without refraction.
c  then only energy is checked.  for pion0 and pion- since they feel no
c  potential their energies are not checked at all.
c changes made 3/2/71 so pions have a potential
35    continue
c if particle cannot get out for all particles except pion+, let
c  particle fall below cutoff. ind(20)=1. for isobars below cutoff
c  ind(22)=1.  pion+ is reflected if it cannot leave the nucleus.
   37 if(ws(5)-ws(4)-refra6) 44,36,36
c part of change made 10/19/70 to force pi+ below cut-off to be
c reflected
44    if((ws(3).lt.3.).or.(ws(3).gt.5.))go to 51
c change made 3/2/71 so all pions below cutoff are reflected
      ipisw=1
      go to 30
   51 ind(20)=1
      if(ws(3)-10.0) 14,14,49
   49 ind(22)=1
      do 50 i=1,8
      ws(i+4)=gws(i+4)
   50 continue
      ws(14)=gws(1)
      ws(24)=gws(2)
      ws(25)=gws(3)
      go to 14
c statement 36 is for particles with enough energy to leave.
c for pion+ its energy inside and outside the nucleus is the same.
c  so just let it out.
c for other particles the binding energy and the fermi kinetic energy
c  have to be subtracted.
c change made 3/2/71 so correct pion energy outside calc.
36    if((ws(3).lt.3.).or.(ws(3).gt.5.))go to 40
c change made 3/23/71 to allow for attractive coulomb field on pi- i-
      if(ws(3).eq.5.0)ws(5)=ws(5)-icb*csa11
      if(ws(3).eq.3.0)ws(5)=ws(5)+icb*csa11
      if(ind(21).lt.4)go to 18
      ws(5)=ws(5)-vpion
      go to 18
   40 ws(5)=ws(5)-csa(5)-ws(24)
c change made 3/23/71 to allow for attractive coulomb field on pi- i-
      if(ws(3).eq.25.0)ws(5)=ws(5)-icb*csa11
      if(ws(3).eq.13.0)ws(5)=ws(5)+icb*csa11
      if(ws(3).gt.10.)ws(5)=ws(5)-vpion
      go to 18
c in the following statements, the radial component of momentum is
c  calculated, and checked to see if the particle can get out.
   30 refra2=0.0
      do4i=1,3
      refra2=refra2+ws(i+5)*ws(i+9)
    4 continue
      ws3=dsqrt(ws(5)*ws(5)-ws(4)*ws(4))
      refra2=(refra2/csa(9))*ws3
      refra3=refra2*refra2
c part of change made 10/19/70 to force pi+ below cut-off to be
c reflected
      if(ipisw.eq.1)go to 5
      if((ws(3).eq.5.).or.(ws(3).eq.25.))go to 31
   29 if(refra6*refra6+refra2*refra2-2.0*ws(5)*refra6) 5,6,6
c statement 5 and the following is for reflection at the surface.
    5 refra4=2.*(refra2/ws3)
      do7i=1,3
      ws(i+5)=ws(i+5)-refra4*ws(i+9)/csa(9)
    7 continue
      ws(25)=0.
      do8i=1,3
      ws(i+9)=ws(i+9)+refra1*ws(i+5)
      ws(25)=ws(25)+ws(i+9)*ws(i+9)
    8 continue
      ws(14)=dsqrt(ws(25))
      call look(ws(14),istep)
c
      if (istep.le.0) then
      write(6,*) ' warning radius error , r= ',ws(14)
      istep = 1
      endif
c
c statement 11 to 13 is to give back the proper fermi energy for each
c  type of particle.
   11 if(ws(3)-10.0)24,25,25
   24 if(ws(3)-2.0) 13,12,46
   46 ws(24)=0.0
      go to 14
   25 if(ws(3)-14.0)13,13,12
   12 ws(24)=ef(2,istep)
      go to 14
   13 ws(24)=ef(1,istep)
   14 return
c statement 6 and on is for refraction at the surface
c change made 3/23/71 to allow for attractive coulomb barrier
6     continue
      go to 15
c change made 3/23/71 to allow for attractive coulomb field on pi- i-
31    continue
      if((ws(5)-ws(4)).le.refra6)go to 5
      refra6=refra6-csa11
      go to 29
   15 call refra(refra6,refra3,refra2)
c statement 33 to 17 energy outside is calculated.
   33 if(ws(3)-10.0) 47,47,39
   47 if(ws(3)-2.0) 17,18,48
   48 if(ws(3)-4.0) 17,18,16
   39 if(ws(3)-14.0) 34,17,32
   34 ws(5)=ws(5)+csa11+csa11
      go to 18
   32 if(ws(3)-24.0) 18,18,16
c fori- refraction at the surface is influenced by the negative coulomb
c  barrier it feels.  to check if i- has enough to get out it only has
c  to have enough energy to overcome the nuclear potential.
c change made 3/23/71 to allow for attractive coulomb field on pi- i-
16    continue
      ws(5)=ws(5)-csa11
      go to 18
   17 ws(5)=ws(5)+csa11
   18 icsb13=icsb(13)
      st(icsb13+1)=-ws(1)
      icsb(1)=icsb(1)-1
c ind(6) is set to one, if a particle is leaving the nucleus.
c when a particle leaves, edit is called.
      ind(6)=1
      st(icsb13+13)=ws(19)-ws(16)
c if isobar leaves the nucleus, aiout is called to make isobar decay at
c  the surface.
      if(ws(3)-10.0) 45,45,38
   38 call aiout
c if ind(18)=1, then decay is not allowed. then ind(22) is set to one so
c  that another partner is chosen and the isobar can either scatter or
c  capture.
c change made 10/19/70 so that isobars unable to escape because they
c can not decay and isobars below cut-off are treated the same.they are
c taken back to where they were before being advanced and forced to
c do something there.
      if(ind(18).eq.1)go to 49
      go to 14
   45 ws(21)=ws(21)+1.0
c     nnst = ws(21)
c     nst(nnst) = ncol
      call edit(ws)
      go to 14
      end subroutine


************************************************************************
*                                                                      *
      subroutine deltas
c ----------------------------------------------------------------------
c     subroutine delta is a bookkeeping routine.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
c ind(1)=1 when particle is in collision mode.
      if(ind(1)-1)1,2,1
c in statement 1 the new x,y,z,r,r*r, and potential of the particle
c  is stored in st.
    1 icsb13=icsb(13)
      st(icsb13+10)=ws(10)
      st(icsb13+11)=ws(11)
      st(icsb13+12)=ws(12)
      st(icsb13+14)=ws(14)
      st(icsb13+15)=ws(25)
      st(icsb13+16)=ws(24)
c ind(4)=1 when there will be no collision, or when rlar is entered.
c  in both instances the energy, direction cosines, etc. have to be
c  stored in st.
      if(ind(4)-1)6,4,6
    4 ind(4)=0
      st(icsb13+5)=ws(5)
      st(icsb13+6)=ws(6)
      st(icsb13+7)=ws(7)
      st(icsb13+8)=ws(8)
      st(icsb13+9)=ws(9)
      st(icsb13+15)=ws(25)
      st(icsb13+16)=ws(24)
      if(ind(6)-1)6,8,6
    8 st(icsb13+2)=-ws(2)
    9 return
    6 if(ws(13)-ws(18))10,9,9
   10 ws(18)=ws(13)
      go to 9
c ind(6)=1 when particle escapes the nucleus.  then as in 13, all
c  quantities are stored into st.
    2 if(ind(6)-1)11,12,11
   12 ws(1)=-ws(1)
      ws(2)=-ws(2)
      go to 13
   11 call amfp
      if(ws(13)-ws(18))14,13,13
   14 ws(18)=ws(13)
   13 icsb13=icsb(13)
      st(icsb13+1)=ws(1)
      st(icsb13+2)=ws(2)
      st(icsb13+3)=ws(3)
      st(icsb13+4)=ws(4)
      st(icsb13+5)=ws(5)
      st(icsb13+6)=ws(6)
      st(icsb13+7)=ws(7)
      st(icsb13+8)=ws(8)
      st(icsb13+9)=ws(9)
      st(icsb13+10)=ws(10)
      st(icsb13+11)=ws(11)
      st(icsb13+12)=ws(12)
      st(icsb13+13)=ws(13)
      st(icsb13+14)=ws(14)
      st(icsb13+15)=ws(25)
      st(icsb13+16)=ws(24)
      icsb(3)=icsb(3)-1
      if(icsb(3))9,16,9
   16 icsb(16)=0
      icsb(11)=0
      icsb(10)=icsb(9)*icsb(11)
      ind(1)=0
      go to 9
      end subroutine


************************************************************************
*                                                                      *
      subroutine ajtest
c ----------------------------------------------------------------------
c     subroutine ajtest loads either from st or stp into ws.  if during
c     one time interval, there was one or more collisions and there was
c     one or more particles above the cutoff, then ajtest loads ws from
c     stp.  otherwise it loads ws from st.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
c icsb(1) is the no. of particles to be considered in next time interval
      if(icsb(1))1,13,1
c icsb(2) is the no. of particles still to be considered in the current
c  time interval.
    1 icsb(2)=icsb(2)-1
      if(icsb(2))4,3,4
    3 icsb(2)=icsb(1)
      ws(17)=ws(18)
      ws(19)=ws(18)+ws(19)
      icsb(7)=0
      go to 5
c if ind(1)=1 then ajtest loads ws from stp.   otherwise load ws from
c  st.
    4 if (ind(1)-1)7,6,7
    7 icsb(7)=icsb(7)+1
    5 icsb6=icsb(7)*icsb(5)+1
      if(st(icsb6))7,8,8
    8 icsb(13)=icsb6-1
      icsb13=icsb(13)
      ws(1)=st(icsb13+1)
      ws(2)=st(icsb13+2)
      ws(3)=st(icsb13+3)
      ws(4)=st(icsb13+4)
      ws(5)=st(icsb13+5)
      ws(6)=st(icsb13+6)
      ws(7)=st(icsb13+7)
      ws(8)=st(icsb13+8)
      ws(9)=st(icsb13+9)
      ws(10)=st(icsb13+10)
      ws(11)=st(icsb13+11)
      ws(12)=st(icsb13+12)
      ws(13)=st(icsb13+13)
      ws(14)=st(icsb13+14)
      ws(25)=st(icsb13+15)
      ws(24)=st(icsb13+16)
      ws(16)=ws(17)
      if(ws(3).lt.10.)go to 10
      ws(26)=st(icsb13+77)
      ws(27)=st(icsb13+78)
      ws(28)=st(icsb13+79)
      ws(29)=st(icsb13+80)
      go to 10
    6 icsb(10)=icsb(9)*icsb(11)
      icsb10=icsb(10)
      ws(1)=stp(icsb10+1)
      ws(2)=stp(icsb10+2)
      ws(3)=stp(icsb10+3)
      ws(4)=stp(icsb10+4)
      ws(5)=stp(icsb10+5)
      ws(6)=stp(icsb10+6)
      ws(7)=stp(icsb10+7)
      ws(8)=stp(icsb10+8)
      ws(9)=stp(icsb10+9)
      ws(10)=stp(icsb10+10)
      ws(11)=stp(icsb10+11)
      ws(12)=stp(icsb10+12)
      ws(13)=stp(icsb10+13)
      ws(14)=stp(icsb10+14)
      ws(25)=stp(icsb10+15)
      ws(24)=stp(icsb10+16)
      iws15=dint(ws(15)/100.)
      ws(15)=stp(icsb10+17)-iws15*100
      icsb(15)=icsb(15)+1
      icsb(13)=icsb(15)*icsb(5)
      icsb13=icsb(13)
      st(icsb13+17)=stp(icsb10+17)
      st(icsb13+18)=stp(icsb10+18)
      st(icsb13+19)=stp(icsb10+19)
      st(icsb13+20)=stp(icsb10+20)
      st(icsb13+21)=stp(icsb10+21)
      icsb(11)=icsb(11)+1
      ws(16)=ws(13)
   10 if(ind(6)-1)13,14,13
   14 ws(18)=ws(13)
      ind(6)=0
   13 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine edit(w)
c ----------------------------------------------------------------------
c     subroutine edit does three things.
c     (1) whenever a particle escapes, its identity, energy and cosine
c         with respect to the z-axis are stored for later use.
c     (2) when a cascade is over edit does some partial editing
c         such as calculating zae.
c     (3) when all the cascades have been calculated, a final editing is
c         done.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/qmri/npcle, nhole

      common/mecc/eray(60),aray(60),bray(60),gray(60),nopart,kind(60)
      common/hiro/st(8000)
      common/palcom/npl(2),monin
      common/angcom/pxut,pyut,pzut,axut,ayut,azut
      common/statis/bm(4),sb(4),sbs(4),mob
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
      common/aedit/data1(5),data4(800),ad4(20),sd4(20),aa(180),
     &             transp,nd4(20),n(5),kid(5000),ket(250),kct(100)
      common/matis/b(7,100),mont
      common/tapes/jkey,kont,b21,mprin
      common/ptype/text(5)
      common /kstep/nst(30),ncol

      dimension w(30)
      dimension kint(100)
      integer   data4,kid
c  ----- for q-value test (only for 12c )
c     data qvpn /18.1203/
c     data qvnp /12.5867/
c ----------------------------------------------------------------------
c
  159 format(36h1e' distribution,a=-1 to 18,e'=25mev)
  161 format(28h av. no. of struck protons =,f10.4/28h av. no. of struck
     1 neutrons=,f10.4)
  152 format(1x,i5,86x,i14  )
c 153 format(      return  ,5(2x,i2),2(2f7.1,8x),5x,i14
  153 format(1x,2i5,i4,f8.2,5(2x,i2),2(2f7.1,8x),5x,3z4)
c when icsb(14)=0, statement 2 starts final editing.
      if(icsb(14))1,2,1
c ind(7)=1 when a cascade is finished.  statement 6 starts partial
c  editing.
    1 if(ind(7)-1)3,6,3
c ind(3)=1 when there has been no interaction.  statement 5 takes care
c  of transparencies.
    3 if(ind(3)-1)4,5,4
c statement 4 stores the identity, energy and cosine w.r.t. z-axis of
c  the escaping particle into aa for later use.
4     m=idint(6.0*ws(21)-6.0)
      aa(m+4)=w(4)
      aa(m+5)=w(6)
      aa(m+6)=w(7)
      aa(m+1)=w(3)
      aa(m+2)=w(5)-w(4)
      aa(m+3)=w(8)
      ppr= dsqrt(w(5)*w(5) - w(4)*w(4))
      pxu =w(6)*ppr
      pyu =w(7)*ppr
      pzu =w(8)*ppr
      pxut=pxut-pxu
      pyut=pyut-pyu
      pzut=pzut-pzu
      axut=axut -w(11)*pzu +w(12)*pyu
      ayut=ayut -w(12)*pxu +w(10)*pzu
      azut=azut -w(10)*pyu +w(11)*pxu
    7 return
c for transparencies there is no partial editing.
    5 transp=transp+1.0

cKN
c     icsb(14)=icsb(14)+1
cKN

      iws1= idint(abs (ws(1)) )
      go to 7
    6 if(ind(3)-1)100,7,100
c for incident nucleons the intiial excitation energy is the sum of in-
c  cident kinetic energy and binding energy.  for incident pions the
c  initial excitation energy is just the kinetic energy.
  100 continue
      k80=icsb(15)*80+80
cprt  if(mprin.ne.0) print 1988,((st(i+j-80),j=1,31),i=80,k80,80)
 1988 format (20(1x,16f8.2/3(1x,5f8.2)/))
      if(csa(1)-3.0) 8,14,14
    8 ex=csa(2)+csa(5)
      imass=idint(csa(1))
      ex=ex-amass(imass)
      go  to 15
   14 ex=csa(2)
c change made 3/23/71 to allow for pi-being accelerated by the coulomb
c barrier
      if(csa(1).eq.5.)ex=ex+csa(11)
15    m=idint(6.*ws(21)-5.)
  203 do 9 i=1,5
      n(i)=0
    9 continue
      if(csa(1)-2.0)10,11,12
c iq and ia are the z and a of compound nucleus.
   10 ia=idint(csa(3))+1.0
      iq=idint(csa(4))+1.0
      go to 38
   11 ia=idint(csa(3))+1.0
      iq=idint(csa(4))
      go to 38
   12 ia=idint(csa(3))
      if(csa(1)-4.0)17,18,19
   17 iq=idint(csa(4))+1.0
      go to 38
   18 iq=idint(csa(4))
      go to 38
   19 iq=idint(csa(4))-1.0
   38 if(ws(21)) 204,205,204
c statement 204 to 16, the z,a.ex of residual nucleus are calculated.
204   do 20 i=1,m,6
c data1(k) is a cumulative sum of all types of emitted particles.  for
c  example, when proton is emitted, data1(1) is increased by one.  data1
c  (k) is for all cascades.
      k=idint(aa(i))
      data1(k)=data1(k)+1.0
c n(k)is similar to data1(k) except it is only for one cascade at a
c  time. at the beginning of a cascade, n(k) are set to zero.
      n(k)=n(k)+1
      if(k-2) 21,22,23
   21 ia=ia-1
      iq=iq-1
      ex=ex-aa(i+1)-csa(5)
      go to 20
   22 ia=ia-1
      ex=ex-aa(i+1)-csa(5)
      go to 20
   23 if(k-4) 24,16,26
   24 iq=iq-1
      ex=ex-aa(i+1)-273.18
      go to 20
   26 iq=iq+1
      ex=ex-aa(i+1)-273.18
c change made 12/2/71 to correct for pi- escaping
      ex=ex-csa(11)
      go to 20
   16 ex=ex-aa(i+1)-264.20
   20 continue
c
c data2 bins the kinetic energy of emitted particles in 20 mev bins for
c each type of particle.each particle can have a maximum of 1000 mev.
c data3 bins the cosine of the outgoing particle w.r.t. z-axis, in 20
c  equal bins,from cosine -1.0 to +1.0.
      nopart=0
      do 27 i=1,m,6
      k=idint((aa(i+1)*.510976)/10.+1.)
      ke=min0(k,50)
      k1=idint((aa(i)-1.0)*1000.0)
c out
      qsign=1.0
      iqaa=  int(aa(i))
      qkq=aa(i+1)*0.510976
      nopart=nopart+1
c
c for step wise
c     if (nst(nopart).eq.1 .and. iqaa.eq.1) iqaa=13
c     if (nst(nopart).eq.2 .and. iqaa.eq.1) iqaa=14
c     if (nst(nopart).eq.3 .and. iqaa.eq.1) iqaa=15
c     write(20) qsign,iqaa,qkq,aa(i+2)
c     write(6,*) 'escape ',  qsign,iqaa,qkq,aa(i+2)
c out
      kind(nopart)=iqaa-1
      eray(nopart)=qkq
      gray(nopart)=aa(i+2)
      aray(nopart)=aa(i+4)
      bray(nopart)=aa(i+5)
c
      if(aa(i+2))29,30,31
   30 icos=10
      go to 32
   29 kcos=idint(dabs(aa(i+2)))
      icos=kcos*10.0
      icos=max0(10-icos,1)
      go to 32
   31 icos=idint(aa(i+2)*10.0)
      icos=min0(11+icos,20)
   32 kndex=k1+(ke-1)*20+icos
      kid(kndex)=kid(kndex)+1
      ke=idint(50.*(aa(i)-1.))+ ke
      ket(ke)=ket(ke)+1
      icos=idint(20.*(aa(i)-1))+ icos
      kct(icos)=kct(icos)+1
   27 continue
c data4 bins the excitation energy of the residual nucleus and its a in
c  25 mev interval.  20 emitted particles, and up to 1500mev excitation
c  energy can be accommodated by this table.
  205 axut=axut*0.511/197.32858
      ayut=ayut*0.511/197.32858
      azut=azut*0.511/197.32858
      pxut=pxut*0.511
      pyut=pyut*0.511
      pzut=pzut*0.511
      bm(3)=dsqrt(pxut**2+pyut**2+pzut**2)
      bm(4)=pzut/bm(3)
      bm(1)=dsqrt(axut**2+ayut**2+azut**2)
      bm(2)=azut/bm(1)
      mob=mob+1
      do 1006 i=1,4
      sb(i)=sb(i)+bm(i)
 1006 sbs(i)=sbs(i)+bm(i)**2
      if(jkey.eq.0) go to 201
      mont=mont+1
      b(1,mont)=ex*0.511
      b(2,mont)=dble(ia)
      b(3,mont)=dble(iq)
      b(4,mont)=bm(3)
      b(5,mont)=bm(4)
      b(6,mont)=bm(1)
      b(7,mont)=bm(2)
c     if(mprin.ne.0) print 1238,mont,(b(i,mont),i=1,7)
c1238 format (' eva tape ',i10,1x,7f12.3/)
      if(mont.lt.100) go to 201
c     write(1) mont,b
      mont=0.
c
  201 ex=ex*0.510976
      iex=idint(ex/25.0+1)
      iex=max0(1,iex)
      iex=min0(iex,20)
      a=csa(3)-ia
      if(a.gt.18.or.a.lt.-1) go to 206
      j=idint(a+2.01)
      nd4(j)=nd4(j)+1
      ad4(j)=ad4(j)+ex
      sd4(j)=sd4(j)+ex*ex
      i=idint(a+1.)*40+iex
      data4(i)=data4(i)+1
c
  206 npc=icsb(15)+1
      jn=0
      do 9990 i=1,npc
      i2=80*(i-1)
      kint(i)=idint(dabs(st(i2+2))+21.1)
      j1=jn+1
      jn=jn+kint(i)
      j2=jn-5
      do 9995 i3=j1,j2
      i4=i3-j1+1
 9995 st(i3)=st(i2+i4)
      st(j2+1)=st(i2+77)
      st(j2+2)=st(i2+78)
      st(j2+3)=st(i2+79)
      st(j2+4)=st(i2+80)
 9990 continue
      iws1= idint(abs(ws(1)))
c     if(mprin.ne.0) print 153,iws1,iq,ia,ex,(n(i),i=1,5),
c    1  (bm(i),i=1,4),irand
c     if(mprin.ne.0) print 153,iws1,iq,ia,ex,(n(i),i=1,5),
      m=idint(ws(21)*6.)
      if(m.eq.0) m=1
c
c     nopart=1
      zpr= dble(iq)
      apr= dble(ia)
      if(ex.lt.0.) nopart=-1111
c****  print final results
c
      nplt=npl(1)+npl(2)  !cKN
      npcle=monin-nplt+1
      nhole=npcle-1
    2 return
c   2 elst=poe1-transp
 1001 qcross=pip/poe1
       write(6,*) 'cross section=  ',qcross,'(mb/interaction)'
      write(15,8881)       80.0, elst, transp, pip, 1.0
8881  format(1h ,5f12.3)
      print 34,poe1,elst,transp
   34 format (1h1,5x,'r e s u l t s'//' number of cascades    =',f10.0//
     * ' number of non elastic cascades=',f10.0//
     * ' number  of  transparencis     =',f10.0/)
      do 33i=1,5
   33 data1(i)=data1(i)/elst
      print 154,(text(i),data1(i),i=1,5)
  154 format ('0 avrage number of ',a7,'=',f10.4)
c
      do 212 i=1,4
      sb(i)=sb(i)/mob
  212 sbs(i)=dsqrt(sbs(i)/mob-sb(i)**2)
      print 213,sb,sbs,mob
  213 format (1h0//10x,'ang. momentum',8x,'lin. momentum'/
     *  10x,'!i!    iz/!i!',8x,'!p!    pz/!p!'//
     *  ' mean',4x,2f7.2,7x,2f7.2//' s.d.',4x,2f7.2,7x,2f7.2//
     *  ' no. ',18x,i6)
      nplt=npl(1)+npl(2)
      print 217,monin,nplt,npl(1)
  217 format ('0 total no of interactions',i10/
     * '  pauli principle prevents ',i9,' interactions   (',i6,
     *       ' are nucleon-nucleon interactions.)'/)
c
      do 56 i=1,5
      do 58 j=1,20
      if(kct(20*(i-1)+j).ne.0) go to 61
   58 continue
      go to 56
   61 continue
      print 9992
 9992 format (1h1,'energy and angular distribution of particles relative
     * to beam direction'/)
      print 510,text(i)
  510 format(1h0,10x,33hcosine interval=0.1,from -1 to +1/10x,39henergy
     1interval=10mev,from 0 to 500 mev/3x,a8,'distribution' )
      print 511,(   ii,ii=1,20)
  511 format(1h0,6x,20(i2,3x),5x,5htotal/)
      j2=i*1000
      j1=j2-999
      l=0
      do 57 j=j1,j2,20
      l=l+1
      jh=j+19
      k=jh/20
      print 512,l,(kid(kod),kod=j,jh),ket(k)
      if(l.eq.25) print 515,(ii,ii=1,20)
   57 continue
   56 continue
c
  512 format(i4,20i5,i10 /)
  515 format(1h1/2x,12hcontinuation/7x,20(i2,3x),5x,5htotal/)
      print 520
  520 format (1h1)
      do 530 i=1,5
      print 518,text(i),(j,j=1,20)
  518 format ('0angular distribution of ',a8,'in cos. interval = 0.1
     *from  -1 to  +1'//1x,20(3h j=,i2))
      print 519,  (kct((i-1)*20+j),j=1,20)
  519 format (1x,20i5/)
  530 continue
c
      print 159
      ka=idint(csa(3)+1.01)
      do 163 i=1,20
      if(nd4(i).eq.0) go to 163
      ad4(i)=ad4(i)/nd4(i)
      sd=sd4(i)/nd4(i)-ad4(i)*ad4(i)
      sd4(i)=dsqrt(sd)
      i2=i*40
      i1=i2-39
      print 160,ka,nd4(i),ad4(i),sd4(i),(data4(j),j=i1,i2)
  160 format('0a=',i3,78x,'n(a)=',i5,'    e*(a)=',f8.2,'    s.d.=',f8.2/
     * / 5x,20i5/5x,20i5/)
  163 ka=ka-1
      ws(22)=ws(22)/elst
      ws(23)=ws(23)/elst
      print 161,ws(22),ws(23)
      ws(30)=ws(30)/elst
      print 1005,ws(30)
1005  format('0avg.no.isobar-nucleon charge exchanges=',f10.4)
c
      if(jkey.eq.0) return
c     if(mont.ne.0) write(1) mont,((b(i,j),i=1,7),j=1,mont)
c     print 1980,mont
 1980 format ('0 b matrix is written on a tape.  mont=',i6)
      mont=1
      b(1,1)=-10000.
      b(2,1)=poe1+b21
c     write (1) mont,(b(i,1),i=1,7)
c     end file 1
      go to 7
      end subroutine


************************************************************************
*                                                                      *
      subroutine psel(gdc9,gdc10)
c ----------------------------------------------------------------------
c     subroutine psel selects the momentum of the partner of collision.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/hiro/st(8000)
      common/switch/reat,vpion,isonsw,idistr,icb,nzman
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),denn(
     19),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40
     2),amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),wwsp(30)
     3,aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
      common/aden1/x(10),dx(8),ro(8),xr(8),eva(14),pt
      common/pinin/tcm,i9,i10,i11
      common/tapes/jkey,kont,b21,mprin
c ----------------------------------------------------------------------
c
   59 format(6e13.6)
      ipisw=0
      ind21=ind(21)
c gws(1) and istep are radius and region number before advancing the
c  particle to its new position.
c given a radius, subroutine look finds the corresponding region.
      call look(gws(1),istep)
c
      if (istep.le.0) then
      write(6,*) ' warning radius error , r= ',gws(1)
      istep = 1
      endif
c
c ws(14) and istep1 are the new radius and the new region number after
c  advancing.
      call look(ws(14),istep1)
c
      if (istep1.le.0) then
      write(6,*) ' warning radius error , r= ',ws(14)
      istep = 1
      endif
c
c statement 7 is for no change of region.
c statement 8 is for going toward the center of the nucleus.
    6 if(istep-istep1)8,7,9
c istep2 is the new region number.
    8 istep2=istep1
      go to 10
    9 istep2=istep1+1
   10 ws(14)=x(istep2)
      ws(25)=ws(14)*ws(14)
      step3=0.
      do11i=1,3
      step3=step3+gws(i+5)*gws(i+9)
   11 continue
      gf=dsqrt(step3*step3-gws(3)+ws(25))
      if(ws(25)-gws(3))12,13,13
c gdc9 is the distance to the new region.  see note in notebook.
c if during this distance traveld, the particle is going to cross into
c  a new region, the particle is brought to the boundary and the
c  distance up to the boundary is the distance it is going to travel
c  during this interval.  at the beginning of the next time interval
c  this particle starts at the boundary.  also during crossing in model
c  with refraction a particle is refracted or reflected at the boundary
c  depending on whether it has enough energy to cross or not.
   12 gdc9=-gf-step3
      go to 14
   13 gdc9=gf-step3
c subroutine dis calculates the new x,y,z, coordinate.
   14 call dis(gdc9,psel1,psel2)
c ind(8)=1 when the particle is going to cross.
      ind(8)=1
    7 wsp(18)=den(istep)
c gdc(10) is a random number between zero and one.
   48 ll=wsp(3)
      wsp(24)=ef(ll,istep)
      wsp(9)=pf(ll,istep)*gdc10**(1./3.)
      wsp(5)=dsqrt(wsp(9)*wsp(9)+amass(ll)*amass(ll))
c in psel if particle is pion then select partner and check if there
c  will be an interaction.  if there will be, then pion makes isobar.
c  otherwise pion advances in the nucleus.  statement 50 starts pion
c  interaction.
   17 if(ws(3)-2.0)18,19,52
   52 if(ws(3)-10.0) 53,51,51
   51 if(ws(3)-23.0)18,19,19
c if pion is going to cross then make sure it crosses by going to
c  statement 27.  then go to statement 50 which takes care of pion
c  interaction.
   53 if(ind(8)-1)50,55,50
   55 ind(8)=0
c change made 3/2/71 to allow for a pi potential
      gdc6=0.0
      go to (25,25,1,25,25,1),ind21
   18 ws(24)=ef(1,istep)
      go to 20
   19 ws(24)=ef(2,istep)
   20 if(ind(8)-1)21,22,21
   22 ind(8)=0
      if(ws(3)-2.0)23,24,74
   74 if(ws(3)-15.0)23,23,24
c gdc(6) is the change of potential between regions.
   23 gdc6=ef(1,istep1)-ws(24)
      go to  (25,25,1,25,25,1),ind21
   24 gdc6=ef(2,istep1)-ws(24)
      go to  (25,25,1,25,25,1),ind21
c for step without refraction when particle is crossing change energy
c  but not directions.
    1 ws(5)=ws(5)+gdc6
c if the energy of a particle (isobar) becomes negative when crossing
c  ind(22)=1.  then isobar is forced to do something in the old region.
      if(ws(5).gt.ws(4))go to 27
      if(ws(3).gt.10.)go to 2
c changes made to reflect a pion that can not cross in a no refraction
c reflection option
      ipisw=1
      ws(5)=ws(5)-gdc6
      go to 26
    2 ind(22)=1
      iws1=abs (ws(1))
    5 format(70h isobar mass became more than its total energy during cr
     1ossing        ,3z4,8x,i10)
c   5 format(70h isobar mass became more than its total energy during cr
c    1ossing        i20,i10)
c change made 2/2/71
      ws(5)=gws(5)
      beta=dsqrt(gws(5)*gws(5)-ws(4)*ws(4))/gws(5)
      ws(20)=.1*dsqrt((gws(10)-ws(10))*(gws(10)-ws(10))+
     1(gws(11)-ws(11))*(gws(11)-ws(11))+
     2(gws(12)-ws(12))*(gws(12)-ws(12)))
      ws(16)=ws(16)-ws(20)/beta
      gdc9=10.*ws(20)
      go to 21
   25 if(gdc6)26,27,28
   27 if(gdc9+step3)29,30,30
   29 if(ws(14)-psel2)31,31,33
c the addition or subtraction of 1.e-6 fermi is to insure that the
c  particle is on the other side of the new region.
   31 gdc9=gdc9+1.e-6
      call dis(gdc9,psel1,psel2)
      go to 27
   30 if(psel2-ws(14))31,31,33
   33 ws(25)=psel1
      ws(14)=psel2
      if(ws(3)-10.0)49,46,46
   49 if(ws(3)-2.0)46,46,50
   46 go to 21
c ind(5)=1 when the change of potential is plus or zero.
   28 ind(5)=1
c statement 26 starts calculating whether the particle has enough radial
c  momentum to cross.
   26 gdc7=0.
      do 34i=1,3
      gdc7=gdc7+ws(i+5)*ws(i+9)
   34 continue
      gdc4=dsqrt(ws(5)*ws(5)-ws(4)*ws(4))
      gdc2=(gdc7/ws(14))*gdc4
      psel4=gdc6*gdc6+2.0*gdc6*ws(5)+gdc2*gdc2
c changes made to reflect a pion that can not cross in a no refraction
c reflection option
      if(ipisw.eq.0)go to 80
      print 150
150   format(' pion had to be reflected in psel')
      go to 35
80    continue
      if(psel4)35,36,36
c statement 36 begins calculating refraction at the boundary.
   36 psel4=dsqrt(psel4)
      if(ind(5)-1)37,38,37
   38 psel4=-psel4
      ind(5)=0
   37 gdc1=psel4-gdc2
      ws(5)=ws(5)+gdc6
      gdc7=dsqrt(ws(5)*ws(5)-ws(4)*ws(4))
      do39i=1,3
      gdc3=gdc1*ws(i+9)/ws(14)
      ws(i+5)=(gdc4*ws(i+5)+gdc3)/gdc7
   39 continue
      go to 27
c statement 35 begins calculating reflection at the boundary.
   35 gdc4=gdc2/gdc4
      gdc4=(gdc4+gdc4)/ws(14)
      do40i=1,3
      ws(i+5)=ws(i+5)-ws(i+9)*gdc4
   40 continue
   41 if(gdc9+step3)43,42,42
   43 if(psel2-ws(14))44,44,33
   44 gdc9=gdc9-1.e-6
      call dis(gdc9,psel1,psel2)
      go to 41
   42 if(ws(14)-psel2)44,44,33
   21 return
c statement 50 starts to calculate pion interaction to make isobar.
c to make isobar, it is not sufficient to choose angle between the two
c  partice.  theta and phi have to be chosen.  see notes.
c  50 irand=irand*513+300000001
   50 random=unirn(dummy)
c  50 rando =ransu(0)
      wsp(6)=random
      wsp(6)=wsp(6)+wsp(6)-1.0
      wsp(8)=dsqrt(1.0-wsp(6)*wsp(6))
c     irand=irand*513+300000001
      random=unirn(dummy)
c     randh =ransu(0)
      randh=random*1.d0
      wsp(7)=6.28318*randh
      wsp(29)=dcos(wsp(7))
      wsp(28)=dsin(wsp(7))
      wsp(7)=wsp(8)*wsp(29)
      wsp(8)=wsp(8)*wsp(28)
      wsp(15)=0.0
      do 56i=1,3
      wsp(15)=wsp(15)+gws(i+5)*wsp(i+5)
   56 continue
      gws3=dsqrt(gws(5)*gws(5)-ws(4)*ws(4))
      wsp3=dsqrt(wsp(5)*wsp(5)-wsp(4)*wsp(4))
c knowing the angles of both particles, e is kinetic energy in center of
c  mass.
      e=(dsqrt(ws(4)*ws(4)+wsp(19)+2.0*gws(5)*wsp(5)-2.0*wsp(15)*gws3*ws
     1p3)-wsp(4)-ws(4))*0.510976
c present ptotal has cross sections up to 660 mev if e greater than
c 660 mev choose another set of theta and phi
      if(e-660.)57,58,58
c change made 11/10/70
   58 print 59,ws(5),ws(4),wsp(5),wsp(4),wsp(15),e
      go to 50
c ind(9)=1 for ii collisions.
c ind(10)=1 for ij collisions.
c ind(11)=1 for 0 collisions.
c ptotal calculates pion interaction cross sections.
   57 if(ws(3)-4.0)60,61,62
   62 if(wsp(3)-2.0)64,63,64
   64 ind(10)=1
      go to 65
   63 ind(9)=1
      go to 65
   60 if(wsp(3)-2.0)63,64,63
   61 ind(11)=1
65    continue
      if(e.le.250.)go to 310
      i9=ind(9)
      i10=ind(10)
      i11=ind(11)
      tcm=e
310   continue
      call ptotal(e,ind(9),ind(10),ind(11),scat,absor)
c probability of interaction is calculated just as in subroutine ainter.
      wsp(20)=scat*wsp(18)/gws3/wsp(5)
      wsp(21)=ws(4)*ws(4)*wsp(19)
      selpt1=-gws3*wsp3*wsp(15)+gws(5)*wsp(5)
      wsp(21)=dsqrt(selpt1*selpt1-wsp(21))*wsp(20)
      selpt1=ws(20)*wsp(21)
c     irand=irand*513+300000001
      random=unirn(dummy)
c     rando =ransu(0)
      wsp(20)=random
c change made 3/23/71 to put d restriction on pis
      if(ifrist(wsp(3),ws(2),icsb(13),wsp(18),gws,ws(19),st).gt.1)goto67
      if(wsp(20)-selpt1)66,66,67
c ind(16)=1 when an isobar is made.
   66 ind(16)=1
      if(e.gt.250.)ind(16)=0
      ws(20)=wsp(20)/wsp(21)
      wsp(21)=ws(20)/(dsqrt(gws(5)*gws(5)-ws(4)*ws(4))/gws(5))
c ws(16) is the time remaining in time interval.
210   continue
      ws(16)=ws(16)-wsp(21)
      wsp(13)=ws(16)
c advance the pion to collision site.
      gdc9=ws(20)*10.0
      ws(25)=0.0
      do 68i=1,3
      ws(i+9)=gws(i+9)+gdc9*gws(i+5)
      ws(25)=ws(25)+ws(i+9)*ws(i+9)
   68 continue
      ws(14)=dsqrt(ws(25))
c give it back its original direction cosines.
      do69i=1,4
      ws(i+4)=gws(i+4)
   69 continue
      if(e.gt.250.)go to 300
c subroutine save is called to save some pion characteristics in wws.
ch    call save
      call saveis
c identity of the isobar is calculated. i+ and i0 forget how they are
c  made.
      ws(3)=10.0*wsp(3)+ws(3)
      if(ws(3)-15.0) 101,100,101
 100  ws(3)=24.0
c delta is the fudge factor when there is a difference in potential
c  after an interaction.  see notes.
      delta=ef(2,istep)-ef(1,istep)
      go to 112
  101 if(ws(3)-23.0) 103,102,103
  102 ws(3)=14.0
      delta=ef(1,istep)-ef(2,istep)
      go to 112
  103 delta=0.0
  112 if(ws(3)-14.0) 104,104,105
  104 ws(24)=ef(1,istep)
      go to 106
  105 ws(24)=ef(2,istep)
c energy is adjusted to take care of change in potential.
c following changes made 10/12/70 to conserve energy and at least
c the total momentum direction
106   continue
      ws3=dsqrt(ws(5)*ws(5)-ws(4)*ws(4))
      wsp3=dsqrt(wsp(5)*wsp(5)-wsp(4)*wsp(4))
c combine two particles into the isobar. see notes.
      pi6=0.0
      do70i=1,3
      ws(i+5)=ws(i+5)*ws3+wsp(i+5)*wsp3
      pi6=pi6+ws(i+5)*ws(i+5)
   70 continue
      ws3=sqrt (pi6)
      do71i=1,3
      ws(i+5)=ws(i+5)/ws3
   71 continue
      ws(5)=ws(5)+wsp(5)
c additional changes made 10/16/70.see notes
      cprime=(ws(5)+delta)/ws(5)
      ws(4)=cprime*sqrt(ws(5)*ws(5)-ws3*ws3)
      ws(5)=ws(5)+delta
c change made 3/23/71 to only follow isobars whose mass is greater
c than 1080 mev.force others to capture
      if(ws(4).lt.2113.6)ind(22)=1
      ind(3)=0
c ws(15) has time of birth of isobar.
      ws(15)=ws(19)-ws(16)
c ws(9)=0.0 because isobar was not selected from the fermi sea.
      ws(9)=0.0
c the nucleon that combined with pion is stored in st.
      icsb(15)=icsb(15)+1
      icsb8=icsb(15)*icsb(5)
      st(icsb8+1)=-ws(1)
      st(icsb8+2)=5.0
      st(icsb8+3)=wsp(3)
      st(icsb8+4)=wsp(4)
      st(icsb8+5)=wsp(5)
      st(icsb8+6)=0.0
      st(icsb8+7)=0.0
      st(icsb8+8)=0.0
      st(icsb8+9)=wsp(9)
      st(icsb8+10)=ws(10)
      st(icsb8+18)=ws(10)
      st(icsb8+11)=ws(11)
      st(icsb8+19)=ws(11)
      st(icsb8+12)=ws(12)
      st(icsb8+20)=ws(12)
      st(icsb8+13)=wsp(13)
      st(icsb8+14)=ws(14)
      st(icsb8+15)=ws(25)
      st(icsb8+16)=wsp(24)
      st(icsb8+17)=ws(15)
      st(icsb8+21)=0.0
c the pion that combined with nucleon is placed in the next free space
c  in st so that isobar can take old place in st that pion occupied.
c change made 2/2/71
      call wwsst
c ind(20)=1 because after an isobar is made, it selects a partner and
c  finishes off the time interval the pion had started.
  108 ind(20)=1
  110 icsb13=icsb(13)
c wssto is called to store isobar characterisics into st.
      call wssto(icsb13)
      st(icsb13+17)=ws(15)
      st(icsb13+18)=ws(10)
      st(icsb13+19)=ws(11)
      st(icsb13+20)=ws(12)
      st(icsb13+21)=ws(5)
      call updeca
      ws(29)=wws(3)
      st(icsb13+77)=ws(26)
      st(icsb13+78)=ws(27)
      st(icsb13+79)=ws(28)
      st(icsb13+80)=ws(29)
      go to 21
c change made 6/22/71 to force a pi which can never escape to form an
c isobar which is below cutoff
67    continue
      picut=ws(4)+vpion
      if(ws(3).ne.4.)picut=csa(11)+picut
      if(gws(5).le.picut)go to 200
c ind(4)=1 when there is no interaction, the next particle from st is
c  loaded into ws.
      ind(4)=1
      go to 21
200   continue
      wsp(20)=0.0
      ws(20)=0.0
      wsp(21)=0.0
      ind(16)=1
      ind(22)=1
      go to 210
c no isobar formation
300   continue
      ind(4)=0
      go to 21
      end subroutine


************************************************************************
*                                                                      *
      subroutine ainter(gdc9,gdc10)
c ----------------------------------------------------------------------
c     subroutine ainter decides whether there will be an interaction or
c     not.  also if there will be a collision, ainter decides where the
c     collision will take place.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c////////////
      real(8)::ptotal  ! to avoid conflict with  sub. name
c////////////
c ----------------------------------------------------------------------
c
   20 format (6e13.6)
   16 format(' ainter u=',f10.5)
   17 format(' ainter v=',f10.5)
c statement 9 starts getting cross sections for isobars.
      if(ws(3)-10.0)6,6,9
c gws(5) is the total energy of the particle before crossing into a
c  new region.  gws(5) should be used because a particle is advanced to
c  the boundary of the new region, and collision takes place in the
c  old region.
    6 gdc4=dsqrt(gws(5)*gws(5)-ws(4)*ws(4))
      wschk=wsp(5)*wsp(5)-ws(4)*ws(4)
      if (wschk.lt.0.) go to 2
      selpt1=gws(5)*wsp(5)-dsqrt(wsp(5)*wsp(5)-wsp(4)*wsp(4))*wsp(15)*g
     1dc4
      e=((selpt1/wsp(4))-wsp(4))*0.510976
c ind(19)=1 when one wants to calculate sigmaii or sigmaij in subroutine
c  totcrs.
c subroutine totcrs calculates sigmaii or sigmaij or total cross section
c c is cross section in unit of mb.
cmedium
c     qqe = (gws(5)-ws(4))*0.510976
      call look(gws(1),nstep)
      if(ws(3).eq.wsp(3)) then
                         c=cii(e)
                    c = ciicug(e)
c       call ciic(cx)
c          c = cx
cii   if (nstep.ge.1 .and. nstep.le.7) c=ciim(e,nstep)
      endif
c
      if(ws(3).ne.wsp(3)) then
                         c=cij(e)
                    c = cijcug(e)
c       call cijc(cx)
c          c = cx
cii   if (nstep.ge.1 .and. nstep.le.7) c=cijm(e,nstep)
      endif
c
c     if(ws(3).eq.wsp(3))c=cii(e)
c     if(ws(3).ne.wsp(3))c=cij(e)
      if (wsp(18) .eq.0.) go to 2222
      if (wsp( 5) .eq.0.) go to 2222
      wsp(20)=c*wsp(18)/gdc4/wsp(5)
   22 wsp(21)=ws(4)*ws(4)*wsp(19)
      wsp(21)=wsp(20)*dsqrt(selpt1*selpt1-wsp(21))
c selpt 1 contains the probability of interaction.  if the random no. is
c  less than selpt1, then there will be a collision.  program then
c  branches to statement 1.
      selpt1=ws(20)*wsp(21)
c check to see wheter collision is allowed by distance restriction
      if(ifrist(wsp(3),ws(2),icsb(13),wsp(18),gws,ws(19),st).gt.1)goto2
c     irand=irand*513+300000001
      random=unirn(dummy)
      randh =random
      if(randh-selpt1)1,1,2
c ind(4)=1 when there will be no collision, then the next particle in
c  st is loaded into ws.
    2 ind(4)=1
    3 return
 2222 ind(4) = -2
      go to 3
c ws(20) contains the distance in unit of 10(-12)cm.
    1 ws(20)=randh/wsp(21)
      ws(16)=ws(16)-ws(20)/(dsqrt(gws(5)*gws(5)-ws(4)*ws(4))/gws(5))
      wsp(13)=ws(16)
c here gdc9 is the distance to the collision site.  the following
c  cards bring the particle to the collision site.
      gdc9=ws(20)*10.
      ws(25)=0.
      do 4i=1,3
      ws(i+9)=gdc9*gws(i+5)+gws(i+9)
      ws(25)=ws(25)+ws(i+9)*ws(i+9)
    4 continue
      ws(14)=dsqrt(ws(25))
      do5i=1,4
      ws(i+4)=gws(i+4)
    5 continue
      ind(4)=0
c***added for pi production
      if(ws(3).gt.10.)go to 3
c     irand=irand*513+300000001
      random=unirn(dummy)
c     randh =ransu(0)
      randh=random*1.d0
      if(ws(3).eq.wsp(3))finel=finlii(e)
      if(ws(3).ne.wsp(3))finel=finlij(e)
      if(randh.gt.finel)go to 3
      ind(28)=1
      go to 3
    9 if(ws(3)-25.0)8,10,8
    8 if(ws(3)-13.0)11,10,11
c ind(12)=1 when ws contains an i++ or i-.
   10 ind(12)=1
   11 gws3=dsqrt(gws(5)*gws(5)-ws(4)*ws(4))
      wsp3=dsqrt(wsp(5)*wsp(5)-wsp(4)*wsp(4))
   21 selpt1=gws(5)*wsp(5)-gws3*wsp3*wsp(15)
      wsp(16)=selpt1/wsp(4)
c e is kinetic energy of the isobar hitting a stationary nucleon in the
c  laboratory system.  since zeev fraenkels table only goes up to
c  tlab=1000mev, if e is greater than 1000 mev, a line is printed, and a
c  new angle of collision(wsp(15)) between isobar and nucleon is chosen.
      e=(wsp(16)-ws(4))*0.510976
   13 call look(gws(1),istep)
      aiso=gws3/gws(5)
      aiso1=den(istep)
      u1=ws(4)*ws(4)+wsp(19)
      u2=gws(5)*wsp(5)-wsp(15)*gws3*wsp3
c u is the total energy in the center of mass in unit of bev.
      u=(dsqrt(u2+u2+u1))*0.000510976
c tables for isobar scattering have a minimun u and v of 2.1 and 1.11
c  and if calculated u and v are less than these then u and v are set
c  equal to 2.1 and 1.11.
      if(u-2.1)14,15,15
   14 u=2.1
c v is the mass of the isobar in unit of bev.
   15 v=0.000510976*ws(4)
      if(v-1.11)18,19,19
   18 v=1.11
   19 call aitota(e,aiso,aiso1,u,v)
      aitot=aiscat+aidec+aicap
c changes made to treat isobar decays independently from isobar-nucleon
c interactions
      wsp(20)=wsp(18)/(gws3*wsp(5))
      pdecay=den(istep)*aidec
      wsp(21)=ws(4)*ws(4)*wsp(19)
      wsp(21)=wsp(20)*dsqrt(selpt1*selpt1-wsp(21))
      pcapt=aicap*wsp(21)
      pscat=aiscat*wsp(21)
      if((ind(22).eq.1).or.(ind(24).eq.1))go to 40
c check to see wheter collision is allowed by distance restriction
c isobars can always decay
      if(ifrist(wsp(3),ws(2),icsb(13),wsp(18),gws,ws(19),st).le.1)goto50
      pcapt=0.0
      pscat=0.0
c  50 irand=irand*513+300000001
   50 random=unirn(dummy)
c  50 randh =ransu(0)
      randh=random*1.d0
      ptotal=ws(20)*(pdecay+pcapt+pscat)
      if(randh.gt.ptotal)go to 2
c there will be an isobar interaction.which type.
      wsp(21)=ptotal/ws(20)
c     irand=irand*513+300000001
      random=unirn(dummy)
c     rant  =ransu(0)
      rant=random*1.d0
      if(rant.gt.(ws(20)*pdecay/ptotal))go to 30
c its a decay
      aitot=aidec
      aicap=0.0
      aiscat=0.0
      go to 1
c there will be an isobar nucleon interaction
30    continue
      aidec=0.0
      aitot=aicap+aiscat
      go to 1
40    continue
c just want the capture and decay probabilities
      aiscat=0.0
      if((aidec+aicap).le.0.0)return
      aidec=pdecay/(pdecay+pcapt)
      aicap=1.-aidec
      aitot=1.
      return
      end subroutine
************************************************************************
*                                                                      *
      subroutine look(arg,ii)
c ----------------------------------------------------------------------
c     subroutine look is to find which region the particle is in.  it is
c     the step number.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/aden1/x(10),dx(8),ro(8),xr(8),eva(14),pt
c ----------------------------------------------------------------------
c
      do 1 i=1,9
      if(arg.gt.x(i)) go to 2
    1 continue
    2 ii=i-1
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine calcu(u3,ws5,ws4,u,v)
c ----------------------------------------------------------------------
c     subroutine calcu calculates the total energy in the center of mass
c     in unit of bev between an isobar and a stationary proton or
c     neutron.  this routine is called only by amfp.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
c u3 contains the mass of the proton or neutron.
c ws5 and ws4 are the total energy and mass of the isobar.
c v is the mass of the isobar in unit of bev.
      u1=ws4*ws4+u3*u3
      u=0.000510976*dsqrt(2.0*u3*ws5+u1)
      if(u-2.1)1,2,2
    2 v=ws4*0.000510976
      if(v-1.11)3,6,6
    6 return
c tables for isobar scattering have a minimum u and v of 2.1 and 1.11
c  and if the calculated u and v are less than these, then u and v are
c  set to 2.1 and 1.11.
    1 u=2.1
      go to 2
    3 v=1.11
      go to 6
      end subroutine


************************************************************************
*                                                                      *
      subroutine aitota(e,beta,dens,u,v)
c ----------------------------------------------------------------------
c     subroutine aitota calculates isobar interaction cross section,
c     that is decay, scattering and capture.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/switch/reat,vpion,isonsw,idistr,icb,nzman
      common/hiro/st(8000)
      common/tapes/jkey,kont,b21,mprin
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
      dimension ener(40),enec(40),scat(21),capa(37),capb(37),capc(37)
      data ener/0.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0
     1,500.0,550.0,600.0,650.0,700.0,750.0,800.0,850.0,900.0,950.0,1000.
     20,19*0.0/,enec/2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.5,2.55,2.6,
     32.65,2.7,2.75,2.8,2.85,2.9,2.95,3.0,3.05,3.1,3.15,3.25,3.3,3.35,3.
     44,3.45,3.5,3.55,3.6,3.65,3.7,3.75,3.85,3.9,3.95,4.0,3*0.0/,scat/
     5245.0,216.1,192.6,172.4,153.5,137.8,124.8,114.3,105.4,97.7,90.9,84
     6.9,79.7,74.9,70.7,67.1,63.7,60.7,57.9,55.3,53.1/,capa/16253.76,14
     7123.25,5729.29,3004.16,1739.6,761.59,217.32,95.61,49.75,28.34,17.1
     87,11.06,7.42,5.36,3.95,2.91,2.42,1.89,1.68,1.42,1.08,1.12,0.97,0.7
     99,0.83,0.72,0.54,0.59,0.58,0.41,0.62,0.47,0.31,0.3,0.26,0.19,0.17/
      data capb/-34794.55,-31778.34,-13042.64,-6933.94,-4068.27,-1778.72
     1,-483.33,-201.23,-97.89,-51.15,-27.65,-15.44,-8.61,-5.15,-2.95,-1.
     238,-1.01,-0.39,-0.4,-0.23,0.21,-0.25,-0.41,-0.17,-0.45,-0.36,-0.07
     3,-0.31,-0.41,-0.06,-0.67,-0.39,-0.06,-0.16,-0.11,0.03,0.02/,capc/
     418664.23,17924.52,7458.38,4028.7,2400.83,1058.04,258.24,119.44,59.
     534,32.35,18.82,11.8,7.85,5.82,4.49,3.48,3.19,2.73,2.64,2.43,2.06,2
     6.25,2.16,1.91,2.0,1.87,1.6,1.68,1.67,1.38,1.7,1.47,1.21,1.17,1.09,
     70.95,0.92/
c ----------------------------------------------------------------------
c
      ind(25)=0
      call look(gws(1),istep)
      if(isonsw.lt.1)go to 2
      aiscat=0.0
      go to 100
    2 if (e-950.)5,5,1
    5 ik=20
      call find(e,ener,ie,ik)
      call binter(e,ener,ie,aint3,aint4,aint5)
c aiscat is the scattering cross section.
      aiscat=scat(ie-1)*aint3+scat(ie)*aint4+scat(ie+1)*aint5
      if(ind(12)-1)6,7,6
    7 if(ws(3)-13.0)9,8,9
    8 if(wsp(3)-1.)11,10,11
c fscat is the branching ratio for exchange scattering.  see zeev
c  fraenkels notes.
c up to statement 12 scattering cross section is calculated.
   10 fscat=1.0
      go to 12
   11 fscat=0.33333
      go to 12
    9 if(wsp(3)-1.)10,11,10
    6 if(ws(3)-14.0) 17,13,17
   13 if(wsp(3)-1.)15,16,15
   16 fscat=0.77778
      go to 12
   15 fscat=0.55555
      go to 12
   17 if(wsp(3)-1.)16,15,16
   12 aiscat=fscat*aiscat
100   continue
c up to statement 49 decay cross section is calculated.
c change made 3/23/71 to force isobars with mass lt 1080 mev to
c capture
      if(ws(4).lt.2113.6)go to 33
      if(ind(12)-1) 3,14,3
   14 if(ws(3)-13.0) 28,27,28
c dmass and pmass are mass of nucleon and pion.
   27 dmass=1836.14
      ws24=ef(1,istep)
      go to 29
   28 dmass=1838.67
      ws24=ef(2,istep)
      go to 29
   29 pmass=273.18
c function en calculates energy of the decay nucleon.  this is done to
c  check if decay is kinematically allowed (the decay nucleon has to
c  stay above the fermi level)
      if(en(dmass,pmass,ws(4),beta)-dmass-ws24) 33,33,34
c ind(25)=1 when i+ to n and pion+ or i0 to n and pion0 not allowed.
c ind(25)=2 when i+ to p and pion0 or i0 to p and pion- not allowed.
c ind(25)=3 when i+ and i0 not allowed to decay into either of the
c  modes.
c ind(25)=0 when all modes are allowed.
c first ask which mode or all modes allowed.  thus ind(25) is set.
c  then calculate the fraction the allowed mode is.  then calculate
c  unreduced decay cross section.  then calculate the reduced decay
c  cross section using the fraction of the allowed mode.
    3 if(ws(3)-14.0) 37,35,37
   35 dmass=1838.67
      pmass=273.18
      ws24=ef(2,istep)
      if(en(dmass,pmass,ws(4),beta)-dmass-ws24) 38,38,39
   38 ind(25)=ind(25)+1
   39 dmass=1836.14
      pmass=264.20
      ws24=ef(1,istep)
      if(en(dmass,pmass,ws(4),beta)-dmass-ws24) 40,40,43
   40 ind(25)=ind(25)+2
   43 if(ind(25)-3) 44,33,44
   44 if(ind(25)) 45,34,45
45    continue
c p1 is cross section for ij cross section.
      p1=1./3.
c p2 is cross section for 0 cross section.
      p2=2./3.
   46 pp1=p1/(p1+p2)
      pp2=p2/(p1+p2)
      if(ind(25)-1) 48,47,48
c p is the fraction of allowed mode.
   48 p=pp1
      go to 34
   47 p=pp2
      go to 34
   37 dmass=1838.67
      pmass=264.20
      ws24=ef(2,istep)
      if(en(dmass,pmass,ws(4),beta)-dmass-ws24) 50,50,51
   50 ind(25)=ind(25)+1
   51 dmass=1836.14
      pmass=273.18
      ws24=ef(1,istep)
      if(en(dmass,pmass,ws(4),beta)-dmass-ws24) 52,52,53
   52 ind(25)=ind(25)+2
   53 if(ind(25)-3) 54,33,33
   54 if(ind(25)) 55,34,55
55    continue
c p1 is cross section for 0 cross section
      p1=2./3.
c p2 is cross section for ij cross section.
      p2=1./3.
      go to 46
c all modes are not allowed.  then decay cross section is set to zero.
   33 aidec=0.0
      go to 42
c the following change was made 1/21/71 because a factor of 10
c was left out of the decay cross section
c i.e.the conversion factor should have been 1.e27/(1.e39*3.e10*.73e-23)
c which is 10./2.19 not 1./2.19 to get the cross section in mb.
c
c energy dependend width for isobars a la ginoccio
c the following change was made 20/8/80
c
   34 an=ws(4)*0.510976/1000.
      r1=0.83
      r2=0.62
      aaa=(an*an-0.9)**2 - 0.0172
      aidec=0.
      if(aaa.lt.0.) go to 59
      q=dsqrt(aaa)/(2.*an)*10.
      q=q/1.97
      width=159.7*(q**3)/(1.+(r1*q)**2+(r2*q)**4)
      aidec=(dsqrt(1.0-beta*beta)*width*10.)/(beta*dens*197.)
c
   59 if(ind(25).lt.1) go to 42
c reduced decay cross section is calculated.
   49 aidec=p*aidec
c fcapt is the factor for capture that cross section obtained by opem
c  programs must be multiplied in order to obtain the contribution of
c  the diagram to the total cross section of the process.
c statement 42 starts calculating capture cross section.
c for i++ and p or i- and n capture is not possible.  cross section for
c  capture is set to zero.
   42 if(ind(12)-1) 18,19,18
   18 fcapt=0.66667
c change made 7/6/71.the capture factors for i+ +n and i0 +p are 4/3
c not 2/3 because the final states for these reactions lead to p + n
c and the factors for the reverse reactions need not be divided by 2.
c see z.fraenkel paper on isobar capture.
      if((ws(3).eq.14.).and.(wsp(3).eq.2.))fcapt=1.33333
      if((ws(3).eq.24.).and.(wsp(3).eq.1.))fcapt=1.33333
   20 ik=36
      if(u-2.1) 30,30,31
   30 a=capa(1)
      b=capb(1)
      c=capc(1)
      go to 32
   31 call find(u,enec,ie,ik)
      call binter(u,enec,ie,aint3,aint4,aint5)
      a=capa(ie-1)*aint3+capa(ie)*aint4+capa(ie+1)*aint5
      b=capb(ie-1)*aint3+capb(ie)*aint4+capb(ie+1)*aint5
      c=capc(ie-1)*aint3+capc(ie)*aint4+capc(ie+1)*aint5
   32 aicap=(a*v*v+b*v+c)*fcapt
   21 aitot=aiscat+aidec+aicap
   22 return
   19 if(ws(3)-13.0)23,24,23
   24 if(wsp(3)-1.)25,26,25
   25 fcapt=2.0
      go to 20
   23 if(wsp(3)-1.)26,25,26
   26 aicap=0.0
      go to 21
    1 if(mprin.ge.9) print 4,e
    4 format (25h energy too large aitota=,f10.4)
      aiscat=0.0
      go to 100
      end subroutine


************************************************************************
*                                                                      *
      subroutine refra(refra6,refra3,refra2)
c ----------------------------------------------------------------------
c     subroutine refra refracts the particle when it enters or leaves
c     the nucleus.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/hiro/ st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),denn(
     19),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40
     2),amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),wwsp(30)
     3,aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
c ind(2)=1 when particle is entering the nucleus.
      if(ind(2)-1)1,2,1
    2 refra2=0.
      do3i=1,3
      refra2=ws(i+5)*ws(i+9)+refra2
    3 continue
      refra2=dsqrt(ws(5)*ws(5)-ws(4)*ws(4))*refra2/csa(9)
      refra3=refra2*refra2
c refra6 is the amount of change of potential.
    1 refra1=refra3-2.0*refra6*ws(5)+refra6*refra6
c refra1 has to be checked, for negative pions this could be negative
c  which corresponds to negative coulomb scattering.  when this happens,
c  this cascade is counted as a transparency.
      if(refra1) 8,8,11
    8 ind(23)=1
      print 12
   12 format(41h pion suffers negative coulomb scattering)
      go to 10
   11 refra1=dsqrt(refra1)
      if(ind(2)-1)4,5,4
    5 ind(2)=0
      refra1=-refra1-refra2
      go to 6
    4 refra1=refra1-refra2
    6 refra4=dsqrt(ws(5)*ws(5)-ws(4)*ws(4))
      ws(5)=ws(5)-refra6
      refra2=dsqrt(ws(5)*ws(5)-ws(4)*ws(4))
      do 7i=1,3
      refra3=refra1*ws(i+9)/csa(9)
      ws(i+5)=(refra3+refra4*ws(i+5))/refra2
    7 continue
   10 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine aiout
c ----------------------------------------------------------------------
c     subroutine aiout is called when an isobar is escaping.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
      common/tapes/jkey,kont,b21,mprin
c ----------------------------------------------------------------------
c
c save quantities pertinent to the isobar(in ws) into gws in case it is
c  not allowed to decay at the surface.
c change made 10/20/70.store the direction cosines and energy
c of isobar inside as well as value of ind(25) before trying to decay
c outside also put energy and direction cosines isobar has outside
c in gws
c
      if(mprin.gt.1) print 1001,ws(3),ws(4),ws(5),ws(24),ws(15)
 1001 format(' aiout ',5f10.4)
      gws2=ws(24)
      gws3=ws(3)
      gws4=ws(4)
      gws5=gws(5)
      gws6=gws(6)
      gws7=gws(7)
      gws8=gws(8)
      ind25=ind(25)
      ind(25)=0
      gws(5)=ws(5)
      gws(6)=ws(6)
      gws(7)=ws(7)
      gws(8)=ws(8)
      wws(4)=ws(4)
c  ind(17)=1 when isobar is escaping
   17 ind(17)=1
      if(ws(3)-25.0)13,12,13
c change made 10/20/70 to recalculate ind(25) for i+ and i0
13    if(ws(3).ne.13.)go to 100
c ind(12)=1 for i++ and i-.
   12 ind(12)=1
   14 call decay(ans)
c ind(18)=1 when in subroutine pih the momentum in center of mass to be
c  divided by the two decay particles is negative.  this of course
c  indicate a not allowed decay.
      if(ind(18)-1) 3,5,3
c change made 10/20/70.if isobar can not decay restore gws,ws,ind(25)
5     continue
      ind(25)=ind25
      ws(3)=gws3
      ws(4)=gws4
      ws(24)=gws2
      ws(5)=gws(5)
      ws(6)=gws(6)
      ws(7)=gws(7)
      ws(8)=gws(8)
      gws(5)=gws5
      gws(6)=gws6
      gws(7)=gws7
      gws(8)=gws8
      go to 6
3     continue
c change made 11/13/70 to restore gws to values that they are assumed
c to have
      gws(2)=gws2
      gws(3)=gws3
      gws(4)=gws4
      if(gws(3)-13.)1,2,1
c aip(1),aip(2),ain(1),ain(2),keep count of no. of i++,i+,i0,i- emitted.
    2 aip(1)=aip(1)+1.0
      go to 11
    1 if(gws(3)-14.0) 7,4,7
    4 aip(2)=aip(2)+1.0
      go to 11
    7 if(gws(3)-24.0) 9,10,9
   10 ain(1)=ain(1)+1.0
      go to 11
    9 ain(2)=ain(2)+1.0
c when an isobar escapes, the time,x,y,z of escape,and energy before
c  decay are stored.
   11 ws(15)=ws(19)-ws(16)
      icsb13=icsb(13)
      st(icsb13+3)=gws(3)
      st(icsb13+4)=gws(4)
      st(icsb13+5)=gws(5)
      st(icsb13+16)=gws(2)
      icsb13=16.0+ws(2)
      icsb13=icsb13+icsb(13)
      ws(2)=-(ws(2)+5.0)
      ws(1)=-ws(1)
      st(icsb13+2)=ws(10)
      st(icsb13+3)=ws(11)
      st(icsb13+4)=ws(12)
      st(icsb13+1)=ws(15)
      st(icsb13+5)=gws(5)
      call wssto(icsb(13))
c the decay of the isobar at the surface is treated in almost exactly
c  way as decay within the the nucleus, except the nucleon from the
c  decay is not checked for pauli forbiddenness.
      call elas(ans)
      call cmb
      call cmc(cm(7),cm(8),cm(9),cm(10))
      ws(5)=cm2(3)
      ws(6)=cm2(4)/cm2(2)
      ws(7)=cm2(5)/cm2(2)
      ws(8)=cm2(6)/cm2(2)
      call cmc(cm(11),cm(12),cm(13),cm(14))
      wsp(5)=cm2(3)
      wsp(6)=cm2(4)/cm2(2)
      wsp(7)=cm2(5)/cm2(2)
      wsp(8)=cm2(6)/cm2(2)
   16 wsp(14)=ws(14)
      do 15i=1,3
      wsp(i+9)=ws(i+9)
   15 continue
      wsp(25)=ws(25)
      wsp(27)=ws(15)
      ws(27)=ws(15)
      ws(2)=5.0
      wsp(2)=5.0
      ws(1)=poe1
      wsp(1)=poe1
c subroutine aiouta is called to store the decay particles into st.
      call aiouta(wsp)
      ws(21)=ws(21)+1.0
c subroutine edit is called to count the decay particle as escapees.
      call edit(wsp)
      ws(21)=ws(21)+1.0
      call edit(ws)
      call aiouta(ws)
    6 return
c change made 10/20/70 to recalculate ind(25) for i+ and i0
100   continue
      if(ws(3).eq.14.)go to 150
c have an i0
      if(ws(4).le.(amass(2)+amass(4)))ind(25)=1
      if(ws(4).le.(amass(1)+amass(5)))ind(25)=2
      go to 14
150   continue
c have an i+
      if(ws(4).le.(amass(2)+amass(3)))ind(25)=1
      if(ws(4).le.(amass(1)+amass(4)))ind(25)=2
      go to 14
      end subroutine


************************************************************************
*                                                                      *
      subroutine dis(gdc9,psel1,psel2)
c ----------------------------------------------------------------------
c     subroutine dis advances the particle by the distance in location
c     gdc9.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      ws(20)=gdc9/10.
      psel1=0.
      do 1i=1,3
      ws(i+9)=gws(i+9)+gdc9*gws(i+5)
      psel1=psel1+ws(i+9)*ws(i+9)
    1 continue
      psel2=dsqrt(psel1)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine saveis
c ----------------------------------------------------------------------
c     subroutine save saves some quantities from ws into wws.  in case
c     the isobar makes a forbidden collision, such as forbidden charge
c     exchange scattering, capture, or decay.  also when a pion tries to
c     interact save ws into wws,so that if an isobar is formed,wws is
c     stored into st via subroutine wwsst.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),denn(
     19),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40
     2),amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),wwsp(30)
     3,aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      do 1 i=3,9
      wws(i)=ws(i)
    1 continue
      wws(24)=ws(24)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine wwsst
c ----------------------------------------------------------------------
c     subroutine wwsst stores some quantities from wws into st.  if a
c     pion interacts forming an isobar, characteristics of pion (stored
c     into wws by subroutine save) are stored into st.  also in allowed
c     isobar exchange scattering, capture and decay, characteristics of
c     the isobar (stored into wws by subroutine save) are stored into
c     st.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/hiro/st(8000)
      common/palcom/npl(2),monin
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      icsb13=icsb(13)
      igf=ws(2)+16.0
      ws(2)=ws(2)+5.0
      ws(1)=-ws(1)
      icsb(15)=icsb(15)+1
      icsb8=icsb(15)*icsb(5)
      st(icsb8+1)=ws(1)
      st(icsb8+2)=ws(2)
      st(icsb8+3)=wws(3)
      st(icsb8+4)=wws(4)
      st(icsb8+5)=wws(5)
      st(icsb8+6)=wws(6)
      st(icsb8+7)=wws(7)
      st(icsb8+8)=wws(8)
      st(icsb8+9)=ws(9)
      st(icsb8+10)=ws(10)
      st(icsb8+11)=ws(11)
      st(icsb8+12)=ws(12)
      st(icsb8+13)=ws(13)
      st(icsb8+14)=ws(14)
      st(icsb8+15)=ws(25)
      st(icsb8+16)=wws(24)
      do 2 i=17,igf
      st(i+icsb8)=st(i+icsb13)
    2 continue
      i8igf=icsb8+igf
      st(i8igf+1)=ws(15)
      st(i8igf+2)=ws(10)
      st(i8igf+3)=ws(11)
      st(i8igf+4)=ws(12)
      st(i8igf+5)=wws(5)
    4 ws(1)=dabs(ws(1))
      ws(2)=5.0
      if(ind(16)-1)6,7,6
    6 if(ind(14)-1)8,10,8
   10 i=idint(wwsp(3))
      go to 11
    8 return
    7 i=idint(wsp(3))
   11 ws(i+21)=ws(i+21)+1.0
      monin=monin+1
      ws(15)=ws(15)+i*100.
      go to 8
      end subroutine


************************************************************************
*                                                                      *
      subroutine wssto(icsb8)
c ----------------------------------------------------------------------
c     subroutine wssto stores some quantities from ws into st.  when an
c     isobar is first made, its characteristics are sotred into st. also
c     when isobar leaves the nucleus, its characteristics (before
c     decaying at the surface) are stored into st.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/hiro/ st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      if(ind(17)-1) 2,1,2
    2 st(icsb8+3)=ws(3)
      st(icsb8+4)=ws(4)
      st(icsb8+5)=ws(5)
    1 st(icsb8+1)=ws(1)
      st(icsb8+2)=ws(2)
      st(icsb8+6)=ws(6)
      st(icsb8+7)=ws(7)
      st(icsb8+8)=ws(8)
      st(icsb8+9)=ws(9)
      st(icsb8+10)=ws(10)
      st(icsb8+11)=ws(11)
      st(icsb8+12)=ws(12)
      st(icsb8+13)=ws(13)
      st(icsb8+14)=ws(14)
      st(icsb8+15)=ws(25)
      st(icsb8+16)=ws(24)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine updeca
c ----------------------------------------------------------------------
c     routine calculate cos(delta), cos(phi), and sin(phi) fo decay of
c     isobars with respect to a specific z-axis in the isobar frame.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c///////////
      real(8)::gamma
c////////////
c ----------------------------------------------------------------------
c
      cosw=0.0
      do 10 i=1,3
      cosw=cosw+ws(i+5)*wws(i+5)
10    continue
      p2=dsqrt(wws(5)*wws(5)-wws(4)*wws(4))
      alphac=ws(5)/ws(4)
      betac=dsqrt(alphac*alphac-1.)
      p2cx=alphac*p2*cosw-betac*wws(5)
      sinw=dsqrt(1.-cosw*cosw)
c  ws(28)=cos(delta)
      ws(28)=p2cx/dsqrt(p2cx*p2cx+p2*p2*sinw*sinw)
      if(sinw.lt.(1.e-7))go to 20
      if(dabs(ws(8)).ge.(.9999999))go to 30
c normal calc. of sin(phi) and cos(phi)
      denom=dsqrt(1.-ws(8)*ws(8))
      beta=(-ws(7)*wws(6)+ws(6)*wws(7))/denom
      gamma=(-ws(8)*(ws(6)*wws(6)+ws(7)*wws(7))+wws(8)*denom*denom)
      gamma=gamma/denom
c ws(26)=cos(phi) ws(27)=sin(phi)
      ws(26)=-gamma/sinw
      ws(27)=beta/sinw
      return
20    continue
c sinw almost zero.hence phi=0
      ws(26)=1.
      ws(27)=0.0
      return
30    continue
c isobar moving almost parallel to lab.z axis
      ws(26)=wws(6)/sinw
      ws(27)=wws(7)/sinw
      if(ws(8).lt.(0.0))ws(26)=-ws(26)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine cma
c ----------------------------------------------------------------------
c     subroutine cma is part of the center of mass transformation as
c     outlined by metropolis, etc, in lams-2360.  the other center of
c     mass routines are cmb, cmc, and elas.  cma brings the particle
c     into the center of mass system.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension dcm(25)
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      gws5=dble(gws(5))
      wsp5=dble(wsp(5))
      wsp4=dble(wsp(4))
      ws4= dble(ws(4))
      dcm(15)=ws4*wsp4
      dcm(16)=gws5*wsp5
      gws3=dsqrt (gws5*gws5-ws4*ws4)
      wsp3=dsqrt(wsp5*wsp5-wsp4*wsp4)
      dcm(25)=(dcm(16)-wsp(15)*gws3*wsp3)/dcm(15)
      dcm(1)=dsqrt (dcm(25)*dcm(25)-1.d0)
      dcm(2)=gws5+wsp5
      dcm(16)=dcm(2)*dcm(2)
      dcm(17)=dcm(25)*dcm(15)
      dcm(18)=ws4*ws4
      dcm(20)=dcm(18)+dcm(17)
      dcm(19)=wsp4*wsp4
      dcm(21)=dcm(19)+dcm(17)
      dcm(22)=dcm(21)+dcm(20)
      dcm(3)=dsqrt (dcm(22))
      dcm(4)=dsqrt (dcm(16)-dcm(22))
      dcm(22)=dcm(4)*dcm(1)
      dcm(16)=dcm(25)*(gws5-wsp5)
      dcm(17)=wsp5*dcm(18)
      dcm(5)=(dcm(16)+(gws5*dcm(19)-dcm(17))/dcm(15))/dcm(22)
      dcm(6)=-dsqrt (1.d0-dcm(5)*dcm(5))
      dcm(16)=dcm(1)*dcm(15)/dcm(3)
      dcm(7)=dcm(16)*dcm(5)
      dcm(11)=-dcm(7)
      dcm(9)=dcm(16)*dcm(6)
      dcm(13)=-dcm(9)
      dcm(10)=dcm(20)/dcm(3)
      dcm(14)=dcm(21)/dcm(3)
      dcm(25)=(wsp(15)*wsp3+gws3)/dcm(4)
      dcm(1)=-dsqrt (1.d0-dcm(25)*dcm(25))
      dcm(8)=0.d0
      dcm(12)=0.d0
      do 15 i=1,25
   15 cm(i)=dcm(i)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine cmb
c ----------------------------------------------------------------------
c     subroutine cmb is part of the center of mass transformation
c     routine.  it sets up the 4 inverse matrices.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      ws(8)=dabs(gws(8))
      ws8=ws(8)
      if(ws8+0.0000001-1.)2,1,1
    2 cm1(7)=gws(8)
      cm1(9)=dsqrt(1.0-gws(8)*gws(8))
      cm1(1)=gws(6)
      cm1(5)=cm1(1)/cm1(9)
      cm1(3)=-cm1(5)*cm1(7)
      cm1(4)=gws(7)
      cm1(2)=-cm1(4)/cm1(9)
      cm1(6)=cm1(2)*cm1(7)
      cm1(8)=0.
8     continue
      if(ind(15).eq.1.)go to 10
c if isobar decay phi and omega rotations are changed
c     irand=irand*513+300000001
      random=unirn(dummy)
c     randh =ransu(0)
      randh=random*1.d0
      cm1(21)=randh*6.28318
      cm1(14)=dcos(cm1(21))
      cm1(17)=dsqrt(1.-cm1(14)*cm1(14))
      if(cm1(21)-3.14159)3,4,4
    4 cm1(17)=-cm1(17)
    3 cm1(13)=cm1(17)*cm(1)
      cm1(12)=cm(1)
      cm1(16)=-cm(1)*cm1(14)
      cm1(10)=cm(25)
      cm1(15)=-cm(25)*cm1(17)
      cm1(18)=cm(25)*cm1(14)
      cm1(11)=0.
      cm1(19)=cm(2)/cm(3)
      cm1(20)=cm(4)/cm(3)
      return
    1 if(gws(8)) 5,6,6
    6 cm1(7)=1.
      cm1(5)=1.
      cm1(3)=-1.
      go to 7
    5 cm1(3)=1.
      cm1(5)=1.
      cm1(7)=-1.
    7 cm1(1)=0.
      cm1(2)=0.
      cm1(4)=0.
      cm1(6)=0.
      cm1(9)=0.
      cm1(8)=0.
      go to 8
10    continue
      if(ws(26).eq.10.)go to 15
      cm1(14)=ws(26)
      cm1(17)=ws(27)
      go to 3
15    continue
      cm1(14)=1.
      cm1(17)=0.0
      go to 3
      end subroutine


************************************************************************
*                                                                      *
      subroutine cmc(cm7,cm8,cm9,cm10)
c ----------------------------------------------------------------------
c     subroutine cmc is part of the center of mass transformation.  the
c     particles after collision are brought back into the laboratory
c     system.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      cm2(8)=cm7*cm1(19)
      cm2(4)=cm10*cm1(20)+cm2(8)
      cm2(5)=cm8
      cm2(6)=cm9
      cm2(8)=cm7*cm1(20)
      cm2(3)=cm10*cm1(19)+cm2(8)
      cm2(9)=0.
      cm2(10)=0.
      cm2(11)=0.
      do1i=1,3
      cm2(9)=cm2(9)+cm2(i+3)*cm1(i+9)
      cm2(10)=cm2(10)+cm2(i+3)*cm1(i+12)
      cm2(11)=cm2(11)+cm2(i+3)*cm1(i+15)
    1 continue
      cm2(4)=0.
      cm2(5)=0.
      cm2(6)=0.
      do4i=1,3
      cm2(4)=cm2(4)+cm2(i+8)*cm1(i)
      cm2(5)=cm2(5)+cm2(i+8)*cm1(i+3)
      cm2(6)=cm2(6)+cm2(i+8)*cm1(i+6)
    4 continue
      cm2(8)=0.0
      do7i=1,3
      cm2(8)=cm2(8)+cm2(i+3)*cm2(i+3)
    7 continue
      cm2(2)=dsqrt(cm2(8))
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine anglen(idp1,idp2,kecm,coscm)
c ----------------------------------------------------------------------
c     program selects the cos in the c.m.system for the elastic
c     scattering of two nucleons.  the cos of the selected angle
c     is coscm,idp1 and idp2 are the identities of the two nucleons
c     colliding.  1.implies protons, 2.implies neutrons.kecm is the
c     kinetic energy in the c.m.of the 2 particles.
c     all selections are based on p.clements and l.winsberg (ucrl-9043)
c     polynomial fits to nucleon-nucleon scattering data.
c     written by george d.harp in march,1969
c     program should be all right for energies greater than 6.15 bev.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      real*8 idp1,idp2,kecm
c////////////
      real(8):: alpha, cii, cij
      eta2 =unirn(dummy)   ! to avoid undefined use 
c///////////
c ----------------------------------------------------------------------
c
c first decide whether part.1 and 2 are different
      if(idp1.ne.idp2)go to 1000
c start i-i selection.first the find the k.e.of a proton striking
c stationary proton that would give the same kecm.this energy,
c ep,will be in units of 939 mev
      ep=kecm*(5.6754243e-7*kecm+2.12992545e-3)
c at low energies (t less than 126.8 mev)have isotropic scattering
      if(ep.gt..135)go to 10
c     irand=irand*513+300000001
      random=unirn(dummy)
c     rando =ransu(0)
      eta=random
      coscm=2.*eta-1.
      return
10    continue
      if(ep.gt.1.065)go to 50
      if(ep.gt..600)go to 40
c now in region t=127 mev to t=1000 mev
      alpha=1.949-0.327*ep
      beta=9.1*ep*ep*ep*ep*ep
      cii=alpha/(alpha+beta)
      dii=beta/(alpha+beta)
      go to 45
40    continue
      alpha=-20.97+ep*(90.52+ep*(-117.25+48.02*ep))
      beta=150.03+ep*(-578.33+ep*(734.67-299.46*ep))
      cii=alpha/(alpha+beta)
      dii=beta/(alpha+beta)
45    continue
c now use von neumann rejection in this region
c 46  irand=irand*513+300000001
  46  random=unirn(dummy)
c 46  rando =ransu(0)
      eta1=random
c     irand=irand*513+300000001
      random=unirn(dummy)
c     rando =ransu(0)
      eta2=random
      coscm=2.*eta2-1.
      if((cii+dii*coscm*coscm*coscm*coscm).ge.eta1)return
      go to 46
50    continue
      if(ep.gt.6.55)go to 60
c now in region t=1000mev to t=6150 mev
      n=-9.68+ep*(16.89+ep*(-4.348+0.4469*ep))
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta1  =ransu(0)
      eta1=random*1.d0
      coscm=eta1**(1./float(n+1))
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta2  =ransu(0)
      eta2=random*1.d0
      if(eta2.lt..5)coscm=-coscm
      return
60    continue
c now in region t greater than 6.15 bev
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta1  =ransu(0)
      eta1=random*1.d0
      coscm=eta1**(1./51.)
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta2  =ransu(0)
      eta2=random*1.d0
      if(eta2.lt..5)coscm=-coscm
      return
c this finishes i-i collisions now start i-j
c
c
c first find the k.e.of a neutron striking a stationary proton that
c give the same kecm.this energy,en,will be in units of 939 mev
1000  continue
      en0=kecm*(5.6754243e-7*kecm+2.1313931e-3)
c for t less than 14 mev scattering is isotropic
      if(en0.gt..0150)go to 70
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta   =ransu(0)
      eta=random*1.d0
      coscm=2.*eta-1.
      return
70    continue
      if(en0.gt..200)go to 90
c now in region from t=14mev to t=188mev
      oneoen=1./en0
      tij=2.10+oneoen*(-.1624+oneoen*(.04011-4.0e-4*oneoen))
      uij=1.57+en0*(251.7+en0*(-5021.+en0*(41192.-1.162e5*en0)))
      if(uij.lt.0.0)uij=0.0
      vij=en0*(-41.0+en0*(2978.+en0*(-33953.+108909.*en0)))
      if(vij.lt.0.0)vij=0.0
      pij=-.90+oneoen*(.377+oneoen*(.0184-1.77e-4*oneoen))
      qij=en0*(182.9+en0*(-1408.8+3066.2*en0))
c the max.occur at cos=1 or -1
      pmax=tij+uij+vij
      if(pmax.lt.(pij+qij))pmax=pij+qij
      aij=tij/pmax
      bij=uij/pmax
      cij=vij/pmax
      dij=pij/pmax
      eij=qij/pmax
      if(uij.eq.0.0)go to 75
      if(vij.eq.0.0)go to 80
c use von neuman rejection for all 3 possible distributions
c 71  irand=irand*513+300000001
  71  random=unirn(dummy)
c 71  eta1  =ransu(0)
      eta1=random*1.d0
      coscm=2.*eta1-1.
       x=dabs(coscm)
c     x=    coscm
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta2  =ransu(0)
      eta2=random*1.d0
      if(coscm.lt.0.0)go to 72
      if((dij+eij*x*x*x*x).ge.eta2)return
      go to 71
72    continue
      if((aij+x*x*(bij+cij*(x**12.))).ge.eta2)return
      go to 71
75    continue
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta1  =ransu(0)
      eta1=random*1.d0
      coscm=2.*eta1-1.
      x=dabs(coscm)
c     x=coscm
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eat2  =ransu(0)
      eta2=random*1.d0
      if(coscm.lt.0.0)go to 77
      if((dij+eij*x*x*x*x).ge.eta2)return
      go to 75
77    continue
      if((aij+cij*(x**14.)).ge.eta2)return
      go to 75
80    continue
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta1  =ransu(0)
      eta1=random*1.d0
      coscm=2.*eta1-1.
c     x=abs(coscm)
      x=    coscm
      if(x.lt.0.0)go to 82
      if((dij+eij*x*x*x*x).ge.eta2)return
      go to 80
82    continue
      if((aij+bij*x*x).ge.eta2)return
      go to 80
90    continue
      if(en0.gt..620)go to 100
c now in region from t=188mev to t=582mev
      wij=-.50+en0*(18.26+en0*(-41.64+25.08*en0))
      xij=37.42+en0*(-227.38+en0*(512.62-380.37*en0))
      yij=-4.42+en0*(35.43+en0*(-78.62+59.52*en0))
      zij=4.0
      rij=1.32+en0*(3.95-6.71*en0)
      sij=4.19+en0*(3.34+en0*(-44.86+62.18*en0))
      pmax=wij+xij+yij+zij
      if(pmax.lt.(rij+sij))pmax=rij+sij
      aij=wij/pmax
      bij=xij/pmax
      cij=yij/pmax
      dij=zij/pmax
      eij=rij/pmax
      fij=sij/pmax
92    continue
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta1  =ransu(0)
      eta1=random*1.d0
      coscm=2.*eta1-1.
      x=dabs(coscm)
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta2  =ransu(0)
      eta2=random*1.d0
      if(coscm.lt.0.0)go to 94
      if((eij+fij*x*x).ge.eta2)return
      go to 92
94    continue
      dijx=0.
      if(x.gt.0.16) dijx=dij*x**96
      alpha=aij+x*x*x*x*(bij+cij*x**10+dijx)
      if(alpha.ge.eta2)return
      go to 92
100   continue
      if(en0.gt.1.065)go to 105
c now t=582mev to t=1000mev
      alpha=-20.84+en0*(88.88+en0*(-114.82+47.08*en0))
      beta=149.74+en0*(-562.82+en0*(704.13-285.13*en0))
      cii=alpha/(alpha+beta)
      dii=beta/(alpha+beta)
c now can go back to routine already written for i-i collisions
      go to 46
105   continue
      if(en0.gt.6.55)go to 110
c now in region t=1000mev to t=6150mev
c the form of the distribution for i-j collisions is the
c same as that for i-i collisions
      ep=en0
      go to 50
110   continue
c now in region t greater than 6.15 bev and the form of distribution
c for i=j collisions same as for i-i
      go to 60
      end subroutine


************************************************************************
*                                                                      *
      subroutine elas(ans)
c ----------------------------------------------------------------------
c     subroutine elas does the elastic collision in the center of mass
c     system as outlined in lams-2360.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/hiro/ st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      ela(10)=ans
      ela(11)=dsqrt(1.d+0-ela(10)*ela(10))
      ela(12)=cm(7)*ela(11)
      ela(19)=cm(9)
      cm(9)=ela(10)*cm(9)+ela(12)
      ela(12)=ela(19)*ela(11)
      cm(7)=ela(10)*cm(7)-ela(12)
      random=unirn(dummy)
      randh=random*1.d0
      ela(19)=randh*6.28318d+0
      ela(13)=dcos(ela(19))
      ela(14)=dsqrt(1.d+0-ela(13)*ela(13))
      if(ela(19)-3.14159265d+0)1,2,2
    2 ela(14)=-ela(14)
    1 ela(4)=ela(14)*cm(6)
      ela(2)=-ela(4)
      ela(8)=ela(14)*cm(5)
      ela(6)=-ela(8)
      ela(5)=ela(13)
      ela(12)=1.d+0-ela(5)
      ela(3)=cm(5)*cm(6)*ela(12)
      ela(7)=ela(3)
      ela(1)=cm(5)*cm(5)*ela(12)+ela(13)
      ela(9)=cm(6)*cm(6)*ela(12)+ela(13)
      ela(16)=0.
      do3i=1,3
      ela(16)=ela(16)+cm(i+6)*ela(i)
    3 continue
      ela(17)=0.d+0
      do4i=1,3
      ela(17)=ela(17)+cm(i+6)*ela(i+3)
    4 continue
      ela(18)=0.d+0
      do 5 i=1,3
      ela(18)=ela(18)+cm(i+6)*ela(i+6)
    5 continue
      cm(7)=ela(16)
      cm(11)=-ela(16)
      cm(8)=ela(17)
      cm(12)=-ela(17)
      cm(9)=ela(18)
      cm(13)=-ela(18)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine resto
c ----------------------------------------------------------------------
c     subroutine resto restores some quantities that have been stored in
c     wws into ws.  for example, if an isobar has made a forbidden
c     charge exchange scattering, capture, or decay, then the
c     characteristics of the original isobar has to be stored back into
c     ws.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),denn(
     19),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40
     2),amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),wwsp(30)
     3,aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      do 1 i=3,9
      ws(i)=wws(i)
    1 continue
      ws(24)=wws(24)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine decay(ans)
c ----------------------------------------------------------------------
c     subroutine decay takes care of the decay of the isobar in its own
c     center of mass.  the argument ans is the angle between the pion
c     and nucleon after decay.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/switch/reat,vpion,isonsw,idistr,icb,nzman
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
c when ind(17)=1, isobar is trying to escape.  therefore it must be on
c  the surface and hence istep=9.  otherwise one has to calculate the
c  region number.
      if(ind(17)-1) 41,40,41
   40 istep=9
      go to 42
   41 call look(gws(1),istep)
c ind(15)=1 for decay.
   42 ind(15)=1
c ind(12)=1 when the decaying particle is i++ or i-.
      if(ind(12)-1)1,2,1
c delta is a fudging quantity used to conserve energy.  when an i++ or
c  i- decays, no fudging is necessary.  i++ feels a proton potential and
c  i- feels a neutron potential.  hence before and after decay there is
c  no change in potential.  however for i+ and i0, there is a
c  possibility of a change of potential because i+ and i0 feel a proton
c  and neutron potential respectively.
    2 delta=0.0
      if(ws(3)-14.0)3,3,4
c statement 3 is for i++ decaying into pion+ and proton.
    3 wsp(3)=1.0
      ws(3)=3.0
      wsp(24)=ef(1,istep)
      go to 15
c statement 4 is for i- decaying into pion- and neutron.
    4 wsp(3)=2.0
      ws(3)=5.0
      wsp(24)=ef(2,istep)
      go to 15
c gf is a random number used later to choose decay mode of i+ and i0.
c   1 irand=irand*513+300000001
    1 random=unirn(dummy)
c   1 gf    =ransu(0)
      gf=random*1.d0
c ind(25)=1 when i+ decaying to n and pion+ or i0 to n and pion0 not
c  allowed kinematically.
c ind(25)=2 when i+ decaying to p and pion0 or i0 to p and pion- not
c  allowed kinematically.
c ind(25)=3 when i+ and i0 are not allowed kinematically to decay into
c  any mode.
c ind(25)=0 when all modes are allowed.
c ind(25) is set in subroutine aitota.
      if(ind(25)-1) 7,10,11
   10 if(ws(3)-14.0) 20,13,20
   11 if(ws(3)-14.0) 21,14,21
c statement 7 is for all modes allowed.
7     continue
c change made 10/21/70.if the c.m.k.e.is less than or equal 140 mev
c than the t=1/2 cross section should be neglible and the probability
c for an i+ or i0 decaying into a pi+ and a neutron or a pi- and a
c proton should be 1./3.
      aisobd=1./3.
      aisob1=ws(3)
      if(ws(3)-14.0)12,9,12
c for i+, if gf is less than aisobd, then goes to pion+ and n.
c  if greater, goes to pion0 and p.
    9 if(gf-aisobd)14,13,13
   14 wsp(3)=2.0
      ws(3)=3.0
      wsp(24)=ef(2,istep)
      delta=ef(2,istep)-ef(1,istep)-icb*csa(11)
      go to 15
   13 wsp(3)=1.0
      ws(3)=4.0
      wsp(24)=ef(1,istep)
      delta=0.0
   15 i=idint(wsp(3))
      wsp(4)=amass(i)
      i=idint(ws(3))
      ws(4)=amass(i)
      ws(24)=0.0
      wsp(9)=0.0
      ws(9)=0.0
      wsp(19)=wsp(4)*wsp(4)
      go to 8
c for i0, if gf is less than aisobd, then goes to pion- and p.
c if greater, goes to pion0 and n.
   12 if(gf-aisobd)20,21,21
   21 ws(3)=4.0
      wsp(3)=2.0
      delta=0.0
      wsp(24)=ef(2,istep)
      go to 15
   20 ws(3)=5.0
      wsp(3)=1.0
      wsp(24)=ef(1,istep)
      delta=ef(1,istep)-ef(2,istep)+icb*csa(11)
      go to 15
c statements 8 to 22 do pseudo center of mass transformation.
c the decaying isobar is already in its own center of mass system.
c  statement 8 to 22 do what subroutine cma does.
    8 cm(2)=gws(5)
      cm(4)=dsqrt(gws(5)*gws(5)-wws(4)*wws(4))
      cm(3)=wws(4)
      do 22i=1,8
      cm(i+6)=0.0
   22 continue
c subroutine pih divides the energy available for decay in c.m. between
c  the two decay particles.  for more details see flow chart for pih.
      cm(25)=1.
      cm(5)=1.
      cm(1)=0.0
      cm(6)=0.0
      if(ws(26).eq.10.)go to 210
      cm(5)=ws(28)
      cm(6)=-dsqrt(1.-ws(28)*ws(28))
210   continue
      call pih(delta)
c ind(18)=1 when the momentum available for decay is negative.  this is
c  an unallowed decay.
      if(ind(18)-1) 6,5,6
c if it is allowed, choose the angle between the two decay partners
c  using subroutine pangle.
6     continue
      if(ws(26).eq.10.)go to 200
      if(ws(29).eq.10.)go to 300
c on isobar exchange force isotropic decay
      if(ws(3).ne.ws(29))go to 30
      if(ws(3)-4.0)26,25,27
   27 if(wsp(3)-2.0)29,28,29
   29 ind(10)=1
      go to 30
   28 ind(9)=1
      go to 30
   26 if(wsp(3)-2.0)28,29,28
   25 ind(11)=1
   30 e=(wws(4)-ws(4)-wsp(4))*0.510976
      call pangle(e,ind(9),ind(10),ind(11),ws(3),ws(29),ans)
    5 return
200   continue
c     irand=irand*513+300000001
      random=unirn(dummy)
c     ans   =ransu(0)
      ans=random*1.d0
      ans=2.*ans-1.
      go to 5
300   continue
c inelastic isobars decay as p2(costheta)
      call angli1(ans)
      go to 5
      end subroutine


************************************************************************
*                                                                      *
      subroutine aiangl(e,u,v,ans,ind13,ind14)
c ----------------------------------------------------------------------
c     aiangle is called for isobar capture and scattering.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension ener(40),enec(40),ctaba(21),ctabb(21),ctabc(21),
     1 taba(21),tabb(21),tabc(21)
      data ener/2.176,2.297,2.317,2.337,2.357,2.377,2.397,2.416,2.435,2.
     1454,2.473,2.492,2.511,2.530,2.549,2.567,2.585,2.603,2.621,2.639,2.
     2656,19*0.0/,enec/2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.55,2.6,2.
     365,2.7,2.75,2.8,2.95,3.1,3.3,3.45,3.6,3.95,4.0,19*0.0/,taba/2.2,
     42.66,3.114,3.573,4.036,4.503,4.974,5.425,5.881,6.34,6.802,7.269,7.
     5739,8.213,8.69,9.146,9.605,10.07,10.53,11.0,11.45/,tabb/1.37,1.85
     62,2.321,2.793,3.269,3.749,4.232,4.695,5.162,5.632,6.105,6.582,7.06
     73,7.547,8.034,8.499,8.967,9.439,9.913,10.39,10.84/,tabc/6.15,6.02
     8,5.837,5.66,5.49,5.327,5.169,5.024,4.885,4.75,4.619,4.493,4.371,4.
     9253,4.139,4.035,3.933,3.835,3.74,3.647,3.563/
      data ctaba/-7.47,-13.13,-8.61,-6.67,-5.43,-4.04,-2.52,-1.84,-1.31,
     1-1.19,-1.04,-0.95,-0.89,-0.84,-0.73,-0.67,-0.62,-0.59,-0.56,-0.56,
     2-0.56/,ctabb/14.07,27.92,18.07,13.94,11.27,8.13,4.55,2.97,1.77,1.
     35,1.17,0.98,0.84,0.72,0.49,0.37,0.25,0.19,0.13,0.12,0.13/,ctabc/-
     46.19,-14.47,-9.0,-6.68,-5.14,-3.24,-1.01,0.03,0.96,1.23,1.56,1.79,
     52.01,2.21,2.77,3.29,3.99,4.53,5.09,6.41,6.60/
c///////////
      real(8)::bb, alpha
c/////////
c ----------------------------------------------------------------------
      ik=0
      ie=0
      if(ind13-1)10,100,10
c for capture coefficients a and b are functions of u and v. see notes.
   10 ik=20
      if(u-2.1) 30,30,31
   30 a=ctaba(1)
      b=ctabb(1)
      c=ctabc(1)
      go to 32
   31 call find (u,enec,ie,ik)
      call binter (u,enec,ie,aint3,aint4,aint5)
      a=ctaba(ie-1)*aint3+ctaba(ie)*aint4+ctaba(ie+1)*aint5
      b=ctabb(ie-1)*aint3+ctabb(ie)*aint4+ctabb(ie+1)*aint5
      c=ctabc(ie-1)*aint3+ctabc(ie)*aint4+ctabc(ie+1)*aint5
   32 aa=(a*v*v)+(b*v)+c
      if(u-2.35)11,12,12
   11 app=-0.4998*u-0.08044
      go to 13
   12 app=-1.25
   13 if(u-3.025)14,15,15
   14 bpp=2.6524*u-4.03749
      go to 16
   15 bpp=3.52817*u-6.72547
   16 bb=app*v+bpp
c having gotten aa and bb (equivalent to a and b in notes) calculate
c following changes made 8/2/71 to symeterize capture differential
c cross sections
      d1=bb+.01949
      dsq=d1*d1
      asq=aa*aa
      alpha=bb*dsq
      beta=asq*(bb-2.*d1)
      y=alpha/(dsq*dsq)
      y1=(alpha+beta)/((dsq-asq)**2)
c calculate maximum y and store in amax.
      if(y-y1)18,17,17
   18 y=y1
c  also calculate y at where first derivative is zero.  that occurs at
17    xsq=dsq*(3.*bb-2.*d1)/(asq*(2.*d1-bb))
      if((xsq.gt.1.).or.(xsq.lt.0.0))go to 19
      y1=(alpha+beta*xsq)/((dsq-asq*xsq)**2)
      if(y1-y)19,20,20
   19 amax=y
      go to 21
   20 amax=y1
c statement 21 starts rejection method.
c  21 irand=irand*513+300000001
   21 random=unirn(dummy)
c  21 randh =ransu(0)
      randh=random*1.d0
      ran1=randh+randh-1.
c     irand=irand*513+300000001
c     ran2  =ransu(0)
      random=unirn(dummy)
      ran2=random*1.d0
      xsq=ran1*ran1
      yr1=(alpha+beta*xsq)/((dsq-asq*xsq)**2)
      yr2=ran2*amax
      if(yr1-yr2)21,22,22
   22 ans=ran1
   23 return
c statement 100 is for scattering.
  100 ik=20
      if(u-2.176) 105,105,104
c coefficients are only function of u.
  105 a=taba(1)
      b=tabb(1)
      c=tabc(1)
      go to 106
  104 call find(u,ener,ie,ik)
      call binter(u,ener,ie,aint3,aint4,aint5)
      a=taba(ie-1)*aint3+taba(ie)*aint4+taba(ie+1)*aint5
      b=tabb(ie-1)*aint3+tabb(ie)*aint4+tabb(ie+1)*aint5
      c=tabc(ie-1)*aint3+tabc(ie)*aint4+tabc(ie+1)*aint5
c y here is y at 180 degrees.
  106 d=(-a+b+0.196)
      y=1.0/(d*d+0.01*c)
c ymax here is y at minimum or maximum (at first derivative=0).
      ymax=1.0/(0.01*c)
      if(y-ymax)101,102,102
  102 ymax=y
c y here is y at 0 degrees.
  101 y=1.0/((a+b+0.196)*(a+b+0.196)+0.01*c)
      if(y-ymax)107,108,108
c amax = is maximum y.
  107 amax=ymax
      go to 103
  108 amax=y
c statement 103 starts rejection method.
c 103 irand=irand*513+300000001
  103 random=unirn(dummy)
c 103 ran1  =ransu(0)
      ran1=random*1.d0
      ran1=ran1+ran1-1.
c     irand=irand*513+300000001
      random=unirn(dummy)
c     ran2  =ransu(0)
      ran2=random*1.d0
      yr1=1./((a*ran1+b+0.196)*(a*ran1+b+0.196)+0.01*c)
      yr2=ran2*amax
      if(yr1-yr2)103,22,22
      end subroutine


************************************************************************
*                                                                      *
      subroutine pih(delta)
c ----------------------------------------------------------------------
c     subroutine pih divides the momentum of the isobar to the two decay
c     particles.
c     pih is also called when isobar does exchange scattering when delta
c     does not equal to zero.  then the momentum of the two scattered
c     particles have to be recalculated.
c     also in isobar capture when delta does not equal to zero, pih is
c     called to recalculate the momentum of the two captured particles.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),denn(
     19),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40
     2),amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),wwsp(30)
     3,aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
c ----------------------------------------------------------------------
c
      e=ws(4)*ws(4)
c cm3 has the energy for decay,that is the mass plus delta(the fugde
c  factor).  cm(3) has to be divided between nucleon and pion requiring
c  the momentum of the two decay particles to be equal.  also see flow
c  chart for pih.
c changes made 10/6/70 to conserve energy and momentum in decay,
c charge exchange scattering,and absorption when a difference in
c potential occurs.
      cm3=cm(3)
c pi6 is the square of the momentum.  if pi6 is negative then decay is
c  not allowed.  then ind(18)=1.
      pi6=(cm3*cm3+wsp(19)-e)/cm3/2.0
      pi6=pi6*pi6-wsp(19)
      if(pi6.gt.0.0)go to 10
    3 ind(18)=1
      go to 5
    4 pi6=sqrt (pi6)
c now pi6 is the momentum
c cm(7),cm(8),cm(9), and cm(10) are the x,y,z,component of momentum in
c  c.m.and energy of particle in ws.  for decay of isobar px=p.  cm(11),
c  cm(12),cm(13),cm(14) are similar quantities for particle in wsp.
    1 cm(7)=cm(5)*pi6
      cm(11)=-cm(7)
      cm(9)=cm(6)*pi6
      cm(13)=-cm(9)
      pi6=pi6*pi6
      cm(10)=dsqrt(pi6+e)
      cm(14)=dsqrt(pi6+wsp(19))
    5 return
10    continue
      if(delta.eq.0.0)go to 4
c now lengthen or shorten the c.m.momentum vector.(see notes)
      c=cm(2)+delta
      cprime=cm(3)*c/cm(2)
      psqc=(cprime*cprime+e-wsp(19))/(2.*cprime)
      psqc=psqc*psqc-e
      pi6=psqc
      if(pi6.gt.0.0)go to 4
      go to 3
      end subroutine


************************************************************************
*                                                                      *
      subroutine onepin(ecm,id1in,id2in,id1out,id2out,misobr,mpart)
c ----------------------------------------------------------------------
c     program determines the mass and type of isobar(misobr and id1out)
c     and the mass and type of nucleon(mpart and id2out) in a nucleon-
c     nucleon collision between the two types of nucleons(id1in and
c     id2in) with total c.m. energy, ecm.  the choice of isobar type is
c     based on the assumption of the formation of a (3/2,3/2)isobar and
c     the conservation of isotopic spin in the collision.  the model is
c     patterned after that discussed by s.lindenbaum and r.sternheimer,
c     phys. rev. 105, p.1874(1957) for single pion production.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      real*8 id1in,id2in,id1out,id2out,misobr,mpart,massn,massp,mimin
     1,mimax,mmax,mn
      data massn,massp,mimin,mimax/939.516,938.223,1080.000,1573./
c ----------------------------------------------------------------------
c
      if(id1in.ne.id2in)go to 30
      if(id1in.eq.1.)go to 20
c
      random=unirn(dummy)
      eta=random*1.d0
      if(eta.le..75)go to 5
      id1out=24.
      id2out=2.
      mpart=massn
      go to 40
5     continue
c
      id1out=25.
      id2out=1.
      mpart=massp
      go to 40
20    continue
      random=unirn(dummy)
      eta=random*1.d0
      if(eta.le..75)go to 25
c
      id1out=14.
      id2out=1.
      mpart=massp
      go to 40
25    continue
c
      id1out=13.
      id2out=2.
      mpart=massn
      go to 40
30    continue
      random=unirn(dummy)
      eta=random*1.d0
      if(eta.le..5)go to 35
c
      id1out=24.
      id2out=1.
      mpart=massp
      go to 40
35    continue
c
      id1out=14.
      id2out=2.
      mpart=massn
40    continue
c
      mn=mpart
      mmax=ecm-mn
      if(mmax.gt.mimax)mmax=mimax
      rangem=mmax-mimin
      pmax=phase(ecm,mn,mimin)*200.
c
50    continue
      random=unirn(dummy)
      eta1=random*1.d0
      random=unirn(dummy)
      eta2=random*1.d0
      misobr=mimin+rangem*eta1
      tcm=misobr-1077.811
csk ---------------------------------------------------11/01/95 start --
      idum0=0
      idum1=1
      call ptotal(tcm,idum1,idum0,idum0,sig33,absor)
csk   call ptotal(tcm,1,0,0,sig33,absor)
csk ---------------------------------------------------11/01/95  end  --
      pm=sig33*phase(ecm,mn,misobr)/pmax
      if(eta2.gt.pm)go to 50
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine angli1(coscm)
c ----------------------------------------------------------------------
c     program selects the cos of the angle in the c.m.system of the
c     isobar created in a nucleon-nucleon collision which produces
c     only one pion or, equivalenty in the c.m.system it creates an
c     isobar and a nucleon.
c     this routine is subject to change because the actual distribution
c     is unknown accept that it is expected to be peaked forward and
c     backward.
c     start with p(coscm)=.25+.75*coscm*coscm
c     note: pmax=1.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c ----------------------------------------------------------------------
c
   10 continue
      random=unirn(dummy)
      eta1=random*1.d0
      random=unirn(dummy)
      eta2=random*1.d0
      coscm=2.*eta1-1.
      p=.25+.75*coscm*coscm
      if(eta2.gt.p)go to 10
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine branch(tcm,i9,i10,i11,finel,finex)
c ----------------------------------------------------------------------
c     program calculates the probability of an inelastic collision
c     (finel) and the probability of an inelastic or charge exchange
c     collision (finex).  program does linear interpolation on equally
c     spaced data.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension finii(36),finij(42),fini0(42),finxij(42),finxi0(42)
      data finii/0.0,.034,.037,.080,.125,.136,.143,.200,.211,.222,.235,
     1.294,.313,.313,.375,.375,.438,.471,.471,.556,.556,.526,.550,.524,
     2.545,.522,.500,.500,.480,.480,.500,.462,.462,.481,.481,.481/
      data finij/0.0,0.0,.037,.074,.111,.111,.148,.148,.143,.207,.200,
     1.194,.273,.265,.278,.316,.357,.400,.426,.438,.447,.444,.429,.450,
     2.447,.447,.474,.487,.500,.500,.489,.510,.481,.446,.431,.407,.414,
     3.421,.426,.420,.447,.467/
      data fini0/0.0,0.0,.014,.029,.047,.050,.069,.089,.089,.148,.167,
     1.167,.222,.241,.250,.286,.317,.371,.391,.406,.422,.419,.431,.448,
     2.446,.482,.500,.500,.517,.500,.500,.514,.487,.463,.440,.429,.441,
     3.429,.438,.436,.459,.472/
      data finxi0/.405,.410,.431,.441,.453,.483,.483,.518,.518,.555,
     1.574,.574,.592,.611,.607,.607,.617,.629,.641,.656,.641,.645,.672,
     2.655,.660,.696,.679,.672,.684,.656,.647,.653,.645,.638,.630,.643,
     3.631,.619,.613,.590,.594,.583/
      data finxij/.607,.593,.593,.593,.592,.592,.592,.592,.572,.586,
     1.567,.549,.576,.559,.556,.553,.571,.578,.596,.605,.596,.600,.596,
     2.600,.605,.605,.606,.615,.625,.619,.600,.612,.596,.571,.569,.560,
     3.552,.561,.556,.540,.553,.556/
c ----------------------------------------------------------------------
c
c for i-i events data are from tcm=310 mev to tcm=660 mev,for i-j and
c i-j0 events data are from tcm=250 mev to tcm=660 mev
      finel=0.0
      finex=0.0
      if(tcm.ge.660.)go to 1000
      if(i9.eq.1)go to 20
c have an i-j or i-j0 interaction
      if(tcm.lt.250.)return
      re=.1*tcm-25.
      i=idint(re)
      i=i+1
      de=re-float(i-1)
      ae=1.-de
      if(i10.eq.1)go to 10
c have i-j0 interaction
      finel=ae*fini0(i)+de*fini0(i+1)
      finex=ae*finxi0(i)+de*finxi0(i+1)
      return
10    continue
c have an i-j interaction
      finel=ae*finij(i)+de*finij(i+1)
      finex=ae*finxij(i)+de*finxij(i+1)
      return
20    continue
c have an i-i interaction
      if(tcm.lt.310.)return
      re=.1*tcm-31.
      i=idint(re)
      i=i+1
      de=re- float(i-1)
      ae=1.-de
      finel=ae*finii(i)+de*finii(i+1)
      finex=finel
      return
1000  continue
c%%      print 100,tcm,i9,i10,i11
c%100 format(1h1,'energy too large in branch.e=',f10.4,' i9,i10,and i11=
c%%     1',3i10)
      write(6,100) tcm,i9,i10,i11
  100 format(/' *** error message from s.branch ***'
     &/' given energy was too large in branch process.'
     &/' e =',1pe13.6,'  i9,i10,i11=',3i10)
      call parastop( 842 )
      end subroutine


************************************************************************
*                                                                      *
      subroutine pangle(e,ind9,ind10,ind11,ws3,ws29,ans)
c ----------------------------------------------------------------------
c     selects the pion or nucleon c.m.angle for isobar decay,etc. from
c     i-j, i-i, or pi0-ij differential elastic cross sections assuming
c     these cross sections are of the form-alpha+beta*x+gamma*x*x where
c     x is the cosine of the c.m.angle which is actually selected
c     e is the c.m.k.e., ans is the cosine of the angle selected, ind9=1
c     means an i-i decay, ind10=1 means an i-j decay, and ind11 means a
c     pi0-nucleon decay.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension scaiia(67),scaiib(67),scaija(67),scaijb(67),scaoa(67)
     1,scaob(67),scea(67),sceb(67)
      data scaiia/0.0,-.969,-1.62,-1.87,-1.85,-1.71,-1.63,-1.41,-1.22,
     1-1.00,-.809,-.583,-.389,-.260,-.150,-.045,.079,.297,.603,.923,
     21.18,1.43,1.74,2.00,2.33,2.64,2.88,3.07,3.22,3.38,3.46,3.43,3.58,
     33.56,3.53,3.56,3.60,3.64,3.69,3.83,4.00,
     x    4.0,4.13,4.29,4.47,4.41,4.38,4.33,4.29,
     14.07,3.70,3.85,3.85,3.85,3.70,3.70,3.57,3.45,3.67,4.00,4.19,4.19,
     24.38,4.69,4.69,4.69,5.16/
      data scaiib/0.0,.111,.380,.677,.922,1.10,1.27,1.39,1.49,1.65,1.80,
     12.08,2.15,2.18,2.26,2.41,2.55,2.69,2.84,3.07,3.20,3.18,3.35,3.48,
     23.67,3.92,4.18,4.40,4.52,4.63,4.59,4.52,4.63,4.50,4.35,4.25,4.07,
     33.93,3.85,3.75,3.82,
     x    3.80,3.91,4.05,4.21,4.41,4.38,4.33,4.64,
     14.44,4.44,4.62,4.62,4.62,4.44,4.44,4.64,4.48,4.67,5.00,5.16,5.48,
     25.63,6.25,6.56,6.88,7.10/
      data scaija/0.0,.414,.802,1.13,1.37,1.52,1.41,1.34,1.19,1.08,
     1.983,.896,.789,.714,.630,.577,.546,.527,.512,.481,.452,.441,.453,
     2.492,.559,.696,.796,.923,1.08,1.25,1.43,1.62,1.84,2.05,2.29,2.62,
     32.89,3.22,3.46,3.81,4.23,
     x    4.94,5.77,6.30,6.72,6.75,6.90,6.71,6.67,
     16.27,6.03,5.52,5.32,5.22,5.11,5.50,6.11,7.36,8.79,10.1,11.5,12.0,
     213.0,13.9,14.1,15.0,15.9/
      data scaijb/0.0,.048,.191,.411,.677,.957,1.14,1.29,1.42,1.61,1.86,
     12.18,2.50,2.88,3.27,3.32,3.16,2.73,2.28,2.00,1.85,1.72,1.64,1.61,
     21.59,1.68,1.78,1.92,2.14,2.42,2.72,3.02,3.47,3.90,4.39,5.05,5.68,
     36.38,6.97,7.75,8.71,
     x    9.71,10.3,10.9,11.3,11.3,11.3,11.0,11.0,
     110.3,10.7,11.0,12.1,13.7,14.8,17.3,19.2,22.4,26.7,30.4,33.5,34.3,
     235.9,37.1,36.7,37.0,37.4/
      data scaoa/0.0,-1.87,-1.33,-.936,-.688,-.521,-.762,-.652,-.576,
     1-.437,-.318,-.178,-.077,-.013,.028,.070,.139,.278,.538,.836,1.10,
     21.35,1.68,2.04,2.50,2.93,3.23,3.41,3.66,3.68,3.59,3.33,3.45,3.12,
     32.90,2.94,2.73,2.75,2.84,2.93,2.98,
     x    3.30,3.87,4.37,4.88,5.18,5.55,5.73,5.92,
     15.79,5.54,5.42,5.39,5.37,5.18,5.42,5.66,6.23,7.07,8.00,8.75,9.02,
     29.46,10.0,10.1,10.5,11.3/
      data scaob/0.0,.871,1.29,1.42,1.47,1.50,1.58,1.56,1.60,1.70,1.84,
     12.14,2.20,2.21,2.29,2.39,2.45,2.64,2.83,3.16,3.42,3.46,3.76,4.07,
     24.52,4.99,5.36,5.59,5.91,5.92,5.62,5.19,5.37,4.87,4.42,4.37,3.99,
     33.92,4.12,4.12,4.31,
     x    4.73,5.14,5.60,6.07,6.68,7.02,7.20,7.63,
     17.23,7.37,7.64,8.09,8.79,9.07,10.2,11.1,12.2,14.0,16.2,17.4,18.3,
     218.9,20.0,20.6,21.1,21.8/
      data scea/0.0,-.365,-.723,-1.05,-1.31,-1.50,-1.47,-1.48,-1.42,
     1-1.41,-1.38,-1.25,-1.06,-.857,-.556,-.222,.083,.486,.765,.939,
     21.03,1.13,1.25,1.26,1.31,1.42,1.50,1.60,1.68,1.88,2.07,2.31,2.50,
     33.00,3.63,4.00,5.60,6.75,6.75,8.67,13.0,
     x    12.4,11.5,10.4,9.2,7.6,6.0,4.6,3.6,2.6,2.0,
     11.6,1.2,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/
      data sceb/0.0,.036,.145,.323,.553,.813,.947,1.14,1.29,1.52,1.72,
     11.91,2.12,2.37,2.64,3.00,3.28,2.89,2.59,2.30,2.09,1.97,1.89,1.78,
     21.69,1.71,1.77,1.90,1.95,2.12,2.33,2.62,2.83,3.40,4.25,4.86,6.80,
     38.50,8.50,11.3,17.0,
     x    17.0,17.0,16.6,15.6,13.6,11.8,10.3,9.4,9.0,
     19.1,9.6,10.4,11.2,12.1,13.1,14.2,15.4,16.2,16.6,16.4,15.5,13.8,
     211.8,9.6,7.2,5.8/
c ----------------------------------------------------------------------
c
      rie=.1*e
      ie=rie
      ie=ie+1
      de=rie- float(ie-1)
      ae=1.-de
      if(ws3.ne.ws29)go to 30
      if (ind9-1)11,10,11
10    continue
      a=ae*scaiia(ie)+de*scaiia(ie+1)
      b=ae*scaiib(ie)+de*scaiib(ie+1)
      go to 14
   11 if (ind10-1) 12,13,12
13    continue
      a=ae*scaija(ie)+de*scaija(ie+1)
      b=ae*scaijb(ie)+de*scaijb(ie+1)
      go to 14
12    continue
      a=ae*scaoa(ie)+de*scaoa(ie+1)
      b=ae*scaob(ie)+de*scaob(ie+1)
14    continue
c find max of 1.+a*x+b*x*x,i.e.-am1
      am1=1.+a+b
      if(a.lt.0.0)am1=1.-a+b
      if(b.ge.0.0)go to 21
      am=-a/(2.*b)
      if(dabs(am).gt.1.)go to 21
      am2=1.+am*(a+b*am)
      if(am1.lt.am2)am1=am2
c use von-neuman rejection to find ans=cos(theta-c.m.)
c  21 irand=irand*513+300000001
  21  random=unirn(dummy)
c 21  randh =ransu(0)
      randh=random*1.d0
      ran1=randh+randh-1.
c if the probability distribution is negative for some value of costheta
c reject that value of costheta
      y12=1.+ran1*(a+b*ran1)
      if(y12.lt.0.0)go to 21
c     irand=irand*513+300000001
      random=unirn(dummy)
c     ran2  =ransu(0)
      ran2=random*1.d0
      y13=am1*ran2
      if(y12-y13)21,23,23
   23 ans=ran1
24    continue
c shut off ind(9),ind(10),ind(11).change made 10/5/70
      ind9=0
      ind10=0
      ind11=0
      return
30    continue
c charge exchange angular dist.
      a=ae*scea(ie)+de*scea(ie+1)
      b=ae*sceb(ie)+de*sceb(ie+1)
      go to 14
      end subroutine


************************************************************************
*                                                                      *
      subroutine pin1pi(ecm,id1in,id2in,id1out,id2out,misobr,mpart,df)
c ----------------------------------------------------------------------
c     program determines the mass and type of isobar(misobr and id1out)
c     and the mass and type of pion(mpart and id2out)in a pion-nucleon
c     collision between a pion of type, id1in, and a nucleon of type,
c     id2in, with total c.m. energy, ecm.  the choice of isobar type is
c     based on the assumption of the formation of a (3/2,3/2)isobar and
c     the conservation of isotopic spin in the collision.  the model is
c     patterned after that discussed by s.lindenbaum and r.sternheimer,
c     phys. rev. 109, p.1723(1958) for single pion production.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      real*8 id1in,id2in,id1out,id2out,misobr,mpart,massn,massp,mimin
     1,mimax,mmax,mn
     2,mpi0,mpichg
      data mpi0,mpichg/135.0,139.588/
      data massn,massp,mimin,mimax/939.516,938.223,1080.000,1573./
c ----------------------------------------------------------------------
c
c first decide what type of collision is occurring and then decide
c what type isobar and pion result from the collision
      df=0.0
      if(id1in.eq.3.)go to 10
      if(id1in.eq.4.)go to 20
      if(id1in.eq.5.)go to 30
c%%      print 100,id1in
c%%  100 format('0the following pi code is not valid',f10.5)
      write(6,100) id1in
  100 format(/' *** error message from s.pin1pi ***'
     &/' this routine was not accepted the type of isobar.'
     &/' id1in =',1pe13.6)
      call parastop( 824 )
10    continue
      if(id2in.eq.2.)go to 15
c this is a pi+proton collision
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta   =ransu(0)
      eta=random*1.d0
      if(eta.gt..6)go to 11
c an i++ and a pi0 have been created
      id1out=13.
      id2out=4.
      mpart=mpi0
      go to 40
11    continue
c an i+ and a pi+ have been created
      id1out=14.
      id2out=3.
      mpart=mpichg
      go to 40
15    continue
c this is a pi+ neutron collision
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta   =ransu(0)
      eta=random*1.d0
      if(eta.gt..2701393886)go to 16
c an i0 and a pi+ have been created
      id1out=24.
      id2out=3.
      mpart=mpichg
      go to 40
16    continue
      if(eta.gt..6222428337)go to 17
c an i+ and a pi0 have been created
      id1out=14.
      id2out=4.
      mpart=mpi0
      df=-1.
      go to 40
17    continue
c an i++ and and a pi- have been created
      id1out=13.
      id2out=5.
      mpart=mpichg
      df=-1.
      go to 40
20    continue
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta   =ransu(0)
      eta=random*1.d0
      if(eta.gt..6643997206)go to 25
      if(id2in.eq.2.)go to 22
c this is a pi0-proton collision giving an i++ and a pi-
      id1out=13.
      id2out=5.
      mpart=mpichg
      go to 40
22    continue
c this is a pi0-neutron collision giving an i- and a pi+
      id1out=25.
      id2out=3.
      mpart=mpichg
      go to 40
25    continue
      if(eta.gt..8708624797)go to 27
      if(id2in.eq.2.)go to 26
c this is a pi0-proton collision giving an i+ and a pi0
      id1out=14.
      id2out=4.
      mpart=mpi0
      go to 40
26    continue
c this is a pi0-neutron collision giving an i0 and a pi0
      id1out=24.
      id2out=4.
      mpart=mpi0
      go to 40
27    continue
      if(id2in.eq.2.)go to 28
c this is a pi0-proton collision giving an i0 and a pi+
      id1out=24.
      id2out=3.
      mpart=mpichg
      df=1.
      go to 40
28    continue
c this is a pi0-neutron collision giving an i+ and a pi-
      id1out=14.
      id2out=5.
      mpart=mpichg
      df=-1.
      go to 40
30    continue
      if(id2in.eq.1.)go to 33
c this is a pi- neutron collision
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta   =ransu(0)
      eta=random*1.d0
      if(eta.gt..6)go to 32
c an i- and a pi0 have been created
      id1out=25.
      id2out=4.
      mpart=mpi0
      go to 40
32    continue
c an i0 and a pi- have been created
      id1out=24.
      id2out=5.
      mpart=mpichg
      go to 40
33    continue
c this is a pi- proton collision
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta   =ransu(0)
      eta=random*1.d0
      if(eta.gt..2701393886)go to 34
c an i+ and a pi- have been created
      id1out=14.
      id2out=5.
      mpart=mpichg
      go to 40
34    continue
      if(eta.gt..6222428337)go to 35
c an i0 and a pi0 have been created
      id1out=24.
      id2out=4.
      mpart=mpi0
      df=1.
      go to 40
35    continue
c an i- and a pi+ have been created
      id1out=25.
      id2out=3.
      mpart=mpichg
      df=1.
      go to 40
40    continue
c now select the mass of the isobar
      mn=mpart
      mmax=ecm-mn
      if(mmax.gt.mimax)mmax=mimax
      rangem=mmax-mimin
      pmax=phase(ecm,mn,mimin)*200.
c the max of p(mi) is less than or equal (max.sig33)*(max.phase)
c for a given ecm.sig33 is the pi+ + proton cross section.
c now select the isobar mass by the von neuman rejection method.
50    continue
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta1  =ransu(0)
      eta1=random*1.d0
c     irand=irand*513+300000001
      random=unirn(dummy)
c     eta2  =ransu(0)
      eta2=random*1.d0
      misobr=mimin+rangem*eta1
      tcm=misobr-1077.811
csk ---------------------------------------------------11/01/95 start --
      idum0=0
      idum1=1
csk   call ptotal(tcm,1,0,0,sig33,absor)
      call ptotal(tcm,idum1,idum0,idum0,sig33,absor)
csk ---------------------------------------------------11/01/95  end  --
      pm=sig33*phase(ecm,mn,misobr)/pmax
      if(eta2.gt.pm)go to 50
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine find(arg,table,ii,ik)
c ----------------------------------------------------------------------
c     subroutine find is a simple table look up routine.  ik is the
c     number of entries in the table and ii is the desired entry.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension table(40)
c ----------------------------------------------------------------------
c
      do 1 i=1,ik
      if(arg-table(i))2,2,1
    1 continue
    2 ii=i
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine binter(e,ener,i,aint3,aint4,aint5)
c ----------------------------------------------------------------------
c     subroutine binter calculates the coefficients for triple
c     interpolation.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension ener(40)
c ----------------------------------------------------------------------
c
      aint0=ener(i-1)-ener(i)
      aint0=(e-ener(i+1))/aint0
      aint1=ener(i-1)-ener(i+1)
      aint1=(e-ener(i))/aint1
      aint2=ener(i)-ener(i+1)
      aint2=(e-ener(i-1))/aint2
      aint5=aint2*aint1
      aint3=aint0*aint1
      aint4=-aint2*aint0
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine aiouta(y)
c ----------------------------------------------------------------------
c     subroutine aiouta stores the decay particles from the escaping
c     isobars into st.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      common/hiro/st(8000)
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
      dimension y(30)
c ----------------------------------------------------------------------
c
      y(1)=-y(1)
      y(2)=-y(2)
      icsb(15)=icsb(15)+1
      icsb8=icsb(5)*icsb(15)
      st(icsb8+1)=y(1)
      st(icsb8+2)=y(2)
      st(icsb8+3)=y(3)
      st(icsb8+4)=y(4)
      st(icsb8+5)=y(5)
      st(icsb8+6)=y(6)
      st(icsb8+7)=y(7)
      st(icsb8+8)=y(8)
      st(icsb8+9)=y(9)
      st(icsb8+10)=y(10)
      st(icsb8+18)=y(10)
      st(icsb8+11)=y(11)
      st(icsb8+19)=y(11)
      st(icsb8+12)=y(12)
      st(icsb8+20)=y(12)
      st(icsb8+21)=y(5)
      st(icsb8+13)=0.0
      st(icsb8+14)=y(14)
      st(icsb8+15)=y(24)
      st(icsb8+16)=y(25)
      st(icsb8+17)=y(27)
      return
      end subroutine

************************************************************************
*                                                                      *
      function ciicug(e)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension ebin(27),pp(27)
      data pmas / 0.938259/
      data ebin / 40., 50., 60., 70., 80., 90., 100., 110., 125., 150.,
     & 175., 200., 225., 250., 275., 300., 350., 400., 450., 500., 600.,
     & 700., 800.,1000., 1200., 1400., 1500./
      data pp /  55., 43.6, 37.7, 34.2, 31.8, 30.2, 28.9, 27.9, 26.8,
     & 25.5, 24.6, 24.0, 23.5, 23.2, 22.9, 22.6, 22.2, 21.9, 21.7, 21.6,
     & 21.3, 21.1, 21.0, 20.8, 20.3, 19.8, 19.6/
c ----------------------------------------------------------------------
c
      celas = 0.
      cinel = 0.
c
      if (e.le.40.) then
c     ciicug = pp(1)
      egev = e*1.e-3
      plab = dsqrt(egev**2 + 2.*egev*pmas  )
      celas  = 23.5+1000*(0.7-plab)**4
c     ciicug = pp(1)
      go to 999
      endif
      if (e.gt.1500.) then
      celas  = pp(27)
      go to 999
      endif
         do 100 i=1,26
         if (e.gt.ebin(i) .and. e.le.ebin(i+1)) then
             m = i
             go to 200
         endif
  100 continue
c%%    write(6,*) ' **** invarid energy in ciicug *****'
      write(6,61) e,ebin(1),ebin(27)
   61 format(/' *** error message from s.ciicug ***'
     &/' given energy was not included the energy table range in ',
     & 'ciicug.'
     &/' e =',1pe13.6,'  energy range =',2e13.6)
      call parastop( 840 )
  200 celas =(pp(m+1)-pp(m))/(ebin(m+1)-ebin(m))*(e-ebin(m)) + pp(m)
c
c inelastic pp cross sections for e > 350 mev
c   see nucl. phys. a489 (1988) p781-802 (j.cugnon and m.-c. remaire )
c
csk ---------------------------------------------------10/17/95 start --
      if (e.gt.350) then
      egev = e*1.e-3
      plab = dsqrt(egev**2 + 2.*pmas*egev)
c
csk     if (plab.le.1.5) cinel=23.5 + 24.6/(1.+dexp(-10.*plab+12.)
csk  &                        -1250./(plab+50.) + 4.*(plab-1.3)**2
        if (plab.le.1.5) cinel=23.5 + 24.6/(1.+dexp(-10.*plab+12.))
     &                        -1250./(plab+50.) + 4.*(plab-1.3)**2
c
        if (plab.gt.1.5 .and. plab.le.2.0)
     &  cinel=41.0 + 60.0*(plab-0.9)*dexp(-1.2*plab)
     &       -1250./(plab+50.) + 4.*(plab-1.3)**2
c
        if (plab.gt.2.0 )
     &  cinel=41.0 + 60.0*(plab-0.9)*dexp(-1.2*plab)
     &       -77.0/(plab+1.5)
c
      endif
csk ---------------------------------------------------10/17/95  end  --
  999   ciicug = celas + cinel
      return
      end function

************************************************************************
*                                                                      *
      function cijcug(e)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension ebin(27),pn(27)
      data pmas / 0.938259/
      data ebin / 40., 50., 60., 70., 80., 90., 100., 110., 125., 150.,
     & 175., 200., 225., 250., 275., 300., 350., 400., 450., 500., 600.,
     & 700., 800.,1000., 1200., 1400., 1500./
      data pn /  55., 44.9, 39.7, 36.6, 34.5, 33.0, 31.9, 31.0, 30.0,
     & 28.9, 28.1, 27.5, 27.1, 26.8, 26.5, 26.3, 26.0, 25.7, 25.5, 25.4,
     & 25.2, 25.0, 24.9, 24.7, 24.1, 23.4, 23.1/
c ----------------------------------------------------------------------
c
      celas = 0.
      cinel = 0.
c
      if (e.le.0.0d0) goto 999
      if (e.le.40.) then
c     cijcug = pp(1)
      egev = e*1.e-3
      plab = dsqrt(egev**2 + 2.*egev*pmas)
      celas  = 23.5+1000*(0.7-plab)**4
c     cijcug = pp(1)
      go to 999
      endif
      if (e.gt.1500.) then
      celas  = pn(27)
      go to 999
      endif
      do 100 i=1,26
      if (e.gt.ebin(i) .and. e.le.ebin(i+1)) then
         m = i
         go to 200
      endif
  100 continue
c%%    write(6,*) ' **** invarid energy in cijcug *****'
      write(6,61) e,ebin(1),ebin(27)
   61 format(/' *** error message from s.cijcug ***'
     &/' given energy was not included the energy table range in ',
     & 'cijcug.'
     &/' e =',1pe13.6,'  energy range =',2e13.6)
      call parastop( 839 )
  200 celas =(pn(m+1)-pn(m))/(ebin(m+1)-ebin(m))*(e-ebin(m)) + pn(m)
c
c inelastic pp cross sections for e > 350 mev
c   see nucl. phys. a489 (1988) p781-802 (j.cugnon and m.-c. remaire )
c
csk ---------------------------------------------------10/17/95 start --
      if (e.gt.350) then
      egev = e*1.e-3
      plab = dsqrt(egev**2 + 2.*pmas*egev)
c
        if (plab.le.1.0)
     & cinel=33.0 + 196.*dsqrt((abs(plab-0.95))**5)
     &      -31.1 /(dsqrt(plab))
c
        if (plab.gt.1.0 .and. plab.le.2.0)
     & cinel=24.2 + 8.9 *plab
     &      -31.1 /(dsqrt(plab))
c
        if (plab.gt.2.0 )
     & cinel=42.0 + 77.0/(plab+1.5)
c
      endif
csk ---------------------------------------------------10/17/95  end  --
  999   cijcug = celas + cinel
      return
      end function


************************************************************************
*                                                                      *
      function cubert(x)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c ----------------------------------------------------------------------
c
      p=1./3.
      cubert=x**p
      return
      end function


************************************************************************
*                                                                      *
      function en(dmass,pmass,ws4,beta)
c ----------------------------------------------------------------------
c     function en is called by aitota to check if the isobar has enough
c     energy to decay, that is to have enough energy to stay above the
c     fermi level.
c     assume parallel decay because this always gives more energy than
c     when the particles come off opposite each other.  for more details
c     see notes and jackson book on classical electrodynamics,page 395.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
c dmass and pmass are the masses of the decaying nucleon and pion.
      dmass2=dmass*dmass
      pmass2=pmass*pmass
c ws4 is the mass of the isobar.
      amass2=ws4*ws4
c betapn is the beta of the nucleon in the center of mass when the
c  isobar is at rest.
c then do lorentz transformation (see jackson p.361 and notes) to bring
c  particle into the moving frame of isobar(or the lab frame).
      betapn=amass2+dmass2-pmass2
      betapn=(betapn*betapn-4.*amass2*dmass2)/(betapn*betapn)
      if(betapn.lt.0.0)betapn=0.0
      betapn=dsqrt(betapn)
c betan and en are the beta and energy of the decay nucleon in the lab
c  frame.
      betan=(betapn+beta)/(1.0+betapn*beta)
      en=dmass/(dsqrt(1.0-betan*betan))
      return
      end function


************************************************************************
*                                                                      *
      function finlii(e)
c ----------------------------------------------------------------------
c     program calculates the ratio of the cross section for producing
c     one pion to the total cross section in i-i(p-p) collisions
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
      if(e.gt.350.)go to 10
      finlii=0.0
      go to 1000
10    continue
      if(e.gt.400.)go to 20
      finlii=-.14+.0004*e
      go to 1000
20    continue
      if(e.gt.550.)go to 30
      finlii=-.3+.0008*e
      go to 1000
30    continue
      if(e.gt.700.)go to 40
      finlii=-.74+.0016*e
      go to 1000
40    continue
      if(e.gt.800.)go to 50
      finlii=-.04+.0006*e
      go to 1000
50    continue
      if(e.gt.950.)go to 60
      finlii=.28+.0002*e
      go to 1000
60    continue
      finlii=.4852-.000016*e
1000  return
      end function


************************************************************************
*                                                                      *
      function finlij(e)
c ----------------------------------------------------------------------
c     program calculates the ratio of the cross section for producing
c     one pion to the total cross section in i-j(p-n) collisions.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
      if(e.gt.300.)go to 10
      finlij=0.0
      go to 1000
10    continue
      if(e.gt.400.)go to 20
      finlij=-.09+.0003*e
      go to 1000
20    continue
      if(e.gt.500.)go to 30
      finlij=-.25+.0007*e
      go to 1000
30    continue
      if(e.gt.650.)go to 40
      finlij=-.63333+.0014667*e
      go to 1000
40    continue
      if(e.gt.700.)go to 50
      finlij=.060+.0004*e
      go to 1000
50    continue
      if(e.gt.800.)go to 60
      finlij=.27+.0001*e
      go to 1000
60    continue
      finlij=.39571-.000057143*e
1000  return
      end function


************************************************************************
*                                                                      *
      function ifrist(type,xn5,iloc,dens,gws,ws19,st)
c ----------------------------------------------------------------------
c     calculates whether a collision with a nucleon is forbidden now
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/switch/reat,vpion,isonsw,idistr,icb,nzman
      dimension gws(15),st(8000)
c ----------------------------------------------------------------------
c
      ifrist=0
      if(idistr.le.0)return
      itloc=dabs(xn5)
      itloc=itloc+iloc+12
      itype=type
      iprev=st(itloc)/100.
      if(iprev.ne.itype) return
      ws15=st(itloc)-iprev*100
      if(ws15.le.0.) return
      delt=ws19-ws15
c that could be the incident particle
      d=0.0
      do 10 i=1,3
      k=itloc+i
      j=9+i
      d=d+(st(k)-gws(j))*(st(k)-gws(j))
10    continue
      d=dsqrt(d)
      go to (20,30,40),idistr
20    continue
c the restriction is .5 f
creat if(delt.ge.reat) return
      if(d.lt..5)ifrist=2
      return
30    continue
c the restriction is v*dens=1
      dcubed=1./(4.18879020*dens)
c fermi velocity is .258 . (10*.258)**3=17.172
      if(delt*delt*delt.gt.(dcubed/17.172)) return
      if(d*d*d.lt.dcubed)ifrist=2
      return
40    continue
c the restriction is 1/cubert(density)
      dcubed=1./dens
      if(delt*delt*delt.gt.(dcubed/17.172)) return
      if(d*d*d.lt.dcubed) ifrist=2
      return
      end function


************************************************************************
*                                                                      *
      function phase(ecm,m3,mi)
c ----------------------------------------------------------------------
c     calculates the two body phase factor
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      real*8 m3,mi
c ----------------------------------------------------------------------
c
      e3=(ecm*ecm+m3*m3-mi*mi)/(2.*ecm)
      p3sq=e3*e3-m3*m3
      e4=ecm-e3
      if((e3.le.m3).or.(e4.le.mi))go to 10
      phase=dsqrt(p3sq)*e3*e4/ecm
      return
10    continue
      phase=0.
      return
      end function
