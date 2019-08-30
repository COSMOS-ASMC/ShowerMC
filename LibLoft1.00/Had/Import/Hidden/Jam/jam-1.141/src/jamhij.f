c***********************************************************************
c***********************************************************************
c                                                                      *
c        PART 5: Calculate Jet cross section from HIJING               *
c                                                                      *
c   List of subprograms in rough order of relevance with main purpose  *
c      (S = subroutine, F = function, B = block data, E = entry)       *
c                                                                      *
c                                                                      *
c b hidata   to give default values to switches and parameters HIJING  *
c s hijcrs   to calculate the hadronic cross sections                  *
c f ftot     to give the function for cal. of total cross section      *
c f fhin     to give the function for cal  of inel. cross section      *
c f ftotjet  to give the function for cal. of avarage jet x section    *
c f ftotrig  to give the function for cal. fo trigged jet x section    *
c f fnjet    to give the function for the prob. of # of hard scatt.    *
c f sgmin    to calculate log(n!)                                      *
c f omg0     to give the eikonal function of nucleon                   *
c f romg     to give the eikonal function from a table                 *
c f bk       to give BK=EXP(-X)*(X**2-X4**2)**2.50/15.0                *
c                                                                      *
c s crsjet   to calculate the jet cross section by using VEGAS         *
c f fjet     to provide func. for jet x-section                        *
c f fjetrig  to provide func. for trigged jet x-section                *
c f ghvq                                                               *
c f gphoton  to provide function for the calculation of direct photon  *
c f g        to provide function for the calculation of QCD processes  *
c f subcrs1  to give function for q q'-> qq' t-channel gluon exchange  *
c f subcrs2  to give function for q q'-> qq' s-channel gluon exchange  *
c f subcrs3  to give function for qq~->qq~  same flavour               *
c f subcrs4  to give function for qq~->gg                              *
c f subcrs5  to give function for gg->qq~                              *
c f subcrs6  to give function for gq->gq                               *
c f subcrs7  to give function for gg->gg                               *
c s parton   to calculate the structure function with nuclear effect   *
c f gmre     to calculate log. of Gamma function                       *
c                                                                      *
c                    Soft part                                         *
c                                                                      *
c f fnkick   to define function for pt**2 kick for soft interaction    *
c f fnkick2  to define the expotnetial pt distribution of sea quarks   *
c f fnstru   to define func. of x dist. of valence quarks for baryons  *
c f fnstrum  to define func. of x dist. of valence quarks for mesons   *
c f fnstrus  to define func. of x dist. of single diffractive scatt    *
c                                                                      *
c s hijwds   to set up histogram IDH by Wood Saxon                     *
c f wdsax    to give three parameter wood saxon                        *
c f rwdsax   to give r**2*woodaxon                                     *
c f wdsax1   to give three parameter Wood Saxon for projectile         *
c f wdsax2   to give three parameter Wood Saxon for target             *
c                                                                      *
c s hifun    to setup probability for any function                     *
c f hirnd    to give a random number according to the given function   *
c f hirnd2   to generate random number between xmin and xmax           *
c s jamfun   to setup probability distribution for any function        *
c f jamrnd2  to generate random number between xmin and xmax           *
c                                                                      *
c            General utilities                                         *
c                                                                      *
c s vegas0   to perform N-dimensional Monte Calro integ'n              *
c e vegas1   to initializes cummulative variables, but not grid        *
c e vegas2   to integrate with no initialization                       *
c e vegas3   to do main integration loop                               *
c s aran9    to generate random number for vegas                       *
c f gauss1   to perform gaussian one-dimensional integration           *
c f gauss2   to perform gaussian one-dimensional integration           *
c                                                                      *
c***********************************************************************
c***********************************************************************

      block data hidata

c...Purpose: to give default values to switches and parameters
c...for HIJING.
      implicit double precision(a-h, o-z)
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn,iprn,ierr
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/hijdat/hidat0(10,10),hidat(10)
      common/hipyint/mint4,mint5,atco(200,20),atxs(0:200)
      save  /bveg1/,/hiparnt/,/hijdat/,/hipyint/

c...Give all the switchs and parameters the default values.
        data hipr1/
     &  1.5d0,  0.35d0, 0.5d0,  0.9d0,  2.0d0,  0.1d0,  1.5d0,  2.0d0, 
     &  -1.0d0, -2.25d0,
     &  2.0d0,  0.5d0,  1.0d0,  2.0d0,  0.2d0,  2.0d0,  2.5d0,  0.3d0, 
     &   0.1d0,  1.4d0,
     &  1.6d0,  1.0d0,  2.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0, 
     &   0.4d0,  57.0d0,
     &  28.5d0, 3.9d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0, 
     &   0.0d0,  
     &  3.141592654d0,
     &  0.0d0,  0.4d0,  0.1d0,  1.5d0,  0.1d0, 0.25d0, 0.0d0,  0.5d0, 
     &   0.0d0,  0.0d0,
     &  50*0.0d0/

        data ihpr2/
     &  1,    3,    0,    1,    1,    1,    1,    10,    0,    0,
     &  1,    1,    1,    1,    0,    0,    1,     0,    0,    1,
     &  30*0/

c...Initialize all the data common blocks.
        data hint1/100*0/
        data ihnt2/50*0/

        data mint4/0/mint5/0/atco/4000*0.0d0/atxs/201*0.0d0/
        data (hidat0(1,i),i=1,10)/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, 
     & 2.25d0,
     &          2.5d0,4.0d0,4.1d0/
        data (hidat0(2,i),i=1,10)/2.0d0,3.0d0,5.0d0,6.0d0,7.0d0,8.0d0, 
     & 8.0d0,10.0d0,
     &          10.0d0,10.0d0/
        data (hidat0(3,i),i=1,10)/1.0d0,0.8d0,0.8d0,0.7d0,0.45d0, 
     & 0.215d0,
     &          0.21d0,0.19d0,0.19d0,0.19d0/
        data (hidat0(4,i),i=1,10)/0.35d0,0.35d0,0.3d0,0.3d0,0.3d0,0.3d0,
     &          0.5d0,0.6d0,0.6d0,0.6d0/
        data (hidat0(5,i),i=1,10)/23.8d0,24.0d0,26.0d0,26.2d0,27.0d0, 
     & 28.5d0,28.5d0,
     &          28.5d0,28.5d0,28.5d0/
        data ((hidat0(j,i),i=1,10),j=6,9)/40*0.0d0/
        data (hidat0(10,i),i=1,10)/5.0d0,20.0d0,53.0d0,62.0d0,100.0d0, 
     & 200.0d0,
     &          546.0d0,900.0d0,1800.0d0,4000.0d0/
        data hidat/10*0.0d0/

c       data num1/30123984/
        data xl/10*0.d0/,xu/10*1.d0/
        data ncall/1000/itmx/100/acc/0.01d0/nprn/0/iprn/8/ierr/0/

        end

c***********************************************************************

      subroutine hijcrs

c...This is to calculate the cross sections of jet production and
c...the total inelastic cross sections.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/njet/n,ip_crs
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      save  /hiparnt/,/njet/,/pjpars/

      external fhin,ftot,fnjet,ftotjet,ftotrig

c...Calculate jet cross section(in mb).
      if(hint1(1).ge.parp(2)) then
        call crsjet
      endif
 
      aphx1=hipr1(6)*(ihnt2(1)**0.3333333d0-1.0d0)
      aphx2=hipr1(6)*(ihnt2(3)**0.3333333d0-1.0d0)

c....The averaged inclusive cross section sigma_jet
      hint1(11)=hint1(14)-aphx1*hint1(15)
     &                  -aphx2*hint1(16)+aphx1*aphx2*hint1(17)

c***Total and Inel cross section are calculated by Gaussian integration.

c...The averaged cross section for jet production.
      hint1(10)=gauss1(ftotjet,0.0d0,20.0d0,0.01d0)

c...the averaged inelastic cross section.
      hint1(12)=gauss1(fhin,0.0d0,20.0d0,0.01d0)

c...total
      hint1(13)=gauss1(ftot,0.0d0,20.0d0,0.01d0)

c...triggered jet.
      hint1(60)=hint1(61)-aphx1*hint1(62)
     &                  -aphx2*hint1(63)+aphx1*aphx2*hint1(64)

      hint1(59)=gauss1(ftotrig,0.0d0,20.0d0,0.01d0)
      if(hint1(59).eq.0.0d0) hint1(59)=hint1(60)

c...The probability of j=0-20 number of hard scatterings per n-n.
      if(hint1(1).ge.10.0d0) then
        do 20 i=0,20
          n=i
          hint1(80+i)=gauss1(fnjet,0.0d0,20.0d0,0.01d0)/hint1(12)
 20     continue
      endif

      hint1(10)=hint1(10)*hipr1(31) ! jet x-section (mb)
      hint1(12)=hint1(12)*hipr1(31) ! inel x-section (mb)
      hint1(13)=hint1(13)*hipr1(31) ! total x-section (mb)
      hint1(59)=hint1(59)*hipr1(31) ! trigged jet x-section (mb)

c***Parametrized cross section for single diffractive
c***reaction(Goulianos)(not used now).
      if(ihpr2(13).ne.0) then
        hipr1(33)=1.36d0*(1.0d0+36.0d0/hint1(1)**2)
     &             *dlog(0.6d0+0.1d0*hint1(1)**2)
        hipr1(33)=hipr1(33)/hint1(12)
      endif

      end

c***********************************************************************

      double precision function ftot(x)

c...Set function to calculate hadronic total cross section.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/

      omg=omg0(x)*(hipr1(30)+hint1(11))/hipr1(31)/2.0d0
      ftot=2.0d0*(1.0d0-exp(-omg))

      return
      end
 
c***********************************************************************
 
      function fhin(x)

c...Set function for the calculation of inelastic cross section.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/
      omg=omg0(x)*(hipr1(30)+hint1(11))/hipr1(31)/2.0d0
      fhin=1.0d0-exp(-2.0d0*omg)

      end
 
c***********************************************************************
 
      function ftotjet(x)

c...Set function for the calculation of Jet cross section.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/

      omg=omg0(x)*hint1(11)/hipr1(31)/2.0d0
      ftotjet=1.0d0-exp(-2.0d0*omg)

      end

c***********************************************************************

      function ftotrig(x)

c...Set function for the calculation of trigged jet cross section.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/

      omg=omg0(x)*hint1(60)/hipr1(31)/2.0d0
      ftotrig=1.0d0-exp(-2.0d0*omg)

      end

c***********************************************************************

      function fnjet(x)

c...Function for the probability of number of hard scatterings per n-n.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/
      common/njet/n,ip_crs
      save  /njet/

      omg1=omg0(x)*hint1(11)/hipr1(31)
      c0=exp(n*dlog(omg1)-sgmin(n+1))
      if(n.eq.0) c0=1.0d0-exp(-2.0d0*omg0(x)*hipr1(30)/hipr1(31)/2.0d0)
      fnjet=c0*exp(-min(50.0d0,omg1))

      end

c***********************************************************************

      function sgmin(n)

c...log(n!).
      implicit double precision(a-h, o-z)

      ga=0.d0
      if(n.le.2) go to 20
      do 10 i=1,n-1
      z=i
      ga=ga+dlog(z)
10    continue
20    sgmin=ga

      end

c***********************************************************************

      function omg0(x)

c...Eikonal function of nucleon.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common /besel/x4
      save  /hiparnt/,/besel/
      external bk

      x4=hipr1(32)*sqrt(x)
      omg0=hipr1(32)**2*gauss2(bk,x4,x4+20.0d0,0.01d0)/96.0d0

      end

c***********************************************************************

      function romg(x)

c...Gives the eikonal function from a table calculated in the irst call.
      implicit double precision(a-h, o-z)
      dimension fr(0:1000)
      data i0/0/
      save fr,i0

      if(i0.ne.0) go to 100
      do 50 i=1,1001
        xr=(i-1)*0.01d0
        fr(i-1)=omg0(xr)
50    continue

100   i0=1
      if(x.ge.10.0d0) then
         romg=0.0d0
         return
      endif
      ix=int(x*100)
      romg=( fr(ix)*((ix+1)*0.01d0-x) + fr(ix+1)*( x-ix*0.01d0 ))/0.01d0

      end

c***********************************************************************

      function bk(x)

      implicit double precision(a-h, o-z)
      common /besel/x4
      save   /besel/

      bk=exp(-x)*(x**2-x4**2)**2.50d0/15.0d0

      end

c***********************************************************************

      subroutine crsjet

c...Purpose: to calculate the jet cross section,
c...the integration is done by using VEGAS.

      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/njet/n,ip_crs
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn,iprn,ierr
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      common/bveg3/f,ti,tsi
      save  /hiparnt/,/njet/,/bveg1/,/bveg2/,/bveg3/
      external fjet,fjetrig
      parameter (eps=1d-2)

c...NCALL give the number of inner-iteration, ITMX 
c...gives the limit of out-iteration. Nprn is an option
c...(1: print the integration process. 0: do not print)

      ndim=3
      hint1(14)=0.0d0
      hint1(15)=0.0d0
      hint1(16)=0.0d0
      hint1(17)=0.0d0

c...Total inclusive jet cross section (pT>p0) 
      ip_crs=0
      call vegas0(fjet,avgi,sd,chi2a)
      if(ierr.ne.0) then
        write(3,*)'vega=',hint1(1),ihnt2(5),ihnt2(6),avgi/2.5682d0
      endif

ccc   call glp(ndim,fjet,eps,avgi2,esterr,nt)

c...Jet cross section without nuclear shadowing effect.
      hint1(14)=avgi/2.5682d0

c...Jet cross section with nuclear shaowing effect.
      if(ihpr2(6).eq.1) then

c...Jet cross section to account for the projectile shadowing effect.
          if(ihnt2(1).gt.1) then
            ip_crs=1
            call vegas0(fjet,avgi,sd,chi2a)
            hint1(15)=avgi/2.5682d0
          endif

c...Jet cross section to account for the target shadowing effect.
          if(ihnt2(3).gt.1) then
            ip_crs=2
            call vegas0(fjet,avgi,sd,chi2a)
            hint1(16)=avgi/2.5682d0
          endif

c...Jet cross section to account for the cross term of shadowing effect.
          if(ihnt2(1).gt.1.and.ihnt2(3).gt.1) then
            ip_crs=3
            call vegas0(fjet,avgi,sd,chi2a)
            hint1(17)=avgi/2.5682d0
          endif

      endif

c...Cross section of trigger jet.
      if(ihpr2(3).ne.0) then
        ip_crs=0
        call vegas0(fjetrig,avgi,sd,chi2a)
        hint1(61)=avgi/2.5682d0
        if(ihpr2(6).eq.1) then
          if(ihnt2(1).gt.1) then
             ip_crs=1
             call vegas0(fjetrig,avgi,sd,chi2a)
             hint1(62)=avgi/2.5682d0
          endif
          if(ihnt2(3).gt.1) then
             ip_crs=2
             call vegas0(fjetrig,avgi,sd,chi2a)
             hint1(63)=avgi/2.5682d0
          endif
          if(ihnt2(1).gt.1.and.ihnt2(3).gt.1) then
             ip_crs=3
             call vegas0(fjetrig,avgi,sd,chi2a)
             hint1(64)=avgi/2.5682d0
          endif
        endif
      endif

      end

c***********************************************************************

        function fjet(x,wgt)

c...Provide function for the calculation of jet cross section.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/
      dimension x(10)

      pt2=(hint1(1)**2/4.0d0-hipr1(8)**2)*x(1)+hipr1(8)**2
      xt=2.0d0*sqrt(pt2)/hint1(1)
      ymx1=dlog(1.0d0/xt+sqrt(1.0d0/xt**2-1.0d0))
      y1=2.0d0*ymx1*x(2)-ymx1
      ymx2=dlog(2.0d0/xt-dexp(y1))
      ymn2=dlog(2.0d0/xt-dexp(-y1))
      y2=(ymx2+ymn2)*x(3)-ymn2
      fjet=2.0d0*ymx1*(ymx2+ymn2)*(hint1(1)**2/4.0d0-hipr1(8)**2)
     &          *g(y1,y2,pt2)/2.0d0

        end

c***********************************************************************

      function fjetrig(x,wgt)

c...Provide function for the calculation of triggered jet cross section.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/
      dimension x(10)

      ptmin=abs(hipr1(10))-0.25d0
      ptmin=max(ptmin,hipr1(8))
      am2=0.d0
      if(ihpr2(3).eq.3) then
         am2=hipr1(7)**2
         ptmin=max(0.0d0,hipr1(10))
      endif
      ptmax=abs(hipr1(10))+0.25d0
      if(hipr1(10).le.0.0d0) ptmax=hint1(1)/2.0d0-am2
      if(ptmax.le.ptmin) ptmax=ptmin+0.25d0
      pt2=(ptmax**2-ptmin**2)*x(1)+ptmin**2
      amt2=pt2+am2
      xt=2.0d0*sqrt(amt2)/hint1(1)
      ymx1=dlog(1.0d0/xt+sqrt(1.0d0/xt**2-1.0d0))
      y1=2.0d0*ymx1*x(2)-ymx1
      ymx2=dlog(2.0d0/xt-dexp(y1))
      ymn2=dlog(2.0d0/xt-dexp(-y1))
      y2=(ymx2+ymn2)*x(3)-ymn2

      if(ihpr2(3).eq.3) then
         gtrig=2.0d0*ghvq(y1,y2,amt2)
c...Only direct photon production.
      else if(ihpr2(3).eq.2) then
         gtrig=2.0d0*gphoton(y1,y2,pt2)
c...QCD processes.
      else
         gtrig=g(y1,y2,pt2)
      endif

      fjetrig=2.0d0*ymx1*(ymx2+ymn2)*(ptmax**2-ptmin**2)
     &          *gtrig/2.0d0

      end

c***********************************************************************

      function ghvq(y1,y2,amt2)

      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/
      dimension f(2,7)

      xt=2.0d0*sqrt(amt2)/hint1(1)
      x1=0.50d0*xt*(dexp(y1)+dexp(y2))
      x2=0.50d0*xt*(dexp(-y1)+dexp(-y2))
      ss=x1*x2*hint1(1)**2
      af=4.0d0
      if(ihpr2(18).ne.0) af=5.0d0
      dlam=hipr1(15)
      aph=12.0d0*3.1415926d0/(33.0d0-2.0d0*af)/dlog(amt2/dlam**2)
 
      call jmparton(f,x1,x2,amt2)
 
      gqq=4.0d0*( cosh(y1-y2)+hipr1(7)**2/amt2 
     &  )/(1.d0+cosh(y1-y2))/9.0d0
     &    *( f(1,1)*f(2,2)+f(1,2)*f(2,1)+f(1,3)*f(2,4)
     &      +f(1,4)*f(2,3)+f(1,5)*f(2,6)+f(1,6)*f(2,5) )
      ggg=(8.d0*cosh(y1-y2)-1.d0)*(cosh(y1-y2)+2.0d0*hipr1(7)**2/amt2
     &    -2.0d0*hipr1(7)**4/amt2**2)/(1.0d0+cosh(y1-y2))/24.d0
     &    *f(1,7)*f(2,7)
 
      ghvq=(gqq+ggg)*hipr1(23)*3.14159d0*aph**2/ss**2

      end

c***********************************************************************

      function gphoton(y1,y2,pt2)

c...Provide function for the calculation of direct photon.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/
      dimension f(2,7)

      xt=2.0d0*sqrt(pt2)/hint1(1)
      x1=0.50d0*xt*(dexp(y1)+dexp(y2))
      x2=0.50d0*xt*(dexp(-y1)+dexp(-y2))
      z=sqrt(min(0.0d0,1.d0-xt**2/x1/x2))
      ss=x1*x2*hint1(1)**2
      t=-(1.0d0-z)/2.0d0
      u=-(1.0d0+z)/2.0d0
      af=3.0d0
      dlam=hipr1(15)
      aph=12.0d0*3.1415926d0/(33.0d0-2.0d0*af)/dlog(pt2/dlam**2)
      aphem=1.0d0/137.0d0
 
      call jmparton(f,x1,x2,pt2)
 
      g11=-(u**2+1.0d0)/u/3.0d0*f(1,7)*(
     $  4.0d0*f(2,1)+4.0d0*f(2,2)+f(2,3)+f(2,4)+f(2,5)+f(2,6)
     $       )/9.0d0
      g12=-(t**2+1.0d0)/t/3.0d0 * f(2,7)*( 4.0d0*f(1,1)+4.0d0*f(1,2)
     &           +f(1,3)+f(1,4)+f(1,5)+f(1,6) )/9.0d0
      g2=8.0d0*(u**2+t**2)/u/t/9.0d0*( 4.0d0*f(1,1)*f(2,2)
     &   +4.0d0*f(1,2)*f(2,1)+f(1,3)*f(2,4)+f(1,4)*f(2,3)
     &   +f(1,5)*f(2,6)+f(1,6)*f(2,5) )/9.0d0
 
      gphoton=(g11+g12+g2)*hipr1(23)*3.14159d0*aph*aphem/ss**2

      end

c***********************************************************************

      function g(y1,y2,pt2)

c...Provide function for the calculation of QCD processes.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/
      dimension f(2,7)

      xt=2.0d0*sqrt(pt2)/hint1(1)
      x1=0.50d0*xt*(exp(y1)+exp(y2))
      x2=0.50d0*xt*(exp(-y1)+exp(-y2))
      z=sqrt(max(0.0d0,1.d0-xt**2/x1/x2))
      ss=x1*x2*hint1(1)**2
      t=-(1.0d0-z)/2.0d0
      u=-(1.0d0+z)/2.0d0

      af=3.0d0
      dlam=hipr1(15)       ! Scale Lambda of alpha_s
      aph=12.0d0*3.1415926d0/(33.0d0-2.0d0*af)/dlog(pt2/dlam**2)
 
c....Structure function
      call jmparton(f,x1,x2,pt2)
 
c...q q'-> qq' t-channel gluon exchange

      g11=( ( f(1,1) + f(1,2) )*( f(2,3) + f(2,4) + f(2,5) + f(2,6) )
     &  + ( f(1,3) + f(1,4) )*( f(2,5) + f(2,6) ) )*subcrs1(t,u)
 
      g12=( ( f(2,1) + f(2,2) )*( f(1,3) + f(1,4) + f(1,5) + f(1,6) )
     &   +( f(2,3) + f(2,4) )*( f(1,5) + f(1,6) ) )*subcrs1(u,t)
 
      g13=( f(1,1)*f(2,1) + f(1,2)*f(2,2) + f(1,3)*f(2,3)
     $    + f(1,4)*f(2,4) + f(1,5)*f(2,5) + f(1,6)*f(2,6)
     $     )*( subcrs1(u,t)+subcrs1(t,u)-8.d0/t/u/27.d0 )
 
c...q q'-> qq' s-channel gluon exchange.
      g2=(af-1)*( f(1,1)*f(2,2) + f(2,1)*f(1,2) + f(1,3)*f(2,4)
     &          + f(2,3)*f(1,4) + f(1,5)*f(2,6) + f(2,5)*f(1,6) )
     $   *subcrs2(t,u)
 
c...qq~->qq~  same flavour.
      g31=( f(1,1)*f(2,2) + f(1,3)*f(2,4) + f(1,5)*f(2,6) )*subcrs3(t,u)
      g32=( f(2,1)*f(1,2) + f(2,3)*f(1,4) + f(2,5)*f(1,6) )*subcrs3(u,t)
 
c...qq~->gg.
      g4=( f(1,1)*f(2,2) + f(2,1)*f(1,2)
     $   + f(1,3)*f(2,4) + f(2,3)*f(1,4)
     1   + f(1,5)*f(2,6) + f(2,5)*f(1,6) )*subcrs4(t,u)
 
c...gg->qq~.
      g5=af*f(1,7)*f(2,7)*subcrs5(t,u)
 
c...gq->gq.
      g61=f(1,7)*(
     $   f(2,1) + f(2,2) + f(2,3) + f(2,4) + f(2,5) +f(2,6)
     $            )*subcrs6(t,u)
      g62=f(2,7)*(
     $   f(1,1) + f(1,2) + f(1,3) + f(1,4) + f(1,5) +f(1,6)
     $           )*subcrs6(u,t)
 
c...gg->gg.
      g7=f(1,7)*f(2,7)*subcrs7(t,u)
 
      g=(g11+g12+g13+g2+g31+g32+g4+g5+g61+g62+g7)*hipr1(17)*
     1   3.14159d0*aph**2/ss**2

      end

c***********************************************************************

      function subcrs1(t,u)
c...q q'-> qq' t-channel gluon exchange.
      implicit double precision(a-h, o-z)
      subcrs1=4.d0/9.d0*(1.d0+u**2)/t**2
      end

c***********************************************************************

      function subcrs2(t,u)
c...q q'-> qq' s-channel gluon exchange.
      implicit double precision(a-h, o-z)
      subcrs2=4.d0/9.d0*(t**2+u**2)
      end

c***********************************************************************

      function subcrs3(t,u)
c...qq~->qq~  same flavour.
      implicit double precision(a-h, o-z)
      subcrs3=4.d0/9.d0*(t**2+u**2+(1.d0+u**2)/t**2
     1                              -2.d0*u**2/3.d0/t)
      end

c***********************************************************************

      function subcrs4(t,u)
c...qq~->gg.
      implicit double precision(a-h, o-z)
      subcrs4=8.d0/3.d0*(t**2+u**2)*(4.d0/9.d0/t/u-1.d0)
      end

c***********************************************************************

      function subcrs5(t,u)
c...gg->qq~.
      implicit double precision(a-h, o-z)
      subcrs5=3.d0/8.d0*(t**2+u**2)*(4.d0/9.d0/t/u-1.d0)
      end

c***********************************************************************

      function subcrs6(t,u)
c...gq->gq.
      implicit double precision(a-h, o-z)
      subcrs6=(1.d0+u**2)*(1.d0/t**2-4.d0/u/9.d0)
      end

c***********************************************************************

      function subcrs7(t,u)
c...gg->gg.
      implicit double precision(a-h, o-z)
      subcrs7=9.d0/2.d0*(3.d0-t*u-u/t**2-t/u**2)
      end

c***********************************************************************

      subroutine jmparton(f,x1,x2,qq)

c...Calculate the structure function with nuclear effect.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/njet/n,ip_crs
      save  /hiparnt/,/njet/
      dimension f(2,7) 
      parameter(ihij=0)
      double precision q2,xr,xpq(-25:25)

      if(ihij.eq.0) then
        q2=dble(qq)
        xr=dble(x1)

c         if(xr.le.0.0.or.x1.ge.1.) then
c           write(6,800) xr,x1
c800  format(' Error: x1 value outside physical range; x =',
c    $       1p,e12.3,1x,d22.15)
c           return
c         endif

          call pjpdfu(ihnt2(5),xr,q2,xpq)
          f(1,1)=xpq(2)           ! u
          f(1,2)=xpq(-2)          ! u~
          f(1,3)=xpq(1)           ! d
          f(1,4)=xpq(-1)          ! d~
          f(1,5)=xpq(3)           ! s
          f(1,6)=xpq(-3)          ! s~
          f(1,7)=xpq(0)           ! g
          q2=dble(qq)
          xr=dble(x2)

c         if(xr.le.0.0.or.x1.ge.1.) then
c           write(6,801) xr,x2
c801  format(' Error: x2 value outside physical range; x =',
c    $       1p,e12.3,1x,d22.15)
c           return
c         endif

          call pjpdfu(ihnt2(6),xr,q2,xpq)
          f(2,1)=xpq(2)           ! u
          f(2,2)=xpq(-2)          ! u~
          f(2,3)=xpq(1)           ! d
          f(2,4)=xpq(-1)          ! d~
          f(2,5)=xpq(3)           ! s
          f(2,6)=xpq(-3)          ! s~
          f(2,7)=xpq(0)           ! g
c         write(70,*)'kf1 kf2 x qq ',ihnt2(5),ihnt2(6),x1,qq
c         write(70,'(7(f9.6,1x))')(f(1,i),i=1,7)
          goto  5000
      endif

      dlam=hipr1(15) ! the scale Lambda
      q0=hipr1(16)   ! initial scale Q_0
      s=dlog(dlog(qq/dlam**2)/dlog(q0**2/dlam**2))

      if(ihpr2(7).eq.2) go to 200

c...P.R.D30 (1984)49
c....Duke Owens Set 1
C*******************************************************
c...xd_v
      at1=0.419d0+0.004d0*s-0.007d0*s**2
      at2=3.460d0+0.724d0*s-0.066d0*s**2
      gmud=4.40d0-4.86d0*s+1.33d0*s**2
      at3=0.763d0-0.237d0*s+0.026d0*s**2
      at4=4.00d0+0.627d0*s-0.019d0*s**2
      gmd=-0.421d0*s+0.033d0*s**2
C*******************************************************
c...xS
      cas=1.265d0-1.132d0*s+0.293d0*s**2
      as=-0.372d0*s-0.029d0*s**2
      bs=8.05d0+1.59d0*s-0.153d0*s**2
      aphs=6.31d0*s-0.273d0*s**2
      btas=-10.5d0*s-3.17d0*s**2
      gms=14.7d0*s+9.80d0*s**2
C********************************************************
c...xc
C     CAC=0.135*S-0.075*S**2
C     AC=-0.036-0.222*S-0.058*S**2
C     BC=6.35+3.26*S-0.909*S**2
C     APHC=-3.03*S+1.50*S**2
C     BTAC=17.4*S-11.3*S**2
C     GMC=-17.9*S+15.6*S**2
C***********************************************************
c...xG
      cag=1.56d0-1.71d0*s+0.638d0*s**2
      ag=-0.949d0*s+0.325d0*s**2
      bg=6.0d0+1.44d0*s-1.05d0*s**2
      aphg=9.0d0-7.19d0*s+0.255d0*s**2
      btag=-16.5d0*s+10.9d0*s**2
      gmg=15.3d0*s-10.1d0*s**2
      go to 300

C********************************************************
c....Duke Owens Set 2
200   at1=0.374d0+0.014d0*s
      at2=3.33d0+0.753d0*s-0.076d0*s**2
      gmud=6.03d0-6.22d0*s+1.56d0*s**2
      at3=0.761d0-0.232d0*s+0.023d0*s**2
      at4=3.83d0+0.627d0*s-0.019d0*s**2
      gmd=-0.418d0*s+0.036d0*s**2
C************************************
      cas=1.67d0-1.92d0*s+0.582d0*s**2
      as=-0.273d0*s-0.164d0*s**2
      bs=9.15d0+0.530d0*s-0.763d0*s**2
      aphs=15.7d0*s-2.83d0*s**2
      btas=-101.0d0*s+44.7d0*s**2
      gms=223.0d0*s-117.0d0*s**2
C*********************************
C     CAC=0.067*S-0.031*S**2
C     AC=-0.120-0.233*S-0.023*S**2
C     BC=3.51+3.66*S-0.453*S**2
C     APHC=-0.474*S+0.358*S**2
C     BTAC=9.50*S-5.43*S**2
C     GMC=-16.6*S+15.5*S**2
C**********************************
      cag=0.879d0-0.971d0*s+0.434d0*s**2
      ag=-1.16d0*s+0.476d0*s**2
      bg=4.0d0+1.23d0*s-0.254d0*s**2
      aphg=9.0d0-5.64d0*s-0.817d0*s**2
      btag=-7.54d0*s+5.50d0*s**2
      gmg=-0.596d0*s+1.26d0*s**2

C*********************************
300   b12=dexp(gmre(at1)+gmre(at2+1.d0)-gmre(at1+at2+1.d0))
      b34=dexp(gmre(at3)+gmre(at4+1.d0)-gmre(at3+at4+1.d0))
      cnud=3.d0/b12/(1.d0+gmud*at1/(at1+at2+1.d0))
      cnd=1.d0/b34/(1.d0+gmd*at3/(at3+at4+1.d0))
C********************************************************
C     FUD=X*(U+D)
C     FS=X*2(UBAR+DBAR+SBAR)  AND UBAR=DBAR=SBAR
C*******************************************************
      fud1=cnud*x1**at1*(1.d0-x1)**at2*(1.d0+gmud*x1)
      fs1=cas*x1**as*(1.d0-x1)**bs*(1.d0+aphs*x1
     &    +btas*x1**2+gms*x1**3)

c...d
      f(1,3)=cnd*x1**at3*(1.d0-x1)**at4*(1.d0+gmd*x1)+fs1/6.0d0
c...u
      f(1,1)=fud1-f(1,3)+fs1/3.0d0

      f(1,2)=fs1/6.d0
      f(1,4)=fs1/6.d0
      f(1,5)=fs1/6.d0
      f(1,6)=fs1/6.d0
c...gluon
      f(1,7)=cag*x1**ag*(1.d0-x1)**bg*(1.d0+aphg*x1
     &       +btag*x1**2+gmg*x1**3)
 
      fud2=cnud*x2**at1*(1.d0-x2)**at2*(1.d0+gmud*x2)
      fs2=cas*x2**as*(1.d0-x2)**bs*(1.d0+aphs*x2
     &    +btas*x2**2+gms*x2**3)
      f(2,3)=cnd*x2**at3*(1.d0-x2)**at4*(1.d0+gmd*x2)+fs2/6.0d0
      f(2,1)=fud2-f(2,3)+fs2/3.d0
      f(2,2)=fs2/6.d0
      f(2,4)=fs2/6.d0
      f(2,5)=fs2/6.d0
      f(2,6)=fs2/6.d0
      f(2,7)=cag*x2**ag*(1.d0-x2)**bg*(1.d0+aphg*x2
     &       +btag*x2**2+gmg*x2**3)

C***********Nuclear effect on the structure function****************
 5000 continue

      if(ihpr2(6).eq.1) then
        if(ihnt2(1).gt.1) then
          aax=1.193d0*dlog(dble(ihnt2(1)))**0.16666666d0
          rrx=aax*(x1**3-1.2d0*x1**2+0.21d0*x1)+1.0d0
     &       +1.079d0*(dble(ihnt2(1))**0.33333333d0-1.0d0)
     &       /dlog(ihnt2(1)+1.0d0)*dsqrt(x1)*dexp(-x1**2/0.01d0)
          if(ip_crs.eq.1 .or.ip_crs.eq.3) rrx=dexp(-x1**2/0.01d0)
          do 400 i=1,7
           f(1,i)=rrx*f(1,i)
 400      continue
        endif
        if(ihnt2(3).gt.1) then
          aax=1.193d0*dlog(dble(ihnt2(3)))**0.16666666d0
          rrx=aax*(x2**3-1.2d0*x2**2+0.21d0*x2)+1.0d0
     &        +1.079d0*(dble(ihnt2(3))**0.33333d0-1.0d0)
     &        /dlog(ihnt2(3)+1.0d0)*dsqrt(x2)*dexp(-x2**2/0.01d0)
          if(ip_crs.eq.2 .or. ip_crs.eq.3) rrx=dexp(-x2**2/0.01d0)
          do 500 i=1,7
              f(2,i)=rrx*f(2,i)
 500      continue
        endif
      endif
 
      end

c***********************************************************************

      function gmre(x)

c...Log of Gamma function.
      implicit double precision(a-h, o-z)

      z=x
      if(x.gt.3.0d0) go to 10
      z=x+3.d0

10    gmre=0.5d0*dlog(2.d0*3.14159265d0/z)+z*dlog(z)-z+dlog(1.d0
     1  +1.d0/12.d0/z+1.d0/288.d0/z**2-139.d0/51840.d0/z**3
     1  -571.d0/2488320.d0/z**4)
      if(z.eq.x) go to 20
      gmre=gmre-dlog(z-1.d0)-dlog(z-2.d0)-dlog(z-3.d0)
20    continue

      end

c***********************************************************************
c
c                  Funtions for soft processes
c
c***********************************************************************

      function fnkick(x)

c...Pupose: to define function for pt**2 kick from soft interaction.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      fnkick=1.0d0/(x+parc(67)**2)/(x+parc(68)**2)
     &          /(1+exp((sqrt(x)-parc(68))/0.4d0))

c     fnkick=x/(x*x+parc(67)**2)/(x*x+parc(68)**2)
c    &          /(1+exp((x-parc(68))/0.4d0))

      end

c***********************************************************************

      function fnkickr(x)

c...Pupose: to define function for pt**2 kick for resonance production.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      fnkickr=1.0d0/(x+parc(69)**2)/(x+parc(70)**2)
     $                       /(1+exp((sqrt(x)-parc(70))/0.3d0))

c     fnkickr=x/(x*x+parc(69)**2)/(x*x+parc(70)**2)
c    $                       /(1+exp((x-parc(70))/0.3d0))
      end

c***********************************************************************

      function fnkick2(x)

c...Define the expotnetial pt distribution of sea quarks.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/
      fnkick2=x*exp(-2.0d0*x/hipr1(42))
      end

c***********************************************************************

      function fnstru(x)

c...Pupose: to define function of x dist. of valence quarks for baryons.
      implicit double precision(a-h, o-z)
        common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
        save  /hiparnt/
        fnstru=(1.0d0-x)**hipr1(44)/
     &          (x**2+hipr1(45)**2/hint1(1)**2)**hipr1(46)
        end

c***********************************************************************

        function fnstrum(x)

c...Pupose: to define function of x dist. of valence quarks for mesons.
      implicit double precision(a-h, o-z)
        common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
        save  /hiparnt/
        fnstrum=1.0d0/((1.0d0-x)**2+hipr1(45)**2/hint1(1)**2)**hipr1(46)
     &          /(x**2+hipr1(45)**2/hint1(1)**2)**hipr1(46)
        end

c***********************************************************************

        function fnstrus(x)

c...Purpose: to define function of x dist. of valence quark
c...for the disassociated excitation in a single diffractive collision.

      implicit double precision(a-h, o-z)
        common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
        save  /hiparnt/
        fnstrus=(1.0d0-x)**hipr1(47)/
     &          (x**2+hipr1(45)**2/hint1(1)**2)**hipr1(48)
        end

c***********************************************************************

      subroutine hijwds(ia,idh,xhigh)

c...Sets up histogram IDH with radii for
c...nucleus IA distributed according to three param wood saxon.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/
      common/wood/r,d,fnorm,w
      save  /wood/
      dimension iaa(20),rr(20),dd(20),ww(20),rms(20)
      external rwdsax,wdsax

c...Parameters of special nuclei from atomic data and nuc data tables
c...VOL 14, 5-6 1974
        data iaa/4,12,16,27,32,40,56,63,93,184,197,208,8*0.d0/
        data rr/.964d0,2.355d0,2.608d0,2.84d0,3.458d0,3.766d0,3.971d0, 
     & 4.214d0,
     1        4.87d0,6.51d0,6.38d0,6.624d0,8*0.d0/
        data dd/.322d0,.522d0,.513d0,.569d0,.61d0,.586d0,.5935d0,.586d0,
     & .573d0,
     1        .535d0,.535d0,.549d0,8*0.d0/
        data ww/.517d0,-0.149d0,-0.051d0,0.d0,-0.208d0,-0.161d0,14*0.d0/
        data rms/1.71d0,2.46d0,2.73d0,3.05d0,3.247d0,3.482d0,3.737d0, 
     & 3.925d0,4.31d0,
     1        5.42d0,5.33d0,5.521d0,8*0.d0/
 
c...Set Wood-Sax params first, as in date et al.
        a=ia
c...D is Wood Sax diffuse param in fm.
        d=0.54d0
c...R is radius param
        r=1.19d0*a**(1.d0/3.d0) - 1.61d0*a**(-1.d0/3.d0)
c...W is the third of three Wood-Sax param.
        w=0.d0
 
c...Check table for special cases.
      do 10 i=1,12
        if (ia.eq.iaa(i)) then
          r=rr(i)
          d=dd(i)
          w=ww(i)
          rs=rms(i)
        end if
10    continue

c...fnorm is the normalize factor.
      fnorm=1.0d0
      xlow=0.d0
c     xhigh=r+ 12.*d
      xhigh=r+ 3.d0*d

      if (w.lt.-0.01d0)  then
        if (xhigh.gt.r/sqrt(abs(w))) xhigh=r/sqrt(abs(w))
      end if
      fgaus=gauss1(rwdsax,xlow,xhigh,0.001d0)
      fnorm=1.d0/fgaus

      if (idh.eq.1) then
         hint1(72)=r
         hint1(73)=d
         hint1(74)=w
         hint1(75)=fnorm/4.0d0/hipr1(40)
      else if (idh.eq.2) then
         hint1(76)=r
         hint1(77)=d
         hint1(78)=w
         hint1(79)=fnorm/4.0d0/hipr1(40)
      endif

c...Now set up hbook functions IDH for  R**2*RHO(R)
c...these histograms are used to generate random radii.
      call hifun(idh,xlow,xhigh,rwdsax)

      end

c***********************************************************************

        function wdsax(x)

c...Three parameter Wood Saxon.
      implicit double precision(a-h, o-z)
      common/wood/r,d,fnorm,w
      save  /wood/

      wdsax=fnorm*(1.d0+w*(x/r)**2)/(1+exp((x-r)/d))
      if (w.lt.0.d0) then
        if (x.ge.r/sqrt(abs(w))) wdsax=0.d0
      endif

      end

c***********************************************************************

      function rwdsax(x)
      implicit double precision(a-h, o-z)
      rwdsax=x*x*wdsax(x)
      end
 
c***********************************************************************

      function wdsax1(x)
c...Three parameter Wood Saxon for  projectile.
c...hint1(72)=r, hint1(73)=d, hint1(74)=w, hint1(75)=fnorm.

      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/

      wdsax1=hint1(75)*(1.d0+hint1(74)*(x/hint1(72))**2)/
     &       (1+exp((x-hint1(72))/hint1(73)))
      if (hint1(74).lt.0.d0) then
        if (x.ge.hint1(72)/sqrt(abs(hint1(74)))) wdsax1=0.d0
      endif
      end

c***********************************************************************

      function wdsax2(x)

c...Three parameter Wood Saxon for  target.
c...hint1(76)=r,hint1(77)=d, hint1(78)=w, hint1(79)=fnorm.

      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/
      wdsax2=hint1(79)*(1.d0+hint1(78)*(x/hint1(76))**2)/
     &     (1+exp((x-hint1(76))/hint1(77)))
      if (hint1(78).lt.0.d0) then
        if (x.ge.hint1(76)/sqrt(abs(hint1(78)))) wdsax2=0.d0
      endif

      end

c***********************************************************************
c...The next three subroutines are for Monte Carlo generation 
c...according to a given function FHB. One calls first HIFUN 
c...with assigned channel number I, low and up limits. Then to 
c...generate the distribution one can call HIRND(I) which gives 
c...you a random number generated according to the given function.
c 
c***********************************************************************

      subroutine hifun(ip,xmin,xmax,fhb)

c...Setup Monte Calro generation of a specific function.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/hijhb/rr(10,201),xx(10,201)
      dimension itag(10)
      character code*20
      save  /hijhb/
      save itag,mstc21
      external fhb
      data itag/10*0/,mstc21/0/

      if(mstc(21).ne.mstc21) then
         do j=1,10
           itag(j)=0
         end do
         mstc21=mstc21+1
      endif

      if(itag(ip).ne.0) then
        write(code,'(i20)')itag(ip)
        call jamerrm(30,0,'(hifun:) This number already used'//code)
      endif
      itag(ip)=1

      fnorm=gauss1(fhb,xmin,xmax,0.001d0)
      do 100 j=1,201
        xx(ip,j)=xmin+(xmax-xmin)*(j-1)/200.0d0
        xdd=xx(ip,j)
        rr(ip,j)=gauss1(fhb,xmin,xdd,0.001d0)/fnorm
100   continue

      end

c***********************************************************************

      function hirnd(i)

c...This generate random number according to a function.
      implicit double precision(a-h, o-z)
      common/hijhb/rr(10,201),xx(10,201)
      save  /hijhb/

      rx=rn(0)
      jl=0
      ju=202
10    if(ju-jl.gt.1) then
        jm=(ju+jl)/2
        if((rr(i,201).gt.rr(i,1)).eqv.(rx.gt.rr(i,jm))) then
          jl=jm
        else
          ju=jm
        endif
      go to 10
      endif
      j=jl
      if(j.lt.1) j=1
      if(j.ge.201) j=200
      hirnd=(xx(i,j)+xx(i,j+1))/2.0d0

      end     

c***********************************************************************

      function hirnd2(i,xmin0,xmax0)

c...This generate random number between xmin and xmax.

      implicit double precision(a-h, o-z)
      common/hijhb/rr(10,201),xx(10,201)
      save  /hijhb/

      xmin=xmin0
      xmax=xmax0
      if(xmin.lt.xx(i,1)) xmin=xx(i,1)
      if(xmax.gt.xx(i,201)) xmax=xx(i,201)
      jmin=1+200*(xmin-xx(i,1))/(xx(i,201)-xx(i,1))
      jmax=1+200*(xmax-xx(i,1))/(xx(i,201)-xx(i,1))
      rx=rr(i,jmin)+(rr(i,jmax)-rr(i,jmin))*rn(0)
      jl=0
      ju=202
10    if(ju-jl.gt.1) then
         jm=(ju+jl)/2
         if((rr(i,201).gt.rr(i,1)).eqv.(rx.gt.rr(i,jm))) then
            jl=jm
         else
            ju=jm
         endif
      go to 10
      endif

      j=jl
      if(j.lt.1) j=1
      if(j.ge.201) j=200
c     rx=rr(i,jmin)+(rr(i,jmax)-rr(i,jmin))*rn(0)
c     hirnd2=(xx(i,j)+xx(i,j+1))/2.0d0
      rx=(rx-rr(i,j))/(rr(i,j+1)-rr(i,j))
      hirnd2=xx(i,j)+(xx(i,j+1)-xx(i,j))*rx

      end     


c***********************************************************************

      subroutine jamfun(ip,xmin0,xmax0,fhb)

c...Setup Monte Calro generation of a specific function.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter(nbin=200)
      common/jamhb/rr(10,0:nbin),alp(10,0:nbin),bet(10,0:nbin)
     $ ,gam(10,0:nbin),xmin,xmax
      save  /jamhb/
      dimension itag(10)
c     character code*20
c     save itag,mstc21
      external fhb
c     data itag/10*0/,mstc21/0/
c     if(mstc(21).ne.mstc21) then
c        do j=1,10
c          itag(j)=0
c        end do
c        mstc21=mstc21+1
c     endif
c     if(itag(ip).ne.0) then
c       write(6,*)'(jamun:) This number already used',itag(ip)
c     endif
c     itag(ip)=1

      xmin=xmin0
      xmax=xmax0
      dx = (xmax-xmin)/nbin

      rr(ip,0)=0d0
      do i=0,nbin-1
        x1 = xmin + i*dx
        ggg=gauss1(fhb,x1,x1+dx,1d-5)
c       ggg=pjgaus(fhb, x1, x1+dx, 1d-5)
        rr(ip,i+1)=rr(ip,i)+ggg
      end do

      total = rr(ip,nbin)
      if(total.le.0.0d0) then
        write(6,*)'jamfun: toal integral <0',total
        stop
      endif

      do i=1,nbin
        rr(ip,i) = rr(ip,i)/total
      end do

      do i=0,nbin-1
        x0 = xmin+i*dx
        r2=rr(ip,i+1)-rr(ip,i)
        r1=gauss1(fhb,x0,x0+0.5*dx,1d-5)/total
c       r1=pjgaus(fhb,x0,x0+0.5*dx,1d-5)/total
        gam(ip,i)=(2*r2-4*r1)/(dx*dx)
        bet(ip,i)=r2/dx-gam(ip,i)*dx
        alp(ip,i)=x0
        gam(ip,i) = gam(ip,i)*2
      end do

      end

c***********************************************************************

      double precision function jamrnd2(ip,xmin0,xmax0)

c...This generate random number between xmin and xmax.

      implicit none
      real*8 xmin0,xmax0
      integer nbin,ip
      parameter(nbin=200)
      real*8 rr,alp,bet,gam,xmin,xmax 
      common/jamhb/rr(10,0:nbin),alp(10,0:nbin),bet(10,0:nbin)
     $ ,gam(10,0:nbin),xmin,xmax
      save  /jamhb/
      integer nbinmin,nbinmax,ii,middle,nabove,nbelow,itry
      real*8 pmin,pmax,r,rr1,xx,dx,x,rn

      
      if(abs(xmax0-xmin0).lt.1e-5) then
        jamrnd2=xmax0
        return
      endif

      dx = (xmax-xmin)/nbin
      nbinmin=int((xmin0-xmin)/dx)
      nbinmax=int((xmax0-xmin)/dx)+2
      if(nbinmax>nbin) nbinmax=nbin
      pmin=rr(ip,nbinmin)
      pmax=rr(ip,nbinmax)
      
      itry=0
200   r=pmin+(pmax-pmin)*rn(0)
      itry=itry+1
      if(itry.ge.100) then
        print *,'jamrnd2 infinit loop?',itry,' ip=',ip
        print *,'pmin pmax',pmin,pmax,xmin0,xmax0
        stop
      endif
      nabove=nbin+1
      nbelow=0
100   middle=(nabove+nbelow)/2
      if(r .eq. rr(ip,middle-1)) then
        ii=middle-1
        goto 110
      endif
      if(r .lt. rr(ip,middle-1)) then
         nabove=middle
      else
         nbelow=middle
      endif
      if(nabove-nbelow .gt.1) goto 100
      ii=nbelow-1
110   continue      
      rr1 = r - rr(ip,ii)
      if(gam(ip,ii).gt.0d0) then
        xx = (-bet(ip,ii)+sqrt(bet(ip,ii)**2+2*gam(ip,ii)*rr1))
     $       /gam(ip,ii)
      else
        xx=rr1/bet(ip,ii)
      endif
      x=alp(ip,ii)+xx
      if(x.lt.xmin0 .or. x.gt.xmax0) goto 200

      jamrnd2=x

      end     





c***********************************************************************
c
c                  General utilities
c
c***********************************************************************

      subroutine vegas0(fxn,avgi,sd,chi2a)

c...Purpose: to perform N-dimensional monte carlo integ'n
c...- by G.P. Lepage   SEPT 1976/(REV)APR 1978.

      implicit double precision(a-h, o-z)
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn,iprn,ierr
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      common/bveg3/f,ti,tsi   
      save  /bveg1/,/bveg2/,/bveg3/
      external fxn
      dimension d(50,10),di(50,10),xin(50),r(50),dx(10),dt(10),x(10)
     1   ,kg(10),ia(10)
      double precision qran(10)
      data ndmx/50/,alph/1.5d0/,one/1.d0/,mds/-1/

      ndo=1
      do 1 j=1,ndim
1     xi(1,j)=one

c***********************************************************************

      entry vegas1(fxn,avgi,sd,chi2a)

c...Initializes cummulative variables, but not grid.

      ierr=0
      it=0
      si=0.d0
      si2=si
      swgt=si
      schi=si

c***********************************************************************

      entry vegas2(fxn,avgi,sd,chi2a)

c...No initialization.

      nd=ndmx
      ng=1
      if(mds.eq.0) go to 2

      ng=(ncall/2.d0)**(1.d0/ndim)
      mds=1

      if((2*ng-ndmx).lt.0) go to 2
      mds=-1
      npg=ng/ndmx+1
      nd=ng/npg
      ng=npg*nd

2     continue

      k=ng**ndim
      npg=ncall/k
      if(npg.lt.2) npg=2
      calls=npg*k
      dxg=one/ng
      dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-one)
      xnd=nd
      ndm=nd-1
      dxg=dxg*xnd
      xjac=one/calls

      do 3 j=1,ndim
c***this is the line 50
        dx(j)=xu(j)-xl(j)
        xjac=xjac*dx(j)
3     continue

c...Rebin preserving bin density.

      if(nd.eq.ndo) go to 8

      rc=ndo/xnd
      do 7 j=1,ndim

        k=0
        xn=0.d0
        dr=xn
        i=k

4       continue
        k=k+1
        dr=dr+one
        xo=xn
        xn=xi(k,j)
5       continue
        if(rc.gt.dr) go to 4

        i=i+1
        dr=dr-rc
        xin(i)=xn-(xn-xo)*dr

        if(i.lt.ndm) go to 5

        do 6 i=1,ndm
          xi(i,j)=xin(i)
6       continue
        xi(nd,j)=one
7     continue
      ndo=nd
8     continue
      if(nprn.ne.0) write(iprn,200) ndim,calls,it,itmx,acc,mds,nd
     1                           ,(xl(j),xu(j),j=1,ndim)

c***********************************************************************

      entry vegas3(fxn,avgi,sd,chi2a)

c...Main integration loop.

9     continue
      it=it+1
      ti=0.d0
      tsi=ti
      do 10 j=1,ndim
        kg(j)=1
        do 10 i=1,nd
          d(i,j)=ti
          di(i,j)=ti
10    continue

11    continue

      fb=0.d0
      f2b=fb
      k=0
12    continue

      k=k+1
      call aran9(qran,ndim)
      wgt=xjac

      do 15 j=1,ndim

        xn=(kg(j)-qran(j))*dxg+one
c*****this is the line 100
        ia(j)=xn

        if(ia(j).gt.1.and.ia(j).le.50) then
          xo=xi(ia(j),j)-xi(ia(j)-1,j)
          rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
        else if(ia(j).eq.1) then
          xo=xi(ia(j),j)
          rc=(xn-ia(j))*xo
        else
          print *,'(vega:) error ia(j)>50',j,ia(j),kg(j),qran(j),dxg
          ierr=ia(j)
          ia(j)=min(50,ia(j))
          xo=xi(ia(j),j)-xi(ia(j)-1,j)
          rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
        endif

        x(j)=xl(j)+rc*dx(j)
        wgt=wgt*xo*xnd

15    continue

      f=wgt
      f=f*fxn(x,wgt)
      f2=f*f
      fb=fb+f
      f2b=f2b+f2
      do 16 j=1,ndim
        di(ia(j),j)=di(ia(j),j)+f
        if(mds.ge.0) d(ia(j),j)=d(ia(j),j)+f2
16    continue
      if(k.lt.npg) go to 12

      f2b=dsqrt(f2b*npg)
      f2b=(f2b-fb)*(f2b+fb)
      ti=ti+fb
      tsi=tsi+f2b
      if(mds.ge.0) go to 18

      do 17 j=1,ndim
        d(ia(j),j)=d(ia(j),j)+f2b
17    continue

18    continue
      k=ndim

19    continue
      kg(k)=mod(kg(k),ng)+1
      if(kg(k).ne.1) go to 11
      k=k-1
      if(k.gt.0) go to 19

c...Final results for this iteration.

      tsi=tsi*dv2g
      ti2=ti*ti
      wgt=ti2/(tsi+1.0d-37)
      si=si+ti*wgt
      si2=si2+ti2
      swgt=swgt+wgt
      swgt=swgt+1.0d-37
      si2=si2+1.0d-37
      schi=schi+ti2*wgt
      avgi=si/(swgt)
      sd=swgt*it/(si2)
      chi2a=sd*(schi/swgt-avgi*avgi)/(it-.999d0)
      sd=dsqrt(one/sd)

c****this is the line 150.
      if(nprn.eq.0) go to 21

      tsi=dsqrt(tsi)
      write(iprn,201) it,ti,tsi,avgi,sd,chi2a
      if(nprn.ge.0) go to 21

      do 20 j=1,ndim
20    write(iprn,202) j,(xi(i,j),di(i,j),d(i,j),i=1,nd)

c...Refine grid.

21    do 23 j=1,ndim
        xo=d(1,j)
        xn=d(2,j)
        d(1,j)=(xo+xn)/2.d0
        dt(j)=d(1,j)

        do 22 i=2,ndm
          d(i,j)=xo+xn
          xo=xn
          xn=d(i+1,j)
          d(i,j)=(d(i,j)+xn)/3.d0
          dt(j)=dt(j)+d(i,j)
22      continue
        d(nd,j)=(xn+xo)/2.d0
        dt(j)=dt(j)+d(nd,j)
23    continue

      do 28 j=1,ndim
      rc=0.d0
      do 24 i=1,nd
      r(i)=0.d0

      if (dt(j).ge.1.0d18) then
       write(6,*) '************** a singularity >1.0d18'
c      WRITE(5,1111)
c1111  FORMAT(1X,'**************IMPORTANT NOTICE***************')
c      WRITE(5,1112)
c1112  FORMAT(1X,'THE INTEGRAND GIVES RISE A SINGULARITY >1.0D18')
c      WRITE(5,1113)
c1113  FORMAT(1X,'PLEASE CHECK THE INTEGRAND AND THE LIMITS')
c      WRITE(5,1114)
c1114  FORMAT(1X,'**************END NOTICE*************')
      end if    

      if(d(i,j).le.1.0d-18) go to 24

      xo=dt(j)/d(i,j)
      r(i)=((xo-one)/xo/dlog(xo))**alph
      rc=rc+r(i)

24    continue

      rc=rc/xnd
      k=0
      xn=0.d0
      dr=xn
      i=k
25    k=k+1
      dr=dr+r(k)
      xo=xn

c****this is the line 200.
      xn=xi(k,j)
26    if(rc.gt.dr) go to 25
      i=i+1
      dr=dr-rc
      xin(i)=xn-(xn-xo)*dr/(r(k)+1.0d-30)
      if(i.lt.ndm) go to 26

      do 27 i=1,ndm
27    xi(i,j)=xin(i)

      xi(nd,j)=one

28    continue

      if(it.lt.itmx.and.acc*dabs(avgi).lt.sd) go to 9

200   format('0input parameters for vegas:  ndim=',i3,'  ncall=',f8.0
     1    /28x,'  it=',i5,'  itmx=',i5/28x,'  acc=',g9.3
     2    /28x,'  mds=',i3,'   nd=',i4/28x,'  (xl,xu)=',
     3    (t40,'( ',g12.6,' , ',g12.6,' )'))
201   format(///' integration by vegas' / '0iteration no.',i3,
     1    ':   integral =',g14.8/21x,'std dev  =',g10.4 /
     2    ' accumulated results:   integral =',g14.8 /
     3    24x,'std dev  =',g10.4 / 24x,'chi**2 per it''n =',g10.4)
202   format('0data for axis',i2 / ' ',6x,'x',7x,'  delt i  ',
     1    2x,' conv''ce  ',11x,'x',7x,'  delt i  ',2x,' conv''ce  '
     2   ,11x,'x',7x,'  delt i  ',2x,' conv''ce  ' /
     2    (' ',3g12.4,5x,3g12.4,5x,3g12.4))
      return
      end

c***********************************************************************

      subroutine aran9(qran,ndim)

c...Random number.
      implicit double precision(a-h, o-z)
      dimension qran(10)
c     common/seedvax/num1
      do 1 i=1,ndim
c   1 qran(i)=ran(num1)
    1 qran(i)=rn(0)
      end

c***********************************************************************

      function gauss1(f,a,b,eps)

c...Gaussian one-dimensional integration program.
      implicit double precision(a-h, o-z)
      external f
      dimension w(12),x(12)
      data const/1.0d-12/
      data w/0.1012285d0,.2223810d0,.3137067d0,.3623838d0,.0271525d0,
     &       .0622535d0,0.0951585d0,.1246290d0,.1495960d0,.1691565d0,
     &       .1826034d0,.1894506d0/
      data x/0.9602899d0,.7966665d0,.5255324d0,.1834346d0,.9894009d0,
     &       .9445750d0,0.8656312d0,.7554044d0,.6178762d0,.4580168d0,
     &       .2816036d0,.0950125d0/

        delta=const*abs(a-b)
        gauss1=0.0d0
        aa=a
5       y=b-aa
        if(abs(y).le.delta) return
2       bb=aa+y
        c1=0.5d0*(aa+bb)
        c2=c1-aa
        s8=0.0d0
        s16=0.0d0
        do 1 i=1,4
        u=x(i)*c2
1       s8=s8+w(i)*(f(c1+u)+f(c1-u))
        do 3 i=5,12
        u=x(i)*c2
3       s16=s16+w(i)*(f(c1+u)+f(c1-u))
        s8=s8*c2
        s16=s16*c2
        if(abs(s16-s8).gt.eps*(1.d0+abs(s16))) goto 4
        gauss1=gauss1+s16
        aa=bb
        goto 5
4       y=0.5d0*y
        if(abs(y).gt.delta) goto 2
        write(6,7)
        gauss1=0.0d0
        return
7       format(1x,'gauss1....too high acuracy required')

        end

c***********************************************************************

      function gauss2(f,a,b,eps)

      implicit double precision(a-h, o-z)
      external f
      dimension w(12),x(12)
      data const/1.0D-12/
      data w/0.1012285d0,.2223810d0,.3137067d0,.3623838d0,.0271525d0,
     &       .0622535d0,0.0951585d0,.1246290d0,.1495960d0,.1691565d0,
     &       .1826034d0,.1894506d0/
      data x/0.9602899d0,.7966665d0,.5255324d0,.1834346d0,.9894009d0,
     &       .9445750d0,0.8656312d0,.7554044d0,.6178762d0,.4580168d0,
     &       .2816036d0,.0950125d0/

        delta=const*abs(a-b)
        gauss2=0.0d0
        aa=a
5       y=b-aa
        if(abs(y).le.delta) return
2       bb=aa+y
        c1=0.5d0*(aa+bb)
        c2=c1-aa
        s8=0.0d0
        s16=0.0d0
        do 1 i=1,4
        u=x(i)*c2
1       s8=s8+w(i)*(f(c1+u)+f(c1-u))
        do 3 i=5,12
        u=x(i)*c2
3       s16=s16+w(i)*(f(c1+u)+f(c1-u))
        s8=s8*c2
        s16=s16*c2
        if(abs(s16-s8).gt.eps*(1.d0+abs(s16))) goto 4
        gauss2=gauss2+s16
        aa=bb
        goto 5
4       y=0.5d0*y
        if(abs(y).gt.delta) goto 2
        write(6,7)
        gauss2=0.0d0
        return
7       format(1x,'gauss2....too high acuracy required')
        end

