c Last modifications 07.04.2008 Compatibility gcc4 by T.Pierog
c 18.01.2007 update to compile with CORSIKA
c 05.09.2005 Add QGSJET-II subroutine definition to compile
c with preprocessing command -D__ALL__ without __QGSJETII__
c 18.05.04 Link file for conex with nexus
c author T. Pierog

#ifdef __STD__
#define __GHEISHA__
#define __QGSJET__
#define __ANALYSIS__
#endif

c-----------------------------------------------------------------------
      subroutine IniEpos
c-----------------------------------------------------------------------
#include "epos.inc"
#include "epos.incems"
      double precision seedp

c Conex common

      double precision airz,aira,airw,airavz,airava,airi
      common/cxair/airz(3),aira(3),airw(3),airavz,airava,airi(3)
      common/xsfiles/ixsfcx
      character*500 xsfnch,xsfnii,xsfnid,xsfnie,xsfnrj,xsfncs
      common/xsfname/  xsfnch, xsfnii, xsfnid, xsfnie, xsfnrj, xsfncs
      common/xsnfname/nxsfnch,nxsfnii,nxsfnid,nxsfnie,nxsfnrj,nxsfncs
      double precision xsegymin,xsegymax,xselab,xsecms,xsekin,xspnll
     *,xsengy
      common/xsenrgy/xsegymin,xsegymax,xselab,xsecms,xsekin,xspnll
     *,xsengy
      parameter(mxnodyxs=20)
      common/xsnodcy/nrnodyxs,nodyxs(mxnodyxs)
      common/xsappli/iapplxs,modelxs
      common/xsevent/neventxs,iframexs
c      data inicnt/0/
      data init/0/
      save init

      if(init.ge.1)return
      init=init+1

#ifdef __CXDEBUG__

      call utisx1('iniepos ',4)
      write(*,'(a)')'initialize EPOS ...'
#endif

      if(ilowegy.ne.1.or.MCleModel.eq.4)xsegymin=dble(0.5*egymin**2)
      if(MCModel.eq.4)xsegymax=min(xsegymax,dble(0.5*egymax**2))
      nrnody=nrnodyxs
      do i=1,nrnody
        nody(i)= nodyxs(i)
      enddo
#ifndef __CXCORSIKA__
      inicnt=inicnt+1
c      isetcs=2                  ! epos cross-section from tabulated calculation (h-A and AA)
      isetcs=3                  ! epos cross-section from tabulated simulations (h-A and A-A)
      isigma=0                  !do not print out the cross section on screen
      ionudi=1

      nfnii=nxsfnii             ! epos file name
      fnii=xsfnii
      nfnid=nxsfnid
      fnid=xsfnid
      nfnie=nxsfnie
      fnie=xsfnie
      nfnrj=nxsfnrj
      fnrj=xsfnrj
      nfncs=nxsfncs
      fncs=xsfncs
      nfnch=nxsfnch
      fnch=xsfnch

c air

      do i=1,3
        airanxs(i)=sngl(aira(i))
        airznxs(i)=sngl(airz(i))
        airwnxs(i)=sngl(airw(i))
      enddo
      airavanxs=sngl(airava)
      airavznxs=sngl(airavz)

      iappl  = iapplxs
      nevent = neventxs
      iframe = iframexs

      if(fnch(1:nfnch).ne.'none')
     &  open(ifcx,file=fnch(1:nfnch),status='unknown')

      call iclass(idproj,iclpro)
      call iclass(idtarg,icltar)
      if(inicnt.eq.1)then
        call ranfgt(seedp)      !not to change the seed ...
        call hdecin(.false.)
        call hnbspd(iospec)
        ktnbod=0
        call hnbpajini
        if(iclegy2.gt.1)then
          egyfac=(egymax*1.0001/egylow)**(1./float(iclegy2-1))
        else
          egyfac=1.
        endif
      endif
      maproj=64               !to set difnuc up to the maximum mass
      call conini
      call psaini
      call ranfst(seedp)        ! ... after this initialization

#endif

#ifdef __CXDEBUG__
      call utisx2
#endif

      end


c-----------------------------------------------------------------------
      subroutine IniEvtEpo
c-----------------------------------------------------------------------
c Initialization for each type of event (for given proj, targ and egy)
c-----------------------------------------------------------------------
#include "epos.inc"
      common/geom/rmproj,rmtarg,bmax,bkmx
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctain/mtain
      double precision rcproj,rctarg
      common/geom1/rcproj,rctarg
c Conex common
#ifdef __CXDEBUG__
      parameter (mxisx=200)
      character*500 subisx,textisx
      common/cxisx/isx,nisx,subisx(mxisx),isxsub(mxisx)   !also in gheisha_nexus
     &      ,isxsave,isxxsave,textisx
      character*500 fnho,fnck,fnwle,fnwhe,fndkz,fndkl,fndks,fndkm
     &,fnilo,fndke,fndkn,fndkg,fnwgh,fnwgl
      common /cxfiles/fnho,ifho,fnck,ifck,fnwle,ifwle,fnwhe,ifwhe
     &,fndkz,ifdkz,fndkl,ifdkl,fndks,ifdks,fndkm,ifdkm,fnilo,ifilo
     &,fndke,ifdke,fndkn,ifdkn,fndkg,ifdkg,fnwgh,ifwgh
     &,fnwgl,ifwgl
#endif
      double precision xsegymin,xsegymax,xselab,xsecms,xsekin,xspnll
     *,xsengy
      common/xsenrgy/xsegymin,xsegymax,xselab,xsecms,xsekin,xspnll
     *,xsengy
      common/xsnucl1/laprojxs,maprojxs,latargxs,matargxs
      common/xshadr2/idprojxs,idtargxs
      double precision xsamproj,xsamtarg,xsypjtl,xsyhaha,xspnullx
      common/xschadron/xsamproj,xsamtarg,xsypjtl,xsyhaha,xspnullx
      double precision xsrmproj,xsrmtarg,xsbmax,xsbminim,xsbmaxim,xsbkmx
      common/xsgeom/xsrmproj,xsrmtarg,xsbmax,xsbminim,xsbmaxim,xsbkmx
      double precision xsrcproj,xsrctarg
      common/xsgeom1/xsrcproj,xsrctarg
      double precision xstpro,xszpro,xsttar,xsztar,xsttaus
     &,xsdetap,xsdetat
      common/xscttaus/xstpro,xszpro,xsttar,xsztar,xsttaus
     &,xsdetap,xsdetat

      maproj=maprojxs
      laproj=laprojxs
      matarg=matargxs
      latarg=latargxs
      idproj=idprojxs
      idtarg=idtargxs
      amproj=xsamproj
      amtarg=xsamtarg

      engy=sngl(xsengy)
      elab=sngl(xselab)
      ecms=sngl(xsecms)
      ekin=sngl(xsekin)
      pnll=sngl(xspnll)
      pnullx=sngl(xspnullx)
      yhaha=sngl(xsyhaha)
      ypjtl=sngl(xsypjtl)
      detap=xsdetap
      detat=xsdetat
      tpro=xstpro
      zpro=xszpro
      ttar=xsttar
      ztar=xsztar

      xsbminim=dble(bminim)
      xsbmaxim=dble(bmaxim)

      call iclass(idproj,iclpro) **
      call iclass(idtarg,icltar) **
      call emsini(engy,idproj,idtarg)
      call paramini(1)
      bkmxndif=conbmxndif()
      bkmx=conbmx()
      xsbkmx=dble(bkmx)

      if(maproj.gt.1.or.matarg.gt.1)then
        xsbmax=xsrmproj+xsrmtarg
      else
        xsbmax=xsbkmx
      endif

      bimevt=-1
      bmax=sngl(xsbmax)
      rmproj=sngl(xsrmproj)
      rmtarg=sngl(xsrmtarg)
      rcproj=xsrcproj
      rctarg=xsrctarg
      call xsigma                          !set some variabkle according to xs
      if(idtarg.eq.0)idtarg=1120           !air = nucleus
#ifdef __CXDEBUG__
      if(isx.ge.2)write(ifck,*)
     &  'Epos used with (E,proj,maproj,matarg,bmax)',elab,idproj
     &  ,maproj,matarg,bmax
#endif
      return
      end

c-----------------------------------------------------------------------
      double precision function cxepocrse(ek,mapro,matar,id)
c-----------------------------------------------------------------------
c Epos cross section in double precision
c-----------------------------------------------------------------------
#include "epos.inc"
      double precision ek
      common/xshad10/iclproxs,icltarxs      !in conex.incnex

      iclpro=iclproxs   **
      icltar=icltarxs   **

      cxepocrse=dble(eposcrse(sngl(ek),mapro,matar,id))

      end

c-----------------------------------------------------------------------
      subroutine EposInput
c-----------------------------------------------------------------------
#include "epos.inc"
      common/cxsubro/isubin                 !also in conex.inc
      nopen=0
#ifdef __CXSUB__
      ifop=isubin
#else
      ifop=5
#endif
      call aread
      end

c-----------------------------------------------------------------------
      subroutine emsepo(iret)
c-----------------------------------------------------------------------
c  call epos to simulate interaction
c-----------------------------------------------------------------------
c Conex common
      implicit double precision (a-h,o-z)
#include "../conex.h"
#include "../conex.incnex"

c epos common
      parameter (mxptl=200000)   !max nr of particles in epos particle list
      common/col3/ncol,kolpt
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),ityptl(mxptl)
      real pptl,xorptl,tivptl,gaumx,wtmini,wtstep
      common/wana2/isphis,ispall,wtmini,wtstep,iwcent,iana,nbdky
      common/othe1/istore,istmax,gaumx,irescl,ntrymx,nclean,iopdg,ioidch


      double precision pfrx(mamxs),pfry(mamxs)
      integer          ityp(mamxs)

      nptlxs=0
      nevtxs=1
      iret=0
      irest = maprojxs*100+abs(laprojxs)
      inew=0

      call utpri('emsepo',ish,ishini,4)

      call emsaaa(iret)

      if(iret.ne.0)goto 1001

      ncolxs=ncol

      if(ish.ge.2)call alist('list before fragmentation&',1,nptl)

      call gakfra(iret)
      if(iret.gt.0)goto 1001
      if(ish.ge.2)
     &     call alist('list after fragmentation&',1,nptl)
      if(irescl.eq.1)then
        call utresc(iret)
        if(ish.ge.2)call alist('list after rescaling&',1,nptl)
        if(iret.gt.0)goto 1001
      endif

      nbdky=nptl
      call bjinta(ier)
      if(ier.eq.1)goto 1001
      if(ish.ge.2)
     &       call alist('list after int/decays&',1,nptl)

      nsxs=iabs(maprojxs)+iabs(matargxs)

      do is=1,nptl

        if(istptl(is).eq.0.or.is.le.nsxs)then
          nptlxs=nptlxs+1
          if(nptlxs.gt.mxptlxs)stop'Epos: mxptlxs too small'
          id=idptl(is)
          call cxidmass(id,am)

          istptlxs(nptlxs)=istptl(is)
          if ( is .le. maprojxs  .and.
     *          iorptl(is) .eq. 0  .and.  istptl(is) .eq. 0 ) then
            if ( ifragm .ne. 0 ) then
c  compose projectile spectators to remaining nucleus
              istptlxs(nptlxs)=1
              idrest = nptlxs
              if     ( id .eq. 1120 ) then
                inew  = inew + 101
                irest = irest - 101
              elseif ( id .eq. 1220 ) then
                inew  = inew + 100
                irest = irest - 100
              endif
            endif
          elseif( is .le. maprojxs  .and.
     *          iorptl(is) .eq. 0  .and.  istptl(is) .eq. 1 ) then
              jorptl(is)=1
          endif
          xsptl(1,nptlxs)=dble(pptl(1,is))
          xsptl(2,nptlxs)=dble(pptl(2,is))
          xsptl(3,nptlxs)=dble(pptl(3,is))
          xsptl(4,nptlxs)=dble(pptl(4,is))
          xsptl(5,nptlxs)=am    !mass
          ityptlxs(nptlxs)=ityptl(is)
          iorptlxs(nptlxs)=iorptl(is)
          jorptlxs(nptlxs)=jorptl(is)
          ifrptlxs(1,nptlxs)=ifrptl(1,is)
          ifrptlxs(2,nptlxs)=ifrptl(2,is)
          xsorptl(1,nptlxs)=dble(xorptl(1,is))
          xsorptl(2,nptlxs)=dble(xorptl(2,is))
          xsorptl(3,nptlxs)=dble(xorptl(3,is))
          xsorptl(4,nptlxs)=dble(xorptl(4,is))
          xstivptl(1,nptlxs)=dble(tivptl(1,is))
          xstivptl(2,nptlxs)=dble(tivptl(2,is))
          idptlxs(nptlxs)=id


#ifdef __CXDEBUG__
          if(isx.ge.5)write(ifck,'(a,i5,a,i5,a,4(e10.4,1x),f6.3)')
     $       ' particle from EPOS ',nptlxs,' id :',idptlxs(nptlxs)
     $  , ' momentum :',(xsptl(k,nptlxs),k=1,5)
#endif

        endif


      enddo
      if(inew.eq.0)goto 1000

      nptlxsa=nptlxs

      if( inew .eq. 100 .or. inew .eq. 101 ) then
c  remaining nucleus is single neutron or proton
        istptlxs(idrest)=0
        goto 1000

      elseif ( ifragm .ge. 2 ) then
c  remaining nucleus is evaporating nucleons and alpha particles
        jfin  = 0
        call cxvapor( maprojxs,inew,jfin,ityp,pfrx,pfry )
        if ( jfin .eq. 0 ) goto 1000
c loop to treat the remanents of the desintegrated fragment
        do  135  j = 1, jfin
          if(mod(ityp(j),100).eq.0)then
            inucl= abs(ityp(j))/100
            idnucl=100*inucl
          else
            inucl=1
            idnucl=ityp(j)
          endif
          ea = dble(inucl) * dble(pptl(4,idrest))
#ifdef __CXDEBUG__
          if (isx.ge.6) write(ifck,*) 'fragment: j,id,ea=',
     *                                       j,idnucl,ea
#endif
c  momenta squared
          call cxidmass(idnucl,am)
          ptm = ( ea - dble(am) ) * ( ea + dble(am) )
          pt2 = pfrx(j)**2 + pfry(j)**2
          if ( pt2 .ge. ptm ) then
#ifdef __CXDEBUG__
            if (isx.ge.2) write(ifck,*) 'emsepo: pt reject particle',j
#endif
            pfrx(j)=0.d0
            pfry(j)=0.d0
          endif
          plong=ptm-pt2
          if ( plong .ge. 0.d0)then
            plong  = sqrt( plong )
          else
#ifdef __CXDEBUG__
            if (isx.ge.2) write(ifck,*) 'emsepo: ptm reject particle',j
#endif
            do is=1,inucl
              if(istptlxs(is).eq.1.and.jorptlxs(is).eq.0)istptlxs(is)=0
            enddo
            goto 135
          endif
          nptlxs=nptlxs+1
          istptlxs(nptlxs)=0
          xsptl(1,nptlxs)=pfrx(j)
          xsptl(2,nptlxs)=pfry(j)
          xsptl(3,nptlxs)=plong
          xsptl(4,nptlxs)=ea
          xsptl(5,nptlxs)=am    !mass
          ityptlxs(nptlxs)=0
          iorptlxs(nptlxs)=1
          jorptlxs(nptlxs)=maprojxs
          ifrptlxs(1,nptlxs)=0
          ifrptlxs(2,nptlxs)=0
          xsorptl(1,nptlxs)=xsorptl(1,idrest)
          xsorptl(2,nptlxs)=xsorptl(2,idrest)
          xsorptl(3,nptlxs)=xsorptl(3,idrest)
          xsorptl(4,nptlxs)=xsorptl(4,idrest)
          xstivptl(1,nptlxs)=xstivptl(1,idrest)
          xstivptl(2,nptlxs)=xstivptl(2,idrest)
          idptlxs(nptlxs)=idnucl
 135    continue

        elseif ( ifragm .eq. 1 ) then
c  remaining nucleus is one fragment
          nptlxs=nptlxs+1
          istptlxs(nptlxs)=0
          xsptl(1,nptlxs)=0.d0
          xsptl(2,nptlxs)=0.d0
          xsptl(4,nptlxs)=0.d0
          inucl=0
          do is=1,maprojxs
            if(iorptl(is).eq.0.and.jorptl(is).eq.0)then
              inucl=inucl+1
              xsptl(4,nptlxs)=xsptl(4,nptlxs)+dble(pptl(4,is))
            endif
          enddo
          idnucl=100*inucl
          call cxidmass(idnucl,am)
          xsptl(5,nptlxs)=am    !mass
          ptot=(xsptl(4,nptlxs)+am)*(xsptl(4,nptlxs)-am)
          xsptl(3,nptlxs)=sqrt(ptot)
          ityptlxs(nptlxs)=0
          istptlxs(nptlxs)=0
          iorptlxs(nptlxs)=1
          jorptlxs(nptlxs)=maprojxs
          ifrptlxs(1,nptlxs)=0
          ifrptlxs(2,nptlxs)=0
          xsorptl(1,nptlxs)=xsorptl(1,idrest)
          xsorptl(2,nptlxs)=xsorptl(2,idrest)
          xsorptl(3,nptlxs)=xsorptl(3,idrest)
          xsorptl(4,nptlxs)=xsorptl(4,idrest)
          xstivptl(1,nptlxs)=xstivptl(1,idrest)
          xstivptl(2,nptlxs)=xstivptl(2,idrest)
          idptlxs(nptlxs)=idnucl
        endif

#ifdef __CXDEBUG__
      do is=nptlxsa+1,nptlxs
          if(isx.ge.5)write(ifck,'(a,i5,a,i5,a,4(e10.4,1x),f6.3)')
     $       ' particle from EPOS ',is,' id :',idptlxs(is)
     $  , ' momentum :',(xsptl(k,is),k=1,5)
      enddo
#endif


 1000 call utprix('emsepo',ish,ishini,4)
      return


1001  iret=1
      goto 1000

      end

#ifndef __CXCORSIKA__
c-----------------------------------------------------------------------
      subroutine ranfini(seed,iseq,iqq)
c-----------------------------------------------------------------------
      double precision    seed,dummy
      integer iseq,iqq,idum
      dummy=seed
      idum=iqq
      idum=iseq
      return
      end

c-----------------------------------------------------------------------
      subroutine ranfcv(seed)
c-----------------------------------------------------------------------
c Convert input seed to EPOS random number seed
c Useless with ranf
c-----------------------------------------------------------------------
      double precision seed,dummy
      dummy=seed

      return
      end
#endif
