*cmz :  3.14/16 01/11/90  22.18.07  by  rene brun
*-- author :
      subroutine sginit(ire,plab,n,ie,amt,amn,ecm,si,itar,inuc)
*      integer*2 ich,ibar,k1,k2
      common/abltis/am(110),ga(110),tau(110),ich(110),
     +              ibar(110),k1(110),k2(110)
      common /redver/ irii(17),ikii(17),ieii(17),
     +        thresh(268)
      common/reacn/siin(296),wkn(5184)
      common /reach/umo(296),plabf(296),siinh(296),wkh(5184),
     +              nrk(2,268),nure(30,2)
      ie=iefun(plab,ire)
      if (ie.le.ieii(ire)) ie=ie+1
      amt=am(itar)
      amn=am(n)
      amn2=amn*amn
      amt2=amt*amt
      ecm=sqrt(amn2+amt2+2.*amt*sqrt(amn2+plab**2))
c*** interpolation preparation
      ecmo=umo(ie)
      ecm1=umo(ie-1)
      decm=ecmo-ecm1
      dec=ecmo-ecm
      iiki=ikii(ire)+1
      eklim=-thresh(iiki)
      if(inuc.eq.0)then
         wok=siinh(ie)
         wdk=wok-siinh(ie-1)
      else
         wok=siin(ie)
         wdk=wok-siin(ie-1)
      endif
      if (ecm.gt.ecmo) wdk=0.
c*** interpolation in channel weights
      ielim=iefun(eklim,ire)
      delim=umo(ielim)+eklim +1.e-16
     
      dete=(ecm-(ecmo-eklim)*.5)*2.
      if (delim*delim-dete*dete) 20 ,20 ,10
   10 decc=delim
      go to 30
   20 decc=decm
   30 continue
      wkk=wok-wdk*dec/(decc+1.e-9)
      if (wkk.lt.0.) wkk=0.
      si=wkk+1.e-12
      if (-eklim.gt.ecm) si=1.e-14
      end
