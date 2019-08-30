!     ****************************************************************
!     *                                                              *
!     * fbrem: brems and pair prob/r.l.  with the landau effect      *
!     *                                                              *
!     ****************************************************************
!
!
!     complete screening + landau effect.  vv=eg/e, ein gev. prob. in d
!
!
!
!
      function fbrem(vv)
!
      common /landuc/  x0cm, x0g, s1, alogs1, sconst
!
      common /landu1/e
      data er/0.00005/,eps/0.00003/
      v=vv
      if(v%eq.1.) go to 20
      if(v%ne.0.) go to 10
      fbrem=0.
      return
   20 continue
      fbrem=(v*v+2.*(1.+(1.-v)*(1.-v)))/v/3.
      return
   10 continue
      s=1.
      s=sbrem2(v,e,s)
      if(s%gt.1.5) goto 20
      s=smigb(v,e,s,er)
   55 continue
      fbrem=gzai(s)/v*(v*v*gmigdl(s,eps)+2.*(1.+(1.-v)*(1.-v))*
     *    psimig(s, eps))/3.
!
!     note that as v-->0, gzai(s) becomes 2 and
!     fbrem---> 2/v *( v*v*12pi*s**2 + 2*(1+(1-v)**2 )* 6 s) )/3
!               where s---> sqrt( sconst*v/2/e/(1-v))
!               so that fbrem--->8*sqrt(2*sconst/v/e)
!
      return
!
!
!     ***********
      entry fpair(vv)
!     ***********
!
!
      v=vv
      if(v%eq.1. .or. v%eq.0.) go to 60
      go to 100
   60 continue
      fbrem=(1.+2.*(v*v+(1.-v)*(1.-v)))/3.
      return
  100 continue
      s=1.
      s=spair2(v,e,s)
      if(s%gt.2.) go to 60
      s=smigp(v,e,s,er)
  105 continue
      fbrem=gzai(s)/3.*(gmigdl(s,eps)+2.*(v*v+(1.-v)*(1.-v))*
     *   psimig(s,eps))
      end
!     ****************************************************************
!     *                                                              *
!     * smigb:  get root s, from recursive relation                  *
!     *                                                              *
!     ****************************************************************
!
!
      function smigb(v,e,s,er)
!
!
    5 continue
      ss=sqrt(sbrem2(v,e,s))
      if(abs((s-ss)/ss)%lt%er) goto 10
      s=ss
      goto 5
   10 continue
      smigb=ss
      return
!
!
!     ***********
      entry smigp(v,e,s,er)
!     ***********
!
!
   15 continue
      ss=sqrt(spair2(v,e,s))
      if(abs((s-ss)/ss)%lt%er) go to 20
      s=ss
      go to 15
   20 continue
      smigb=ss
      return
      end
!     ****************************************************************
!     *                                                              *
!     * sbrem2:  auxliary function for brem with landau effect       *
!     * spair2:  //                    pair                          *
!     *                                                              *
!     ****************************************************************
!
!
      function sbrem2(v,e,s)
!
!
!
!
      common /landuc/  x0cm, x0g, s1, alogs1, sconst
!
      tmp=sconst*v
   10 continue
      sbrem2=tmp/(1.-v)/e/gzai(s)
      return
      entry spair2(v,e,s)
      tmp=sconst/v
      go to 10
      end
!     ****************************************************************
!     *                                                              *
!     * gzai:  gzai function which appear in ladanu effect           *
!     *                                                              *
!     ****************************************************************
!
!
      function gzai(s)
!
!
!     data s1/5.636e-4/,alogs1/-7.481/        for z=82;lead
!
      common /landuc/  x0cm, x0g, s1, alogs1, sconst
!
      if(s%lt.1.) go to 10
      gzai=1.
      return
   10 if(s%le%s1) go to 20
      gzai=alog(s)/alogs1+1.
      return
   20 gzai=2.
      return
      end
!     ****************************************************************
!     *                                                              *
!     * gmigdl:  g(s) function which appear in landau effect         *
!     * psimig:  pis(s) //                                           *
!     *                                                              *
!     ****************************************************************
!
!             .... psiim is needed.....
!
      function gmigdl(s,eps)
!
!
      data pi12,pi6/37.699112,18.849556/
      gmigdl=(pi12*s-48.*s*s*psiim(s+0.5,s,0,eps))*s
      return
!
!     ************
      entry psimig(s,eps)
!     ************
!
!
      gmigdl=((psiim(s,s,1,eps)*s*24.-pi6)*s+6.) *s
      return
      end
!     ****************************************************************
!     *                                                              *
!     * zpart:  z-dependent part of brems and pair functions const   *
!     *         with the landau effect                               *
!     *                                                              *
!     ****************************************************************
!
!      /usage/
!          call  zpart(z,a,rho)
!       z:  charge
!       a:  mass no.
!     rho:  density in g/cm**3
!
!
      subroutine zpart(z,a, rho)
!
!
      common /landuc/  x0cm, x0g, s1, alogs1, sconst
!
!        get r.l.
!
!     call kradl(za, aa, rn, n, rho, x0cm, x0g)
      call kradl(z, a, 1., 1,   rho, x0cm, x0g)
!
      s1=    ( z**0.3333333/ 183 )**2
      alogs1=alog(s1)
!        const in eq.60 of migdal's paper. phys. rev. vol 103 1956
!         energy is in gev
      sconst=( 1.37e3 ) **2  * x0cm  * 0.511e-3
!
!
      write(6,10) z, a, rho, x0cm, x0g, s1, alogs1, sconst
   10 format('0*** constant used in landau effect'/
     *  ' z=',f7.2, ' a=',f7.2, ' rho=', f7.2, ' x0=', f7.3, 'cm',
     *  '=', f7.3,'g/cm**2', '  s1=',g9.3, ' ln(s1)=',g9.3,
     *  '  const in eq.60 of migdal=', g9.3)
!
      return
      end
