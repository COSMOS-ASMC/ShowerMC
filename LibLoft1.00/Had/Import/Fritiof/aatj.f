c***  header  ********************************************************
c                                                                    *
c  the lund monte carlo for jet fragmentation and e+e- physics       *
c  jetset version 6.3, october 1986                                  *
c  author: torbjorn sjostrand, department of theoretical physics,    *
c          university of lund, solvegatan 14a, s-223 62 lund, sweden *
c          bitnet/earn user thep node seldc51                        *
c  lushowC is written together with mats bengtsson, address as above  *
c  please report any errors to the author                            *
c                                                                    *
c*********************************************************************
c***  jetset version 6.3, general part  ******************************

      subroutine lupartC(ip,kf,pe,the,phi)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      external ludataC, luedatC, luhdatC

c...fill one particle (or jet if kf>=500)
      if(mst(19).ge.1) call lulistC(-1)
      if(ip.ge.mst(30)-5-mst(31)) mst(26)=1
      ir=max(ip,1)
      k(ir,1)=0
      k(ir,2)=kf
      if(mst(9).eq.0) p(ir,5)=ulmassC(1,kf)
      p(ir,4)=max(pe,p(ir,5))
      pa=sqrt(p(ir,4)**2-p(ir,5)**2)
      p(ir,1)=pa*sin(the)*cos(phi)
      p(ir,2)=pa*sin(the)*sin(phi)
      p(ir,3)=pa*cos(the)
      n=ir
      if(ip.eq.0) call luexecC

      return
      end

c*********************************************************************

      subroutine lu1jetC(ip,ifl,iflj,ifli,pe,the,phi)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)

c...fill one jet: leading quark, gluon or diquark
      if(iabs(ifl).gt.100) mst(26)=2
      ir=max(iabs(ip),1)
      call lupartC(ir,ifl+isign(500,2*ifl+1),pe,the,phi)
      if(ip.lt.0) k(n,1)=10000

      if(iflj.ne.0.or.ifli.ne.0) then
c...extra line for order in diquark and/or last quark in hadron jet
        if(iabs(iflj).gt.10.or.iabs(ifli).gt.10.or.(iflj.ne.0.and.
     &  iabs(ifl).lt.10).or.ifl*iflj.lt.0.or.(ifli.ne.0.and.ifl.eq.0)
     &  .or. ifl*(10-iabs(ifl))*ifli.gt.0) mst(26)=2
        n=ir+1
        k(n,1)=60000+ir
        k(n,2)=isign(600+10*iabs(iflj)+iabs(ifli),iflj+ifli)
        do 100 j=1,5
  100   p(n,j)=0.
      endif
      if(ip.eq.0) call luexecC

      return
      end

c*********************************************************************

      subroutine lu2jetC(ip,ifl1,ifl2,ecm)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)

c...flavour checks, fill two jets in cm frame, ifl1 along +z axis
      if(iabs(ifl1).gt.100.or.iabs(ifl2).gt.100.or.(ifl1.eq.0.and.
     &ifl2.ne.0).or.(ifl1.ne.0.and.ifl2.eq.0).or.ifl1*(10-
     &iabs(ifl1))*ifl2*(10-iabs(ifl2)).gt.0) mst(26)=2
      ir=max(iabs(ip),1)
      if(mst(9).eq.0) p(ir,5)=ulmassC(2,ifl1)
      if(mst(9).eq.0) p(ir+1,5)=ulmassC(2,ifl2)
      if(ecm.le.p(ir,5)+p(ir+1,5)) mst(26)=4
      pe1=0.5*(ecm+(p(ir,5)**2-p(ir+1,5)**2)/ecm)
      mst(9)=mst(9)+1
      call lupartC(ir,ifl1+isign(500,2*ifl1+1),pe1,0.,0.)
      k(n,1)=10000
      call lupartC(ir+1,ifl2+isign(500,2*ifl2+1),ecm-pe1,par(71),0.)
      mst(9)=mst(9)-1
      if(ip.eq.0) call luexecC

c...rearrange jets to prepare for shower evolution (optional)
      if(ip.lt.0) then
        k(ir+2,1)=k(ir+1,1)
        k(ir+2,2)=k(ir+1,2)
        do 100 j=1,5
  100   p(ir+2,j)=p(ir+1,j)
        do 110 i=ir+1,ir+3,2
        k(i,1)=70000+i-1
        k(i,2)=1000+i-1
        p(i,3)=0.
        p(i,4)=0.
  110   p(i,5)=0.
        p(ir+1,1)=ir+2
        p(ir+1,2)=ir+2
        p(ir+3,1)=ir
        p(ir+3,2)=ir
        n=n+2
      endif

      return
      end

c*********************************************************************

      subroutine lu3jetC(ip,ifl1,ifl3,ecm,x1,x3)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)

c...flavour checks, calculate masses and momenta
      if(iabs(ifl1).gt.100.or.iabs(ifl3).gt.100.or.(ifl1.eq.0.and.
     &ifl3.ne.0).or.(ifl1.ne.0.and.ifl3.eq.0).or.ifl1*(10-
     &iabs(ifl1))*ifl3*(10-iabs(ifl3)).gt.0) mst(26)=2
      ir=max(ip,1)
      if(mst(9).eq.0) p(ir,5)=ulmassC(2,ifl1)
      if(mst(9).eq.0) p(ir+1,5)=ulmassC(2,0)
      if(mst(9).eq.0) p(ir+2,5)=ulmassC(2,ifl3)
      if(0.5*x1*ecm.le.p(ir,5).or.0.5*(2.-x1-x3)*ecm.le.p(ir+1,5).or.
     &0.5*x3*ecm.le.p(ir+2,5)) mst(26)=4
      pa1=sqrt((0.5*x1*ecm)**2-p(ir,5)**2)
      pa2=sqrt((0.5*(2.-x1-x3)*ecm)**2-p(ir+1,5)**2)
      pa3=sqrt((0.5*x3*ecm)**2-p(ir+2,5)**2)

c...fill three jets in cm frame, ifl1 along +z axis, ifl3 in xz
c...plane with x>0; ifl2 is automatically gluon
      cthe2=(pa3**2-pa1**2-pa2**2)/(2.*pa1*pa2)
      if(abs(cthe2).ge.1.001) mst(26)=4
      if(abs(cthe2).le.1.001) cthe2=max(-1.,min(1.,cthe2))
      the2=-acos(cthe2)
      cthe3=(pa2**2-pa1**2-pa3**2)/(2.*pa1*pa3)
      if(abs(cthe3).ge.1.001) mst(26)=4
      if(abs(cthe3).le.1.001) cthe3=max(-1.,min(1.,cthe3))
      the3=acos(cthe3)
      mst(9)=mst(9)+1
      call lupartC(ir,ifl1+isign(500,2*ifl1+1),0.5*x1*ecm,0.,0.)
      k(n,1)=10000
      call lupartC(ir+1,500,0.5*(2.-x1-x3)*ecm,the2,0.)
      k(n,1)=10000
      call lupartC(ir+2,ifl3+isign(500,2*ifl3+1),0.5*x3*ecm,the3,0.)
      mst(9)=mst(9)-1
      if(ip.eq.0) call luexecC

c...rearrange jets to prepare for shower evolution (optional)
      if(ip.lt.0) then
        do 100 i=2,1,-1
        k(ir+2*i,1)=k(ir+i,1)
        k(ir+2*i,2)=k(ir+i,2)
        do 100 j=1,5
  100   p(ir+2*i,j)=p(ir+i,j)
        do 110 i=ir+1,ir+5,2
        k(i,1)=70000+i-1
        k(i,2)=1000+i-1
        p(i,3)=0.
        p(i,4)=0.
  110   p(i,5)=0.
        isg=1
        if((ifl1.lt.0.and.ifl1.gt.-10).or.ifl1.gt.10) isg=2
        p(ir+1,isg)=ir+2
        p(ir+1,3-isg)=ir+4
        p(ir+3,isg)=ir+4
        p(ir+3,3-isg)=ir
        p(ir+5,isg)=ir
        p(ir+5,3-isg)=ir+2
        n=n+3
      endif

      return
      end

c*********************************************************************

      subroutine lu4jetC(ip,ifl1,ifl2,ifl3,ifl4,ecm,x1,x2,x4,x12,x14)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)

c...flavour checks, calculate masses and momenta
      if(iabs(ifl1).gt.100.or.iabs(ifl2).gt.100.or.iabs(ifl3).gt.100
     & .or. iabs(ifl4).gt.100) mst(26)=2
      if(ifl2.eq.0.and.(ifl3.ne.0.or.(ifl1.eq.0.and.ifl4.ne.0).or.
     &(ifl1.ne.0.and.ifl4.eq.0).or.ifl1*(10-iabs(ifl1))*ifl4*
     &(10-iabs(ifl4)).gt.0)) mst(26)=2
      if(ifl2.ne.0.and.(ifl1.eq.0.or.ifl1*(10-iabs(ifl1))*ifl2*
     &(10-iabs(ifl2)).gt.0.or.(ifl3.eq.0.and.ifl4.ne.0).or.
     &(ifl3.ne.0.and.ifl4.eq.0).or.ifl3*(10-iabs(ifl3))*ifl4*
     &(10-iabs(ifl4)).gt.0)) mst(26)=2
      ir=max(ip,1)
      if(mst(9).eq.0) p(ir,5)=ulmassC(2,ifl1)
      if(mst(9).eq.0) p(ir+1,5)=ulmassC(2,ifl2)
      if(mst(9).eq.0) p(ir+2,5)=ulmassC(2,ifl3)
      if(mst(9).eq.0) p(ir+3,5)=ulmassC(2,ifl4)
      if(0.5*x1*ecm.le.p(ir,5).or.0.5*x2*ecm.le.p(ir+1,5).or.
     &0.5*(2.-x1-x2-x4)*ecm.le.p(ir+2,5).or.0.5*x4*ecm.le.
     &p(ir+3,5)) mst(26)=4
      pa1=sqrt((0.5*x1*ecm)**2-p(ir,5)**2)
      pa2=sqrt((0.5*x2*ecm)**2-p(ir+1,5)**2)
      pa3=sqrt((0.5*(2.-x1-x2-x4)*ecm)**2-p(ir+2,5)**2)
      pa4=sqrt((0.5*x4*ecm)**2-p(ir+3,5)**2)

c...kinematics for four jets in cm frame, ifl1 along +z axis,
c...ifl4 in xz plane with x>0, ifl2 with y>0 and y<0 equally often
      x24=x1+x2+x4-1.-x12-x14+(p(ir+2,5)**2-p(ir,5)**2-
     &p(ir+1,5)**2-p(ir+3,5)**2)/ecm**2
      cthe4=(x1*x4-2.*x14)*ecm**2/(4.*pa1*pa4)
      if(abs(cthe4).ge.1.002) mst(26)=4
      if(abs(cthe4).le.1.002) cthe4=max(-1.,min(1.,cthe4))
      the4=acos(cthe4)
      cthe2=(x1*x2-2.*x12)*ecm**2/(4.*pa1*pa2)
      if(abs(cthe2).ge.1.002) mst(26)=4
      if(abs(cthe2).le.1.002) cthe2=max(-1.,min(1.,cthe2))
      the2=acos(cthe2)
      cthe3=-(pa1+pa2*cthe2+pa4*cthe4)/pa3
      if(abs(cthe3).ge.1.002) mst(26)=4
      if(abs(cthe3).le.1.002) cthe3=max(-1.,min(1.,cthe3))
      the3=acos(cthe3)
      sgn=(-1.)**int(rluC(0)+0.5)
      cphi2=((x2*x4-2.*x24)*ecm**2-4.*pa2*cthe2*pa4*cthe4)/
     &(4.*pa2*sin(the2)*pa4*sin(the4))
      if(abs(cphi2).ge.1.05) mst(26)=4
      if(abs(cphi2).le.1.05) cphi2=max(-1.,min(1.,cphi2))
      phi2=sgn*acos(cphi2)
      cphi3=-(pa2*sin(the2)*cphi2+pa4*sin(the4))/(pa3*sin(the3))
      if(abs(cphi3).ge.1.05) mst(26)=4
      if(abs(cphi3).le.1.05) cphi3=max(-1.,min(1.,cphi3))
      phi3=-sgn*acos(cphi3)

c...fill jets, two separate systems if ifl2 not gluon
      mst(9)=mst(9)+1
      call lupartC(ir,ifl1+isign(500,2*ifl1+1),0.5*x1*ecm,0.,0.)
      k(n,1)=10000
      call lupartC(ir+1,ifl2+isign(500,2*ifl2+1),0.5*x2*ecm,the2,phi2)
      if(ifl2.eq.0) k(n,1)=10000
      call lupartC(ir+2,ifl3+isign(500,2*ifl3+1),0.5*(2.-x1-x2-x4)*ecm,
     &the3,phi3)
      k(n,1)=10000
      call lupartC(ir+3,ifl4+isign(500,2*ifl4+1),0.5*x4*ecm,the4,0.)
      mst(9)=mst(9)-1
      if(ip.eq.0) call luexecC

c...rearrange jets to prepare for shower evolution (optional)
      if(ip.lt.0) then
        do 100 i=3,1,-1
        k(ir+2*i,1)=k(ir+i,1)
        k(ir+2*i,2)=k(ir+i,2)
        do 100 j=1,5
  100   p(ir+2*i,j)=p(ir+i,j)
        do 110 i=ir+1,ir+7,2
        k(i,1)=70000+i-1
        k(i,2)=1000+i-1
        p(i,3)=0.
        p(i,4)=0.
  110   p(i,5)=0.
        if(ifl2.eq.0) then
          isg=1
          if((ifl1.lt.0.and.ifl1.gt.-10).or.ifl1.gt.10) isg=2
          do 120 i=1,4
          p(ir+2*i-1,isg)=ir+2*i
  120     p(ir+2*i-1,3-isg)=ir+2*i-4
          p(ir+1,3-isg)=ir+6
          p(ir+7,isg)=ir
        else
          do 130 j=1,2
          p(ir+1,j)=ir+2
          p(ir+3,j)=ir
          p(ir+5,j)=ir+6
  130     p(ir+7,j)=ir+4
        endif
        n=n+4
      endif

      return
      end

c*********************************************************************

      function kluC(i,j)
      common/lujetsC/n,k(2000,2),p(2000,5)

      kluC=0
      if(i.lt.0.or.j.le.0) return
c...number of lines, number of stable particles/jets, total charge,
c...number of jets
      if(i.eq.0.and.j.le.1) then
        kluC=n
      elseif(i.eq.0) then
        do 100 i1=1,n
        if(j.eq.2.and.k(i1,1).lt.20000) kluC=kluC+1
        if(j.eq.3.and.k(i1,1).lt.20000) kluC=kluC+luchgeC(k(i1,2))
        if(j.ne.3.or.k(i1,1)/10000.ne.6) goto 100
        if(k(i1-1,1).lt.20000) kluC=kluC+luchgeC(k(i1,2))
  100   if(j.eq.4.and.k(i1,1).lt.40000.and.iabs(k(i1,2)).ge.500) kluC=
     &  kluC+1

c...direct readout of k matrix or charge
      elseif(j.le.2) then
        kluC=k(i,j)
      elseif(j.le.4) then
        if(j.eq.3) kluC=luchgeC(k(i,2))

c...particle history: parent, generation, ancestor, rank
        if(j.eq.4) kluC=mod(k(i,1),10000)
      elseif(j.le.7) then
        i2=i
        i1=i
  110   kluC=kluC+1
        i3=i2
        i2=i1
        i1=mod(k(i1,1),10000)
        if(i1.gt.0.and.k(i1,1).lt.40000) goto 110
        if(j.eq.6) kluC=i2
        if(j.eq.7) then
          kluC=0
          do 120 i1=i2+1,i3
  120     if(mod(k(i1,1),10000).eq.i2.and.k(i1,1).lt.40000) kluC=kluC+1
        endif

c...particle code or ifl jet code, else 0 or 1000
      elseif(j.le.9) then
        if(j.eq.8.and.k(i,1).lt.60000.and.iabs(k(i,2)).lt.500) kluC=
     &  k(i,2)
        if(j.eq.9) kluC=1000
        if(j.eq.9.and.k(i,1).lt.60000.and.iabs(k(i,2)).ge.500) kluC=
     &  mod(k(i,2),500)

c...particle or jet code after cuts, else 0
      elseif(j.le.13) then
        if(k(i,1).lt.60000) kluC=k(i,2)
        if(j.ge.11.and.k(i,1).ge.20000) kluC=0
        kfa=iabs(k(i,2))
        if(j.ge.12.and.(kfa.eq.8.or.kfa.eq.10.or.kfa.eq.12.or.kfa.eq.
     &  14)) kluC=0
        if(j.ge.13.and.luchgeC(kfa).eq.0) kluC=0

c...heaviest flavour in hadron, 0 for non-hadron
      elseif(j.eq.14) then
        call luiflvC(k(i,2),ifla,iflb,iflc,ksp)
        if(ksp.ge.0) kluC=ifla
      endif

      return
      end

c*********************************************************************

      function pluC(i,j)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      dimension psum(4)

      pluC=0.
      if(i.lt.0.or.j.le.0.or.(i.eq.0.and.j.gt.6)) return
c...sum of momenta or charges (if i=0) or direct readout of p matrix
      if(i.eq.0.and.j.le.4) then
        do 100 i1=1,n
  100   if(k(i1,1).lt.20000) pluC=pluC+p(i1,j)
      elseif(i.eq.0.and.j.eq.5) then
        do 110 j1=1,4
        psum(j1)=0.
        do 110 i1=1,n
  110   if(k(i1,1).lt.20000) psum(j1)=psum(j1)+p(i1,j1)
        pluC=sqrt(max(0.,psum(4)**2-psum(1)**2-psum(2)**2-psum(3)**2))
      elseif(i.eq.0) then
        do 120 i1=1,n
        if(k(i1,1)/10000.ne.6) goto 120
        if(k(i1-1,1).lt.20000) pluC=pluC+luchgeC(k(i1,2))/3.
  120   if(k(i1,1).lt.20000) pluC=pluC+luchgeC(k(i1,2))/3.
      elseif(j.le.5) then
        pluC=p(i,j)

c...charge, total momentum, transverse momentum, transverse mass
      elseif(j.le.12) then
        if(j.eq.6) pluC=luchgeC(k(i,2))/3.
        if(j.eq.7.or.j.eq.8) pluC=p(i,1)**2+p(i,2)**2+p(i,3)**2
        if(j.eq.9.or.j.eq.10) pluC=p(i,1)**2+p(i,2)**2
        if(j.eq.11.or.j.eq.12) pluC=p(i,5)**2+p(i,1)**2+p(i,2)**2
        if(j.eq.8.or.j.eq.10.or.j.eq.12) pluC=sqrt(pluC)

c...theta and phi in radians or degrees
      elseif(j.le.16) then
        if(j.le.14) pluC=ulanglC(p(i,3),sqrt(p(i,1)**2+p(i,2)**2))
        if(j.ge.15) pluC=ulanglC(p(i,1),p(i,2))
        if(j.eq.14.or.j.eq.16) pluC=pluC*180./par(71)

c...true rapidity, rapidity with pion mass, pseudorapidity
      elseif(j.le.19) then
        pmr=0.
        if(j.eq.17) pmr=p(i,5)
        if(j.eq.18) pmr=ulmassC(0,17)
        pr=max(1e-20,pmr**2+p(i,1)**2+p(i,2)**2)
        pluC=sign(alog(min((sqrt(pr+p(i,3)**2)+abs(p(i,3)))/sqrt(pr),
     &  1e20)),p(i,3))

c...energy and momentum fractions (only to be used in cm frame)
      elseif(j.le.25) then
        if(j.eq.20) pluC=2.*sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)/par(75)
        if(j.eq.21) pluC=2.*p(i,3)/par(75)
        if(j.eq.22) pluC=2.*sqrt(p(i,1)**2+p(i,2)**2)/par(75)
        if(j.eq.23) pluC=2.*p(i,4)/par(75)
        if(j.eq.24) pluC=(p(i,4)+p(i,3))/par(75)
        if(j.eq.25) pluC=(p(i,4)-p(i,3))/par(75)
      endif

      return
      end

c*********************************************************************

      subroutine luroboC(the,phi,bex,bey,bez)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      dimension rot(3,3),pv(3)
      double precision dp(4),dbex,dbey,dbez,dga,dbep,dgabep

      imax=n
      if(mst(2).gt.0) imax=mst(2)
      if(the**2+phi**2.gt.1e-20) then
c...rotate (typically from z axis to direction theta,phi)
        rot(1,1)=cos(the)*cos(phi)
        rot(1,2)=-sin(phi)
        rot(1,3)=sin(the)*cos(phi)
        rot(2,1)=cos(the)*sin(phi)
        rot(2,2)=cos(phi)
        rot(2,3)=sin(the)*sin(phi)
        rot(3,1)=-sin(the)
        rot(3,2)=0.
        rot(3,3)=cos(the)
        do 120 i=max(1,mst(1)),imax
        if(mod(k(i,1)/10000,10).ge.6) goto 120
        do 100 j=1,3
  100   pv(j)=p(i,j)
        do 110 j=1,3
  110   p(i,j)=rot(j,1)*pv(1)+rot(j,2)*pv(2)+rot(j,3)*pv(3)
  120   continue
      endif

      if(bex**2+bey**2+bez**2.gt.1e-20) then
c...lorentz boost (typically from rest to momentum/energy=beta)
        dbex=bex
        dbey=bey
        dbez=bez
        dga=1d0/dsqrt(1d0-dbex**2-dbey**2-dbez**2)
        do 140 i=max(1,mst(1)),imax
        if(mod(k(i,1)/10000,10).ge.6) goto 140
        do 130 j=1,4
  130   dp(j)=p(i,j)
        dbep=dbex*dp(1)+dbey*dp(2)+dbez*dp(3)
        dgabep=dga*(dga*dbep/(1d0+dga)+dp(4))
        p(i,1)=dp(1)+dgabep*dbex
        p(i,2)=dp(2)+dgabep*dbey
        p(i,3)=dp(3)+dgabep*dbez
        p(i,4)=dga*(dp(4)+dbep)
  140   continue
      endif

      return
      end

c*********************************************************************

      subroutine lueditC(medit)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludat2C/ktyp(120),pmas(120),pwid(60),kfr(80),cfr(40)

      if(medit.ge.0.and.medit.le.3) then
c...throw away unwanted jets and particles
        imax=n
        if(mst(2).gt.0) imax=mst(2)
        mnot=0
        i1=max(1,mst(1))-1
        do 120 i=max(1,mst(1)),imax
        if(mnot.eq.1.and.k(i,1)/20000.eq.3) goto 100
        mnot=0
        if(k(i,1).ge.40000) goto 120
        if(medit.ge.1.and.k(i,1).ge.20000) goto 120
        kfa=iabs(k(i,2))
        if(medit.ge.2.and.(kfa.eq.8.or.kfa.eq.10.or.kfa.eq.12.or.
     &  kfa.eq.14)) goto 120
        if(medit.ge.3.and.kfa.le.499.and.luchgeC(kfa).eq.0) goto 120
        if(kfa.ge.500) mnot=1

c...pack remaining jets and particles, origin no longer known
  100   i1=i1+1
        k(i1,1)=10000*(k(i,1)/10000)
        k(i1,2)=k(i,2)
        do 110 j=1,5
  110   p(i1,j)=p(i,j)
  120   continue
        n=i1

      elseif(medit.eq.-1) then
c...save top entries at bottom of lujetsC
        if(2*n.ge.mst(30)) then
          mst(26)=1
          return
        endif
        do 130 i=1,n
        k(mst(30)-i,1)=k(i,1)
        k(mst(30)-i,2)=k(i,2)
        do 130 j=1,5
  130   p(mst(30)-i,j)=p(i,j)
        mst(31)=n

      elseif(medit.eq.-2) then
c...restore bottom entries of lujetsC to top
        do 140 i=1,mst(31)
        k(i,1)=k(mst(30)-i,1)
        k(i,2)=k(mst(30)-i,2)
        do 140 j=1,5
  140   p(i,j)=p(mst(30)-i,j)
        n=mst(31)

      elseif(medit.eq.-3) then
c...mark primary entries in top of lujetsC as untreated
        i1=0
        do 150 i=1,n
        kh=mod(k(i,1),10000)
        if(kh.ge.1) then
          if(k(kh,1)/20000.eq.2) kh=0
        endif
        if(k(i,1).ge.60000) kh=0
        if(kh.ne.0) goto 160
        i1=i1+1
  150   if(k(i,1)/20000.eq.1) k(i,1)=k(i,1)-20000
  160   n=i1
      endif

      return
      end

c*********************************************************************

      subroutine lulistC(mlist)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludat2C/ktyp(120),pmas(120),pwid(60),kfr(80),cfr(40)
      common/ludat3C/dpar(20),idb(120),cbr(400),kdp(1600)
      common/ludat4C/chag(50),chaf(100)
      character chag*4,chaf*4,chap*8,chan*8,chad(4)*8
      dimension ps(6)

      if((mlist.ge.0.and.mlist.le.2).or.(mlist.ge.10.and.mlist.le.12))
     &then
c...list event data
        if(mlist.le.2) write(mst(20),1000)
        if(mlist.ge.10) write(mst(20),1100)
        imax=n
        if(mst(2).gt.0) imax=mst(2)
        do 100 i=max(1,mst(1)),imax
        call lunameC(k(i,2),chap)
        mlc=0
        if(k(i,1)/20000.eq.1) mlc=1
        if(mlc.eq.1.and.iabs(k(i,2)).ge.500) mlc=2
        if(k(i,1)/20000.eq.2) mlc=k(i,1)/10000-1
        if(mlc.ne.0) chap(8:8)=chag(36)(mlc:mlc)
        if(k(i,1).ge.70000) mlc=10
        if(mlist.le.2.and.mlc.lt.10) write(mst(20),1200) i,
     &  mod(k(i,1),10000),chap,(p(i,j),j=1,5)
        if(mlist.ge.10.and.mlc.lt.10) write(mst(20),1300) i,k(i,1),
     &  k(i,2),chap,(p(i,j),j=1,5)
        if(mlist.le.2.and.mlc.eq.10) write(mst(20),1400) i,k(i,1),
     &  k(i,2),(p(i,j),j=1,5)
  100   if(mlist.ge.10.and.mlc.eq.10) write(mst(20),1500) i,k(i,1),
     &  k(i,2),(p(i,j),j=1,5)

c...sum of charges and momenta or extra lines after particles
        if(mlist.eq.1.or.mlist.eq.11) then
          do 110 j=1,6
  110     ps(j)=pluC(0,j)
          if(mlist.eq.1) write(mst(20),1600) ps(6),(ps(j),j=1,5)
          if(mlist.eq.11) write(mst(20),1700) ps(6),(ps(j),j=1,5)
        elseif(mlist.eq.2.or.mlist.eq.12) then
          do 120 i=n+1,n+mst(3)
          if(mlist.eq.2) write(mst(20),1400) i,k(i,1),k(i,2),
     &    (p(i,j),j=1,5)
  120     if(mlist.eq.12) write(mst(20),1500) i,k(i,1),k(i,2),
     &    (p(i,j),j=1,5)
        endif

      elseif(mlist.eq.3) then
c...list particle data table
        write(mst(20),1800)
        kf=max(1,mst(1))-1
  130   kf=kf+1
        write(mst(20),1900)

c...particle number, name, type, mass, width
  140   call lunameC(kf,chap)
        call lunameC(-kf,chan)
        kfa=kf
        if(kf.gt.100) call luiflvC(kf,ifla,iflb,iflc,ksp)
        if(kf.gt.100) kfa=100+ifla
        pm=ulmassC(0,kf)
        kty=ktyp(kfa)
        if(kty.lt.10) write(mst(20),2000) kf,chap,chan,kty,pm
        if(kty.ge.10) write(mst(20),2000) kf,chap,chan,kty,pm,
     &  pwid(2*(kty/10)-1),pwid(2*(kty/10))

        if(kf.gt.100.and.(mst(2).le.0.or.kf.lt.mst(2))) then
c...for heavy hadrons decay data only by group
          call luiflvC(kf+1-50*(kf/392),ifla1,iflb1,iflc1,ksp1)
          if(ifla1.eq.ifla) kf=kf+1
          if(ifla1.eq.ifla) goto 140
          kfa=76+5*ifla+ksp
        endif

c...particle decay: channel number, matrix element, branching
c...ratio, decay products
        if(idb(kfa).eq.0) goto 170
        idc=idb(kfa)-1
  150   idc=idc+1
        mmat=iabs(kdp(4*idc-3))/1000
        if(idc.eq.idb(kfa)) br=100.*cbr(idc)
        if(idc.ne.idb(kfa)) br=100.*(cbr(idc)-cbr(idc-1))
        do 160 j=1,4
  160   call lunameC(mod(kdp(4*idc-4+j),1000),chad(j))
        write(mst(20),2100) idc,mmat,br,(chad(j),j=1,4)
        if(cbr(idc).le.0.99999) goto 150
  170   if((mst(2).le.0.and.kf.lt.392).or.(mst(2).gt.0.and.kf.lt.
     &  mst(2))) goto 130

      elseif(mlist.eq.4) then
c...list parton/jet data table
        write(mst(20),2200)
        ifl=max(0,mst(1))-1
  180   ifl=ifl+1
        if(ifl.gt.0.and.mod(ifl-1,10).ge.8) goto 180
        call lunameC(ifl+500,chap)
        call lunameC(-ifl-500,chan)
        pmc=ulmassC(2,ifl)
        pma=ulmassC(3,ifl)
        kty=ktyp(100+max(ifl/10,mod(ifl,10)))
        if(kty.lt.10) write(mst(20),2300) ifl+500,ifl,chap,chan,kty,
     &  pmc,pma
        if(kty.ge.10) write(mst(20),2300) ifl+500,ifl,chap,chan,kty,
     &  pmc,pma,pwid(2*(kty/10)-1),pwid(2*(kty/10))
        if((mst(2).le.0.and.ifl.lt.88).or.(mst(2).gt.0.and.ifl.lt.
     &  mst(2).and.(mod(ifl,10).ne.8.or.mst(2)-ifl.ge.3))) goto 180

      elseif(mlist.eq.5) then
c...list parameter value table
        write(mst(20),2400)
        do 190 l=1,20
  190   write(mst(20),2500) l,mst(l),mst(l+20),par(l),par(l+20),
     &  par(l+40),par(l+60),dpar(l)

      elseif(mlist.eq.-1) then
c...initialization printout (monte carlo version number and date)
        write(mst(20),2600)
        mst(19)=0
      endif

c...format statements for output on unit mst(20) (default 6)
 1000 format(///20x,'event listing'//5x,'i     ori   part/jet',7x,
     &'px',9x,'py',9x,'pz',9x,'e',10x,'m'/)
 1100 format(///20x,'event listing (extended)'//5x,'i  k(i,1)  k(i,2)',
     &3x,'part/jet',7x,'p(i,1)',7x,'p(i,2)',7x,'p(i,3)',7x,'p(i,4)',
     &7x,'p(i,5)'/)
 1200 format(2x,i4,1x,i7,3x,a8,5(1x,f10.3))
 1300 format(2x,i4,2(1x,i7),3x,a8,5(1x,f12.5))
 1400 format(2x,i4,1x,i7,4x,i7,5(1x,f10.3))
 1500 format(2x,i4,2(1x,i7),11x,5(1x,f12.5))
 1600 format(10x,'sum:',6(1x,f10.3))
 1700 format(16x,'sum:',6(1x,f12.5))
 1800 format(///20x,'particle data table'//4x,'kf    particle   ',
     &'antipart  ktyp         mass       width       w-cut'/18x,
     &'idc    mat    b.r.   decay products')
 1900 format(10x)
 2000 format(1x,i5,4x,a8,3x,a8,1x,i5,1x,f12.5,1x,f11.5,1x,f11.5)
 2100 format(16x,i5,3x,'(',i2,')',1x,f7.1,4(3x,a8))
 2200 format(///20x,'parton/jet data table'//4x,'kf   ifl     parton',
     &'    antipar  ktyp    m-cons    m-c.a.     width     w-cut')
 2300 format(/1x,i5,1x,i5,4x,a8,3x,a8,1x,i4,4(1x,f9.3))
 2400 format(///20x,'parameter value table'//5x,'l',4x,'mst(l)',
     &3x,'&(l+20)',7x,'par(l)',6x,'&(l+20)',6x,'&(l+40)',6x,
     &'&(l+60)',6x,'dpar(l)'/)
 2500 format(1x,i5,2(1x,i9),5(1x,f12.4))
 2600 format(///20x,'the lund monte carlo - jetset version 6.3'/
     &          20x,'  last date of change:  17 october 1986  ')

      return
      end

c*********************************************************************

      subroutine luupdaC(mupda,lfn)
      common/ludat1C/mst(40),par(80)
      common/ludat2C/ktyp(120),pmas(120),pwid(60),kfr(80),cfr(40)
      common/ludat3C/dpar(20),idb(120),cbr(400),kdp(1600)
      common/ludat4C/chag(50),chaf(100)
      character chag*4,chaf*4,cli*72,cut*12,cwr*12,csa*12,cre*12

      if(mupda.eq.1) then
c...write information on file for editing
        do 120 kf=1,120
        if(idb(kf).eq.0) then
          ndc=0
        else
          idc=idb(kf)-1
  100     idc=idc+1
          if(cbr(idc).le.0.99999) goto 100
          ndc=idc+1-idb(kf)
          if(kf.ge.2) then
            if(idb(kf).eq.idb(kf-1)) ndc=-1
          endif
        endif
        kty=ktyp(kf)-10*(ktyp(kf)/10)
        pwi=0.
        pcu=0.
        if(ktyp(kf).ge.10) pwi=pwid(2*(ktyp(kf)/10)-1)
        if(ktyp(kf).ge.10) pcu=pwid(2*(ktyp(kf)/10))
        if(kf.le.100) write(lfn,1000) kf,ndc,kty,pmas(kf),pwi,pcu,
     &  chaf(kf)
        if(kf.gt.100) write(lfn,1000) kf,ndc,kty,pmas(kf),pwi,pcu
        do 110 idc=idb(kf),idb(kf)+ndc-1
        mmat=iabs(kdp(4*idc-3))/1000
  110   write(lfn,1100) cbr(idc),mmat,(mod(kdp(4*idc-4+j),1000),j=1,4)
  120   continue

      elseif(mupda.eq.2) then
c...read information from editing
        do 130 i=1,60
  130   pwid(i)=0.
        do 140 i=1,400
  140   cbr(i)=0.
        do 150 i=1,1600
  150   kdp(i)=0
        iwis=0
        idbs=0
        do 170 kf=1,120
        if(kf.le.100) read(lfn,1000) kfa,ndc,ktyp(kf),pmas(kf),pwi,pcu,
     &  chaf(kf)
        if(kf.gt.100) read(lfn,1000) kfa,ndc,ktyp(kf),pmas(kf),pwi,pcu
        if(pwi.ge.0.0005) then
          pwid(2*iwis+1)=pwi
          pwid(2*iwis+2)=pcu
          iwis=iwis+1
          ktyp(kf)=ktyp(kf)+10*iwis
        endif
        if(ndc.eq.0) then
          idb(kf)=0
        elseif(ndc.eq.-1) then
          idb(kf)=idb(kf-1)
        else
          idb(kf)=idbs+1
          do 160 idc=idbs+1,idbs+ndc
          read(lfn,1100) cbr(idc),mmat,(kdp(4*idc-4+j),j=1,4)
  160     kdp(4*idc-3)=kdp(4*idc-3)+isign(1000*mmat,kdp(4*idc-3))
          idbs=idbs+ndc
        endif
  170   continue

      elseif(mupda.eq.3) then
c...write information for inclusion in program
        do 220 ic=1,12
        ne=120
        if(ic.eq.3) ne=60
        if(ic.eq.5.or.ic.eq.6) ne=200
        if(ic.ge.7.and.ic.le.11) ne=320
        if(ic.eq.12) ne=100
        cli=' '
        if(ic.eq.1) cli(7:16)='data ktyp/'
        if(ic.eq.2) cli(7:16)='data pmas/'
        if(ic.eq.3) cli(7:16)='data pwid/'
        if(ic.eq.4) cli(7:15)='data idb/'
        if(ic.eq.5) cli(7:28)='data (cbr(j),j=1,200)/'
        if(ic.eq.6) cli(7:30)='data (cbr(j),j=201,400)/'
        if(ic.eq.7) cli(7:28)='data (kdp(j),j=1,320)/'
        if(ic.eq.8) cli(7:30)='data (kdp(j),j=321,640)/'
        if(ic.eq.9) cli(7:30)='data (kdp(j),j=641,960)/'
        if(ic.eq.10) cli(7:31)='data (kdp(j),j=961,1280)/'
        if(ic.eq.11) cli(7:32)='data (kdp(j),j=1281,1600)/'
        if(ic.eq.12) cli(7:16)='data chaf/'
        lct=16
        if(ic.eq.4) lct=15
        if(ic.eq.5.or.ic.eq.7) lct=28
        if(ic.eq.6.or.ic.eq.8.or.ic.eq.9) lct=30
        if(ic.eq.10) lct=31
        if(ic.eq.11) lct=32
        csa='start'
        do 210 ie=1,ne
        if(ic.eq.1) write(cut,1200) ktyp(ie)
        if(ic.eq.2) write(cut,1300) pmas(ie)
        if(ic.eq.3) write(cut,1300) pwid(ie)
        if(ic.eq.4) write(cut,1200) idb(ie)
        if(ic.eq.5) write(cut,1300) cbr(ie)
        if(ic.eq.6) write(cut,1300) cbr(200+ie)
        if(ic.eq.7) write(cut,1200) kdp(ie)
        if(ic.eq.8) write(cut,1200) kdp(320+ie)
        if(ic.eq.9) write(cut,1200) kdp(640+ie)
        if(ic.eq.10) write(cut,1200) kdp(960+ie)
        if(ic.eq.11) write(cut,1200) kdp(1280+ie)
        if(ic.eq.12) cut=chaf(ie)
        cwr=' '
        la=1
        lb=1
        do 180 ll=1,12
        if(cut(13-ll:13-ll).ne.' ') la=13-ll
  180   if(cut(ll:ll).ne.' ') lb=ll
        lon=1+lb-la
        cwr(1:lon)=cut(la:lb)
        if(ic.eq.12) then
          do 190 ll=lon,1,-1
          if(cwr(ll:ll).eq.'''') then
            cwr=cwr(1:ll)//''''//cwr(ll+1:11)
            lon=lon+1
          endif
  190     continue
          cut=cwr
          cwr(1:lon+2)=''''//cut(1:lon)//''''
          lon=lon+2
        elseif(ic.eq.2.or.ic.eq.3.or.ic.eq.5.or.ic.eq.6) then
          lon=lon+1
  200     lon=lon-1
          if(cwr(lon:lon).eq.'0') goto 200
          if(lon.eq.1) cwr(1:2)='0.'
          if(lon.eq.1) lon=2
        endif
        if(cwr.ne.csa) then
          iag=1
          csa=cwr
        else
          lex=lon+1
          if(iag.ge.2) lex=lon+3
          if(iag.ge.10) lex=lon+4
          if(iag.ge.100) lex=lon+5
          lct=lct-lex
          iag=iag+1
          write(cre,1200) iag
          lex=1
          if(iag.ge.10) lex=2
          if(iag.ge.100) lex=3
          cut=cwr
          cwr(1:lex+1+lon)=cre(13-lex:12)//'*'//cut(1:lon)
          lon=lon+lex+1
        endif
        if(lct+lon.gt.70) then
          cli(lct+1:72)=' '
          write(lfn,1400) cli
          cli=' '
          cli(6:6)='&'
          lct=6
        endif
        cli(lct+1:lct+lon)=cwr(1:lon)
        lct=lct+lon+1
        if(ie.lt.ne) cli(lct:lct)=','
  210   if(ie.eq.ne) cli(lct:lct)='/'
  220   write(lfn,1400) cli
      endif

c...formats for reading and writing particle data
 1000 format(3i5,3f12.5,2x,a4)
 1100 format(5x,f12.5,5i5)
 1200 format(i12)
 1300 format(f12.5)
 1400 format(a72)

      return
      end

c*********************************************************************

      subroutine luexecC
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludat3C/dpar(20),idb(120),cbr(400),kdp(1600)
      dimension psum(2,5)

c...reset and initialize, sum up energy of original jets/particles
      if(mst(19).ge.1) call lulistC(-1)
      mst(2)=0
      mst(3)=0
      nerr=mst(24)
      mst(25)=0
      mst(32)=0
      mst(34)=mst(34)+1
      par(75)=0.
      do 100 i=1,n
  100 if(mod(k(i,1),10000).eq.0) par(75)=par(75)+p(i,4)

c...sum up momentum, energy and charge for starting entries
      do 110 i=1,2
      do 110 j=1,5
  110 psum(i,j)=0.
      icon=0
      do 130 i=1,n
      if(k(i,1).ge.20000) goto 130
      do 120 j=1,4
  120 psum(1,j)=psum(1,j)+p(i,j)
      psum(1,5)=psum(1,5)+luchgeC(k(i,2))
      if(i.eq.n) goto 130
      if(k(i+1,1)/10000.eq.6) psum(1,5)=psum(1,5)+luchgeC(k(i+1,2))
  130 continue

c...check and prepare system for subsequent fragmentation/decay
      call luprepC
      mst(1)=0
      if(mst(23).eq.1.and.mst(26).ne.0.and.mst(35).lt.5) then
        mst(35)=mst(35)+1
        write(mst(20),1000) mst(34)
        if(mst(26).eq.1) write(mst(20),1100)
        if(mst(26).eq.2) write(mst(20),1200)
        if(mst(26).eq.3) write(mst(20),1300)
        if(mst(26).eq.4) write(mst(20),1400)
        mst(26)=0
      endif

c...administrate jet fragmentation and particle decay chain
      ip=0
  140 ip=ip+1
      if(k(ip,1).ge.20000) then

      elseif(iabs(k(ip,2)).lt.500) then
c...particle decay if unstable
        kfa=iabs(k(ip,2))
        if(kfa.gt.100) call luiflvC(kfa,ifla,iflb,iflc,ksp)
        if(kfa.gt.100) kfa=76+5*ifla+ksp
        if(mst(7).ge.1.and.idb(kfa).ge.1) call ludecyC(ip)

      elseif(iabs(k(ip,2)).le.600) then
c...jet fragmentation: one jet or system
        mos=min(mst(5),2)
        if(mos.eq.2.and.mst(6).gt.0) mos=3
        if(mst(5).ge.1.and.k(ip,1).lt.10000) mos=2
        if(mst(7).ge.2.and.k(ip,1).ge.10000.and.n.gt.ip) then
          kh=mod(k(i,1),10000)
          if(k(ip+1,1).lt.10000.and.kh.gt.0.and.kh.lt.ip) then
            if(k(kh,1).lt.40000.and.iabs(k(kh,2)).lt.400) mos=1
          endif
        endif
        if(mos.eq.1) call lusysjC(ip)
        if(mos.eq.2) call luonejC(ip)
        if(mos.eq.3) call luconsC(ip)
        if(mos.eq.2) icon=1
        if(mos.eq.3.and.(mst(6).le.0.or.mod(mst(6),5).eq.0)) icon=1
      endif

c...error checks and printout
      if(n.ge.mst(30)-20-mst(31).and.ip.lt.n.and.mst(24).eq.nerr)
     &then
        mst(24)=mst(24)+1
        mst(25)=1
      endif
      if((mst(23).eq.1.and.mst(24).gt.nerr).or.(mst(23).ge.2.and.
     &mst(24).ge.2)) then
        write(mst(20),1500) mst(24),mst(34)
        if(mst(25).eq.1) write(mst(20),1100)
        if(mst(25).eq.2) write(mst(20),1200)
        if(mst(25).eq.3) write(mst(20),1300)
        if(mst(25).eq.4) write(mst(20),1600)
        if(mst(25).eq.5) write(mst(20),1700)
      endif
      if((mst(23).eq.1.and.mst(24).ge.5).or.(mst(23).ge.2.and.
     &mst(24).ge.2)) then
        write(mst(20),1800)
        mst(1)=0
        mst(2)=0
        call lulistC(11)
        stop
      elseif(mst(23).ge.1.and.mst(24).gt.nerr) then
        return
      endif
      if(ip.lt.n) goto 140

c...check that momentum, energy and charge were conserved
      do 160 i=1,n
      if(k(i,1).ge.20000) goto 160
      do 150 j=1,4
  150 psum(2,j)=psum(2,j)+p(i,j)
      psum(2,5)=psum(2,5)+luchgeC(k(i,2))
      if(i.eq.n) goto 160
      if(k(i+1,1)/10000.eq.6) psum(2,5)=psum(2,5)+luchgeC(k(i+1,2))
  160 continue
      pdev=(abs(psum(2,1)-psum(1,1))+abs(psum(2,2)-psum(1,2))+
     &abs(psum(2,3)-psum(1,3))+abs(psum(2,4)-psum(1,4)))/
     &(1.+abs(psum(2,4))+abs(psum(1,4)))
      if(icon.eq.0.and.(pdev.gt.par(74).or.abs(psum(2,5)-psum(1,5))
     &.gt.0.25)) then
        mst(24)=mst(24)+1
        mst(25)=6
      endif
      if((mst(23).eq.1.and.mst(24).gt.nerr).or.(mst(23).ge.2.and.
     &mst(24).ge.2)) then
        write(mst(20),1500) mst(24),mst(34)
        write(mst(20),1900) ((psum(i,j),j=1,4),psum(i,5)/3.,i=1,2)
      endif
      if((mst(23).eq.1.and.mst(24).ge.5).or.(mst(23).ge.2.and.
     &mst(24).ge.2)) then
        write(mst(20),1800)
        call lulistC(11)
        stop
      endif

c...format statements for error warnings
 1000 format(/5x,'warning! mst(26) flag was set at luexecC ',
     &'call no',i8,'; error type is')
 1100 format(5x,'1: not enough memory available in commonblock lujetsC')
 1200 format(5x,'2: unphysical flavour setup of jet system')
 1300 format(5x,'3: not enough energy available in jet system ',
     &'(string fragmentation)')
 1400 format(5x,'4: inconsistent kinematics for definition of jet ',
     &'configuration')
 1500 format(/5x,'warning! error no',i2,' has occured in luexecC ',
     &'call no',i8,'; error type is')
 1600 format(5x,'4: not enough energy available in jet system ',
     &'(independent fragmentation)')
 1700 format(5x,'5: no kinematically allowed decays are found for ',
     &'this particle')
 1800 format(5x,'execution will be stopped after printout of ',
     &'event listing')
 1900 format(5x,'6: momentum, energy and/or charge were not conserved'/
     &5x,'sum of',9x,'px',11x,'py',11x,'pz',11x,'e',8x,'charge'/
     &5x,'before',2x,4(1x,f12.5),1x,f8.2/5x,'after',3x,4(1x,f12.5),1x,
     &f8.2)

      return
      end

c*********************************************************************

      subroutine luprepC
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludat3C/dpar(20),idb(120),cbr(400),kdp(1600)
      dimension ps(5),pc(5),ue(3)

c...rearrange parton shower product listing along strings: begin loop
      ns=n
      do 120 iqg=1,2
      do 120 i=1,ns-1
      if(k(i+1,1)/10000.ne.7.or.k(i,1).ge.20000.or.iabs(k(i,2)).lt.
     &500.or.(iqg.eq.1.and.iabs(k(i,2)).eq.500)) goto 120

c...pick up loose string end, copy undecayed parton
      kcs=(3-isign(1,k(i,2)*(510-iabs(k(i,2)))))/2
      ia=i
      nl=0
  100 nl=nl+1
      if(nl.gt.2*ns) then
        mst(26)=2
        return
      endif
      if(k(ia,1).lt.20000) then
        n=n+1
        if(n.ge.mst(30)-5-mst(31)) then
          mst(26)=1
          return
        endif
        k(n,1)=10000
        if(nl.ge.2.and.iabs(k(ia,2)).gt.500) k(n,1)=0
        k(n,1)=k(n,1)+max(0,k(ia+1,2)-1000)
        k(n,2)=k(ia,2)
        do 110 j=1,5
  110   p(n,j)=p(ia,j)
        k(ia,1)=k(ia,1)+20000
        if(k(n,1).lt.10000) goto 120
      endif

c...go to next parton in colour space
      ib=ia
      if(p(ib+1,kcs+2).gt.0.5) then
        ia=nint(p(ib+1,kcs+2))
        p(ib+1,kcs+2)=-p(ib+1,kcs+2)
        mm=0
      else
        if(p(ib+1,kcs).lt.0.5) kcs=3-kcs
        ia=nint(p(ib+1,kcs))
        p(ib+1,kcs)=-p(ib+1,kcs)
        mm=1
      endif
      if(ia.le.0.or.ia.gt.min(ns,mst(30)-mst(31))) then
        mst(26)=2
        return
      endif
      if(nint(p(ia+1,1)).eq.ib.or.nint(p(ia+1,2)).eq.ib) then
        if(mm.eq.1) kcs=3-kcs
        if(nint(p(ia+1,kcs)).ne.ib) kcs=3-kcs
        p(ia+1,kcs)=-p(ia+1,kcs)
      else
        if(mm.eq.0) kcs=3-kcs
        if(nint(p(ia+1,kcs+2)).ne.ib) kcs=3-kcs
        p(ia+1,kcs+2)=-p(ia+1,kcs+2)
      endif
      if(ia.ne.i) goto 100
      k(n,1)=k(n,1)-10000
  120 continue

      if(mst(21).ge.1) then
c...delete unnecessary parton shower evolution information
        k(n+1,1)=0
        i1=0
        do 140 i=1,n
        ks=k(i,1)/10000
        if(ks.ge.7.or.(mst(21).ge.3.and.ks.ge.2.and.ks.le.5)) goto 140
        if(ks.ge.2.and.i.lt.n.and.k(i+1,1)/10000.eq.7) then
          if(mst(21).ge.2.and.ks.ne.6) goto 140
          if(ks.le.3.and.i.gt.mst(1)) goto 140
          if(ks.le.3) k(i,1)=40000
        endif
        i1=i1+1
        k(i1,1)=k(i,1)
        if(i.lt.n.and.k(i+1,1)/10000.eq.7) k(i1,1)=
     &  10000*(k(i1,1)/10000)+max(0,k(i+1,2)-1000)
        k(i1,2)=k(i,2)
        do 130 j=1,5
  130   p(i1,j)=p(i,j)
  140   continue
        n=i1
      endif

c...find lowest-mass colour singlet jet system, ok if above threshold
      if(mst(12).le.0) goto 310
      ns=n
  150 nsin=n-ns
      pdm=1.+par(22)
      ic=0
      do 200 i=1,ns
      if(k(i,1).ge.20000.or.iabs(k(i,2)).lt.500) goto 200
      if(k(i,1).ge.10000.and.ic.eq.0) then
        nsin=nsin+1
        ic=i
        do 160 j=1,4
  160   ps(j)=p(i,j)
        ps(5)=ulmassC(0,k(i,2))
      elseif(k(i,1).ge.10000) then
        do 170 j=1,4
  170   ps(j)=ps(j)+p(i,j)
      elseif(ic.ne.0) then
        do 180 j=1,4
  180   ps(j)=ps(j)+p(i,j)
        ps(5)=ps(5)+ulmassC(0,k(i,2))
        pd=sqrt(max(0.,ps(4)**2-ps(1)**2-ps(2)**2-ps(3)**2))-ps(5)
        if(pd.lt.pdm) then
          pdm=pd
          do 190 j=1,5
  190     pc(j)=ps(j)
          icl=ic
          icu=i
        endif
        ic=0
      endif
  200 continue
      if(pdm.ge.par(22).or.nsin.eq.1) goto 310

c...form two particles from flavours of lowest-mass system, if feasible
      pcm=sqrt(max(0.,pc(4)**2-pc(1)**2-pc(2)**2-pc(3)**2))
      k(n+1,1)=icl
      k(n+2,1)=icu
      if(k(icl+1,1)/10000 .eq. 6 .or. (icu.lt.n.and.k(icu+1,1)/10000
     & .eq.6)) then
        goto 310
      elseif(iabs(k(icl,2)).gt.500) then
        if(mod(k(icl,2),500)*mod(k(icu,2),500)*(510-iabs(k(icl,2)))*
     &  (510-iabs(k(icu,2))).ge.0) goto 310
  210   call luifldC(mod(k(icl,2),500),0,0,ifln,k(n+1,2))
        if(iabs(ifln).ge.100.or.(iabs(ifln).gt.10.and.iabs(k(icu,2))
     &  .gt.510)) goto 210
        call luifldC(mod(k(icu,2),500),0,-ifln,ifldmp,k(n+2,2))
      else
        if(iabs(k(icu,2)).ne.500) goto 310
  220   call luifldC(1+int((2.+par(2))*rluC(0)),0,0,ifln,kdump)
        if(iabs(ifln).ge.100) goto 220
        call luifldC(ifln,0,0,iflm,k(n+1,2))
        if(iabs(iflm).ge.100) goto 220
        call luifldC(-ifln,0,-iflm,ifldmp,k(n+2,2))
      endif
      p(n+1,5)=ulmassC(1,k(n+1,2))
      p(n+2,5)=ulmassC(1,k(n+2,2))
      if(p(n+1,5)+p(n+2,5)+dpar(14).ge.pcm) goto 260

c...perform two-particle decay of jet system, if possible
      if(pcm.ge.0.02*pc(4)) then
        pa=sqrt((pcm**2-(p(n+1,5)+p(n+2,5))**2)*(pcm**2-
     &  (p(n+1,5)-p(n+2,5))**2))/(2.*pcm)
        ue(3)=2.*rluC(0)-1.
        phi=par(72)*rluC(0)
        ue(1)=sqrt(1.-ue(3)**2)*cos(phi)
        ue(2)=sqrt(1.-ue(3)**2)*sin(phi)
        do 230 j=1,3
        p(n+1,j)=pa*ue(j)
  230   p(n+2,j)=-pa*ue(j)
        p(n+1,4)=sqrt(pa**2+p(n+1,5)**2)
        p(n+2,4)=sqrt(pa**2+p(n+2,5)**2)
        mst1s=mst(1)
        mst(1)=n+1
        n=n+2
        call luroboC(0.,0.,pc(1)/pc(4),pc(2)/pc(4),pc(3)/pc(4))
        mst(1)=mst1s
      else
        np=0
        do 240 i=icl,icu
  240   if(k(i,1).lt.20000) np=np+1
        ha=p(icl,4)*p(icu,4)-p(icl,1)*p(icu,1)-p(icl,2)*p(icu,2)-
     &  p(icl,3)*p(icu,3)
        if(np.ge.3.or.ha.le.1.25*p(icl,5)*p(icu,5)) goto 260
        hd1=0.5*(p(n+1,5)**2-p(icl,5)**2)
        hd2=0.5*(p(n+2,5)**2-p(icu,5)**2)
        hr=sqrt(max(0.,((ha-hd1-hd2)**2-(p(n+1,5)*p(n+2,5))**2)/
     &  (ha**2-(p(icl,5)*p(icu,5))**2)))-1.
        hc=p(icl,5)**2+2.*ha+p(icu,5)**2
        hk1=((p(icu,5)**2+ha)*hr+hd1-hd2)/hc
        hk2=((p(icl,5)**2+ha)*hr+hd2-hd1)/hc
        do 250 j=1,4
        p(n+1,j)=(1.+hk1)*p(icl,j)-hk2*p(icu,j)
  250   p(n+2,j)=(1.+hk2)*p(icu,j)-hk1*p(icl,j)
        n=n+2
      endif
      goto 290

c...else form one particle from the flavours available, if possible
  260 if(iabs(k(icl,2)).gt.510.and.iabs(k(icu,2)).gt.510) then
        goto 310
      elseif(iabs(k(icl,2)).gt.500) then
        call luifldC(mod(k(icl,2),500),0,mod(k(icu,2),500),
     &  ifldmp,k(n+1,2))
      else
        ifln=1+int((2.+par(2))*rluC(0))
        call luifldC(ifln,0,-ifln,ifldmp,k(n+1,2))
      endif
      p(n+1,5)=ulmassC(1,k(n+1,2))

c...find parton/particle which combines to largest extra mass
      ir=0
      ha=0.
      do 270 i=1,n
      if(k(i,1).ge.20000.or.(i.ge.icl.and.i.le.icu).or.
     &(iabs(k(i,2)).lt.500.and.i.le.ns)) goto 270
      pcr=pc(4)*p(i,4)-pc(1)*p(i,1)-pc(2)*p(i,2)-pc(3)*p(i,3)
      if(pcr.gt.ha) then
        ir=i
        ha=pcr
      endif
  270 continue

c...shuffle energy and momentum to put new particle on mass shell
      hb=pcm**2+ha
      hc=p(n+1,5)**2+ha
      hd=p(ir,5)**2+ha
      hk2=0.5*(hb*sqrt(((hb+hc)**2-4.*(hb+hd)*p(n+1,5)**2)/
     &(ha**2-(pcm*p(ir,5))**2))-(hb+hc))/(hb+hd)
      hk1=(0.5*(p(n+1,5)**2-pcm**2)+hd*hk2)/hb
      do 280 j=1,4
      p(n+1,j)=(1.+hk1)*pc(j)-hk2*p(ir,j)
  280 p(ir,j)=(1.+hk2)*p(ir,j)-hk1*pc(j)
      n=n+1

c...mark collapsed system, iterate
  290 do 300 i=icl,icu
  300 if(k(i,1).le.20000.and.iabs(k(i,2)).ge.500) k(i,1)=k(i,1)+20000
      if(n.lt.mst(30)-5-mst(31)) goto 150

c...check flavours and invariant masses in string systems
  310 np=0
      kfn=0
      kfs=0
      do 320 j=1,5
  320 ps(j)=0.
      do 350 i=1,n
      if(k(i,1).ge.20000.or.iabs(k(i,2)).lt.500) goto 350
      np=np+1
      if(iabs(k(i,2)).gt.500) then
        kfn=kfn+1
        kfs=kfs+isign(1,k(i,2)*(510-iabs(k(i,2))))
        if(n.gt.i.and.k(i+1,1)/10000.eq.6) kfs=kfs+isign(1,
     &  mod(k(i+1,2),10))
        ps(5)=ps(5)+ulmassC(0,k(i,2))
      endif
      do 330 j=1,4
  330 ps(j)=ps(j)+p(i,j)
      if(k(i,1).lt.10000) then
        if(np.ne.1.and.(kfn.eq.1.or.kfn.ge.3.or.kfs.ne.0)) mst(26)=2
        if(np.ne.1.and.ps(4)**2-ps(1)**2-ps(2)**2-ps(3)**2.lt.(par(22)+
     &  ps(5))**2) mst(26)=3
        np=0
        kfn=0
        kfs=0
        do 340 j=1,5
  340   ps(j)=0.
      endif
  350 continue

      return
      end

c*********************************************************************

      subroutine luconsC(ip)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      dimension ps2(4),nfl(3),ifet(3),iflf(3),te(3),td(3,3)
      double precision dps1(4),dp(4),dbe(3),dga,dbep,dgabep

c...reset counters, identify parton system, boost to cm frame
      ntry=0
  100 ntry=ntry+1
      if(ntry.gt.200) then
        mst(24)=mst(24)+1
        mst(25)=4
        if(mst(23).ge.1) return
      endif
      do 110 j=1,3
      nfl(j)=0
      ifet(j)=0
  110 iflf(j)=0
      if(ntry.eq.1) then
        do 120 j=1,4
  120   dps1(j)=0.
        in=ip-1
        njet=0
  130   in=in+1
        if(in.gt.min(n,mst(30)-mst(31))) then
          mst(24)=mst(24)+1
          mst(25)=2
          if(mst(23).ge.1) return
        endif
        if(k(in,1).ge.20000.or.iabs(k(in,2)).lt.500) goto 130
        njet=njet+1
        do 140 j=1,4
  140   dps1(j)=dps1(j)+p(in,j)
        if(k(in,1).ge.10000.or.(mst(6).le.4.and.n.gt.in.and.
     &  k(in+1,1)/10000.eq.1)) goto 130
        nsys=1+in-ip
        mst(1)=ip
        mst(2)=in
        do 150 j=1,3
  150   dbe(j)=dps1(j)/dps1(4)
        dga=1d0/dsqrt(1d0-dbe(1)**2-dbe(2)**2-dbe(3)**2)
        do 180 i=ip,in
        if(mod(k(i,1)/10000,10).ge.6) goto 180
        do 160 j=1,4
  160   dp(j)=p(i,j)
        dbep=-(dbe(1)*dp(1)+dbe(2)*dp(2)+dbe(3)*dp(3))
        dgabep=dga*(dga*dbep/(1d0+dga)+dp(4))
        do 170 j=1,3
  170   p(i,j)=dp(j)-dgabep*dbe(j)
        p(i,4)=dga*(dp(4)+dbep)
  180   continue

c...take or restore spare copy of partons before treatment
        if(n+nsys.ge.mst(30)-5-mst(31)) then
          mst(24)=mst(24)+1
          mst(25)=1
          if(mst(23).ge.1) return
        endif
        ecm=0.
        do 190 i=ip,in
        if(k(i,1).lt.20000.and.iabs(k(i,2)).ge.500) ecm=ecm+p(i,4)
        do 190 j=1,5
  190   p(n+1+i-ip,j)=p(i,j)
        n=n+nsys
        nsav=n
      else
        n=nsav
        do 200 i=ip,in
        if(k(i,1).ge.100000) k(i,1)=k(i,1)-120000
        do 200 j=1,5
  200   p(i,j)=p(nsav+i-in,j)
      endif

      if(mst(6).ge.10.and.ntry.eq.1.and.njet.ge.3) then
c...boost to frame where string tensions balance for montvay scheme
        phi=ulanglC(p(ip,1),p(ip,2))
        call luroboC(0.,-phi,0.,0.,0.)
        the=ulanglC(p(ip,3),p(ip,1))
        call luroboC(-the,0.,0.,0.,0.)
        chi=ulanglC(p(ip+1,1),p(ip+1,2))
        call luroboC(0.,-chi,0.,0.,0.)
        nbal=0
  210   nbal=nbal+1
        do 220 j1=1,3
        te(j1)=0.
        do 220 j2=1,3
  220   td(j1,j2)=0.
        do 240 i=ip,in
        if(k(i,1).ge.20000.or.iabs(k(i,2)).lt.500) goto 240
        pa=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
        ten=min(1.,pa/par(18))
        if(iabs(k(i,2)).eq.500) ten=par(17)*ten
        do 230 j1=1,3
        te(j1)=te(j1)+ten*p(i,j1)/pa
        td(j1,j1)=td(j1,j1)+ten*p(i,4)/pa
        do 230 j2=1,3
  230   td(j1,j2)=td(j1,j2)-ten*p(i,4)*p(i,j1)*p(i,j2)/pa**3
  240   continue
        if(te(1)**2+te(2)**2+te(3)**2.lt.1e-3) goto 260
        if(nbal.ge.mst(13)) goto 100
        do 250 jl=1,2
        do 250 j1=jl+1,3
        te(j1)=te(j1)-(td(j1,jl)/td(jl,jl))*te(jl)
        do 250 j2=jl+1,3
  250   td(j1,j2)=td(j1,j2)-(td(j1,jl)/td(jl,jl))*td(jl,j2)
        te(3)=te(3)/td(3,3)
        te(2)=(te(2)-td(2,3)*te(3))/td(2,2)
        te(1)=(te(1)-td(1,2)*te(2)-td(1,3)*te(3))/td(1,1)
        ter=1.+sqrt(te(1)**2+te(2)**2+te(3)**2)
        call luroboC(0.,0.,-te(1)/ter,-te(2)/ter,-te(3)/ter)
        goto 210
      endif

c...sum and check jet flavours, fragment jets independently
  260 mst(1)=0
      mst(2)=0
      kfsum=0
      do 270 i=ip,in
      kfa=iabs(k(i,2))
      if(k(i,1).ge.20000.or.kfa.lt.500) goto 270
      if(kfa.ge.501) kfsum=kfsum+isign(1,k(i,2)*(510-kfa))
      ifl=mod(kfa,10)
      if(ifl.ne.0.and.ifl.le.3) nfl(ifl)=nfl(ifl)+isign(1,k(i,2))
      ifl=mod(kfa,100)/10
      if(ifl.ne.0.and.ifl.le.3) nfl(ifl)=nfl(ifl)+isign(1,k(i,2))
      ifl=mod(iabs(k(i+1,2)),10)
      if(n.gt.i.and.k(i+1,1)/10000.eq.6.and.ifl.ne.0) then
        kfsum=kfsum+isign(1,ifl)
        nfl(ifl)=nfl(ifl)+isign(1,k(i+1,2))
      endif
      call luonejC(i)
      k(i,1)=k(i,1)+100000
  270 continue
      if(kfsum.ne.0) then
        mst(24)=mst(24)+1
        mst(25)=2
        if(mst(23).ge.1) return
      endif
      if(mod(mst(6),5).ne.0.and.n-nsav.lt.2) goto 100

      if(mst(6).ge.10.and.ntry.eq.1.and.njet.ge.3) then
c...boost back to cm frame for montvay scheme
        do 280 j=1,4
        ps2(j)=0.
        do 280 i=ip,in
  280   if(k(i,1).ge.100000) ps2(j)=ps2(j)+p(i,j)
        mst(1)=ip
        mst(2)=in
        call luroboC(0.,0.,-ps2(1)/ps2(4),-ps2(2)/ps2(4),-ps2(3)/ps2(4))
        phir=ulanglC(p(ip,1),p(ip,2))
        call luroboC(0.,-phir,0.,0.,0.)
        ther=ulanglC(p(ip,3),p(ip,1))
        call luroboC(-ther,0.,0.,0.,0.)
        chir=ulanglC(p(ip+1,1),p(ip+1,2))
        call luroboC(0.,-chir,0.,0.,0.)
        call luroboC(0.,chi,0.,0.,0.)
        call luroboC(the,phi,0.,0.,0.)
        mst(1)=nsav+1
        mst(2)=0
        call luroboC(0.,0.,-ps2(1)/ps2(4),-ps2(2)/ps2(4),-ps2(3)/ps2(4))
        call luroboC(0.,-phir,0.,0.,0.)
        call luroboC(-ther,0.,0.,0.,0.)
        call luroboC(0.,-chir,0.,0.,0.)
        call luroboC(0.,chi,0.,0.,0.)
        call luroboC(the,phi,0.,0.,0.)
        mst(1)=0
      endif

      if(mod(mst(6),10).ne.0) then
c...subtract off produced hadron flavours, finished if zero
        do 290 i=nsav+1,n
        call luiflvC(k(i,2),ifla,iflb,iflc,ksp)
        if(iabs(ifla).le.3) nfl(iabs(ifla))=nfl(iabs(ifla))-isign(1,
     &  ifla)
        if(iabs(iflb).le.3) nfl(iabs(iflb))=nfl(iabs(iflb))-isign(1,
     &  iflb)
  290   if(iflc.ne.0) nfl(iabs(iflc))=nfl(iabs(iflc))-isign(1,iflc)
        nreq=(iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+
     &  nfl(2)+nfl(3)))/2+iabs(nfl(1)+nfl(2)+nfl(3))/3
        if(nreq.eq.0) goto 370

c...take away flavour of low-momentum particles until enough freedom
        nrem=0
  300   irem=0
        p2min=ecm**2
        do 310 i=nsav+1,n
        p2=p(i,1)**2+p(i,2)**2+p(i,3)**2
        if(k(i,1).lt.100000.and.p2.lt.p2min) irem=i
  310   if(k(i,1).lt.100000.and.p2.lt.p2min) p2min=p2
        if(irem.eq.0) goto 100
        k(irem,1)=k(irem,1)+100000
        call luiflvC(k(irem,2),ifla,iflb,iflc,ksp)
        if(iabs(ifla).ge.4) k(irem,1)=k(irem,1)+100000
        if(iabs(ifla).ge.4) goto 300
        nfl(iabs(ifla))=nfl(iabs(ifla))+isign(1,ifla)
        nfl(iabs(iflb))=nfl(iabs(iflb))+isign(1,iflb)
        if(iflc.ne.0) nfl(iabs(iflc))=nfl(iabs(iflc))+isign(1,iflc)
        nrem=nrem+1
        nreq=(iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+
     &  nfl(2)+nfl(3)))/2+iabs(nfl(1)+nfl(2)+nfl(3))/3
        if(nreq.gt.nrem) goto 300
        do 320 i=nsav+1,n
  320   if(k(i,1).ge.200000) k(i,1)=k(i,1)-200000

c...give low-momentum particles new flavours via random combinations
  330   nfet=2
        if(nfl(1)+nfl(2)+nfl(3).ne.0) nfet=3
        if(nreq.lt.nrem) nfet=1
        if(iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3)).eq.0) nfet=0
        do 340 j=1,nfet
        ifet(j)=1+(iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3)))*rluC(0)
        iflf(j)=isign(1,nfl(1))
        if(ifet(j).gt.iabs(nfl(1))) iflf(j)=isign(2,nfl(2))
  340   if(ifet(j).gt.iabs(nfl(1))+iabs(nfl(2))) iflf(j)=isign(3,nfl(3))
        if(nfet.eq.2.and.(ifet(1).eq.ifet(2).or.iflf(1)*iflf(2).gt.0))
     &  goto 330
        if(nfet.eq.3.and.(ifet(1).eq.ifet(2).or.ifet(1).eq.ifet(3).or.
     &  ifet(2).eq.ifet(3).or.iflf(1)*iflf(2).lt.0.or.iflf(1)*iflf(3)
     &  .lt.0.or.iflf(1)*(nfl(1)+nfl(2)+nfl(3)).lt.0)) goto 330
        if(nfet.eq.0) iflf(1)=1+int((2.+par(2))*rluC(0))
        if(nfet.eq.0) iflf(2)=-iflf(1)
        if(nfet.eq.1) iflf(2)=isign(1+int((2.+par(2))*rluC(0)),-iflf(1))
        if(nfet.le.2) iflf(3)=0
        call luifldC(iflf(1),iflf(3),iflf(2),ifldmp,kf)
        if(kf.eq.0) goto 330
        do 350 j=1,max(2,nfet)
  350   nfl(iabs(iflf(j)))=nfl(iabs(iflf(j)))-isign(1,iflf(j))
        npos=min(1+int(rluC(0)*nrem),nrem)
        do 360 i=nsav+1,n
        if(k(i,1).ge.100000) npos=npos-1
        if(k(i,1).lt.100000.or.npos.ne.0) goto 360
        k(i,1)=k(i,1)-100000
        k(i,2)=kf
        p(i,5)=ulmassC(1,k(i,2))
        p(i,4)=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
  360   continue
        nrem=nrem-1
        nreq=(iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+
     &  nfl(2)+nfl(3)))/2+iabs(nfl(1)+nfl(2)+nfl(3))/3
        if(nrem.gt.0) goto 330
      endif

  370 if(mod(mst(6),5).ne.0.and.mod(mst(6),5).ne.4) then
c...compensate for missing momentum in global scheme (3 options)
        do 380 j=1,3
        ps2(j)=0.
        do 380 i=nsav+1,n
  380   ps2(j)=ps2(j)+p(i,j)
        ps2(4)=ps2(1)**2+ps2(2)**2+ps2(3)**2
        pds=0.
        do 390 i=nsav+1,n
        if(mod(mst(6),5).eq.1) pds=pds+p(i,4)
        if(mod(mst(6),5).eq.2) pds=pds+sqrt(p(i,5)**2+(ps2(1)*p(i,1)+
     &  ps2(2)*p(i,2)+ps2(3)*p(i,3))**2/ps2(4))
  390   if(mod(mst(6),5).eq.3) pds=pds+1.
        do 410 i=nsav+1,n
        if(mod(mst(6),5).eq.1) pdm=p(i,4)
        if(mod(mst(6),5).eq.2) pdm=sqrt(p(i,5)**2+(ps2(1)*p(i,1)+
     &  ps2(2)*p(i,2)+ps2(3)*p(i,3))**2/ps2(4))
        if(mod(mst(6),5).eq.3) pdm=1.
        do 400 j=1,3
  400   p(i,j)=p(i,j)-ps2(j)*pdm/pds
  410   p(i,4)=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)

      elseif(mod(mst(6),5).eq.4) then
c...compensate for missing momentum within each jet
        do 420 i=n+1,n+nsys
        k(i,1)=0
        do 420 j=1,5
  420   p(i,j)=0.
        do 440 i=nsav+1,n
        ir1=k(i,1)
        ir2=n+1+ir1-ip
        k(ir2,1)=k(ir2,1)+1
        pls=(p(i,1)*p(ir1,1)+p(i,2)*p(ir1,2)+p(i,3)*p(ir1,3))/
     &  (p(ir1,1)**2+p(ir1,2)**2+p(ir1,3)**2)
        do 430 j=1,3
  430   p(ir2,j)=p(ir2,j)+p(i,j)-pls*p(ir1,j)
        p(ir2,4)=p(ir2,4)+p(i,4)
  440   p(ir2,5)=p(ir2,5)+pls
        hss=0.
        do 450 i=n+1,n+nsys
  450   if(k(i,1).ne.0) hss=hss+p(i,4)/(ecm*(0.8*p(i,5)+0.2))
        do 470 i=nsav+1,n
        ir1=k(i,1)
        ir2=n+1+ir1-ip
        pls=(p(i,1)*p(ir1,1)+p(i,2)*p(ir1,2)+p(i,3)*p(ir1,3))/
     &  (p(ir1,1)**2+p(ir1,2)**2+p(ir1,3)**2)
        do 460 j=1,3
  460   p(i,j)=p(i,j)-p(ir2,j)/k(ir2,1)+(1./(p(ir2,5)*hss)-1.)*pls*
     &  p(ir1,j)
  470   p(i,4)=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
      endif

      if(mod(mst(6),5).ne.0) then
c...scale momenta for energy conservation
        pms=0.
        pes=0.
        pqs=0.
        do 480 i=nsav+1,n
        pms=pms+p(i,5)
        pes=pes+p(i,4)
  480   pqs=pqs+p(i,5)**2/p(i,4)
        if(pms.ge.ecm) goto 100
        neco=0
  490   neco=neco+1
        fac=(ecm-pqs)/(pes-pqs)
        pes=0.
        pqs=0.
        do 510 i=nsav+1,n
        do 500 j=1,3
  500   p(i,j)=fac*p(i,j)
        p(i,4)=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
        pes=pes+p(i,4)
  510   pqs=pqs+p(i,5)**2/p(i,4)
        if(neco.lt.10.and.abs(ecm-pes).gt.2e-6*ecm) goto 490
      endif

c...boost back jets and particles, remove spare copy
      do 540 i=ip,in
      if(k(i,1).ge.100000) k(i,1)=k(i,1)-100000
      if(mod(k(i,1)/10000,10).ge.6) goto 540
      do 520 j=1,4
  520 dp(j)=p(i,j)
      dbep=dbe(1)*dp(1)+dbe(2)*dp(2)+dbe(3)*dp(3)
      dgabep=dga*(dga*dbep/(1d0+dga)+dp(4))
      do 530 j=1,3
  530 p(i,j)=dp(j)+dgabep*dbe(j)
      p(i,4)=dga*(dp(4)+dbep)
  540 continue
      do 580 i=nsav+1,n
      k(i-nsys,1)=k(i,1)
      k(i-nsys,2)=k(i,2)
      do 550 j=1,5
  550 p(i-nsys,j)=p(i,j)
      if(mod(k(i,1)/10000,10).ge.6) goto 580
      do 560 j=1,4
  560 dp(j)=p(i,j)
      dbep=dbe(1)*dp(1)+dbe(2)*dp(2)+dbe(3)*dp(3)
      dgabep=dga*(dga*dbep/(1d0+dga)+dp(4))
      do 570 j=1,3
  570 p(i-nsys,j)=dp(j)+dgabep*dbe(j)
      p(i-nsys,4)=dga*(dp(4)+dbep)
  580 continue
      n=n-nsys

      return
      end

c*********************************************************************

      subroutine luonejC(ip)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      dimension iflf(3),iflo(2),pxo(2),pyo(2),zo(2),wo(2)

c...initial flavour and momentum values, jet along +z axis
      mst(1)=n+1
      iflf(1)=mod(k(ip,2),500)
      iflf(2)=0
      iflf(3)=0
      iflo(2)=0
      wf=p(ip,4)+sqrt(p(ip,1)**2+p(ip,2)**2+p(ip,3)**2)
      zj=0.
      zi=0.

  100 if(iflf(1).eq.0.and.mst(5).le.2) then
c...initial values and first rank hadron in gluon jet (lund model)
        ns=2
        i=n+1
        k(i,1)=ip
        call luifldC(int(1.+(2.+par(2))*rluC(0)),0,0,iflo(1),kdump)
        call luifldC(-int(1.+(2.+par(2))*rluC(0)),0,0,iflo(2),kdump)
        if(iabs(iflo(2)).lt.100) iflo(2)=-mod(iflo(1),100)
        if(iabs(iflo(2)).gt.100) iflo(1)=-mod(iflo(2),100)
        call luifldC(iflo(1),0,0,ifl1,k(i,2))
        iflo(1)=-ifl1
        p(i,5)=ulmassC(1,k(i,2))
        call luptdiC(iflo(1),pxo(1),pyo(1))
        call luptdiC(iflo(2),pxo(2),pyo(2))
        pr=p(i,5)**2+(pxo(1)+pxo(2))**2+(pyo(1)+pyo(2))**2
        prdiv=rluC(0)*pr
        do 110 jt=1,2
        call luzdisC(iflo(jt),0,0.6*((2-jt)*pr+(2*jt-3)*prdiv),zo(jt))
  110   wo(jt)=0.5*(1.-zo(jt))*wf

c...four-momentum for hadron, optionally skip it if moving backwards
        p(i,1)=-(pxo(1)+pxo(2))
        p(i,2)=-(pyo(1)+pyo(2))
        p(i,3)=0.25*(zo(1)+zo(2))*wf-pr/((zo(1)+zo(2))*wf)
        p(i,4)=0.25*(zo(1)+zo(2))*wf+pr/((zo(1)+zo(2))*wf)
        if(mst(5).ge.2.and.mst(6).ge.0.and.p(i,3).lt.0.) i=i-1
        n=i

      elseif(iflf(1).eq.0.and.(mst(5).eq.3.or.mst(5).eq.4)) then
c...gluon treated like random quark (or antiquark) jet
        ns=1
        if(mst(5).eq.4) mst(32)=1
        iflo(1)=int(1.+(2.+par(2))*rluC(0))*(-1)**int(rluC(0)+0.5)
        call luptdiC(93,pxo(1),pyo(1))
        wo(1)=wf

      elseif(iflf(1).eq.0.and.mst(5).ge.5) then
c...gluon treated like quark-antiquark jet pair, sharing energy
c...according to altarelli-parisi splitting function
        ns=2
        if(mst(5).eq.6) mst(32)=1
        iflo(1)=int(1.+(2.+par(2))*rluC(0))*(-1)**int(rluC(0)+0.5)
        iflo(2)=-iflo(1)
        call luptdiC(93,pxo(1),pyo(1))
        pxo(2)=-pxo(1)
        pyo(2)=-pyo(1)
        wo(1)=wf*rluC(0)**(1./3.)
        wo(2)=wf-wo(1)

      else
c...initial values for quark, diquark or hadron jet
        ns=1
        iflo(1)=iflf(1)
        call luptdiC(93,pxo(1),pyo(1))
        wo(1)=wf

        if(mod(mst(10),2).eq.1.and.iabs(iflf(1)).gt.10) then
c...order and position of quarks in diquark (l and j quarks)
          ifla=iflf(1)/10
          iflb=iflf(1)-10*ifla
          iflf(2)=ifla+int(rluC(0)+0.5)*(iflb-ifla)
          if(n.gt.ip.and.k(ip+1,1)/10000.eq.6.and.iabs(k(ip+1,2)).ge.
     &    610) iflf(2)=mod(k(ip+1,2)/10,10)
          iflo(1)=ifla+iflb-iflf(2)
          call luzdisC(0,1,0.,zj)
          if(n.gt.ip.and.k(ip+1,1)/10000.eq.6) p(ip+1,1)=zj
        endif

c...flavour and position of extra quark in hadron jet (i quark)
        if(n.gt.ip.and.k(ip+1,1)/10000.eq.6) iflf(3)=mod(k(ip+1,2),10)
        if(iflf(3).ne.0) then
          call luzdisC(0,2+(90+iabs(iflf(1)))/100,0.,zi)
          if(iabs(iflf(1)).gt.10.and.mod(mst(10),2).eq.1) zi=zi*zj
          p(ip+1,3)=zi
        endif
      endif

c...initial values for rank, flavour, pt (including relative pt
c...in diquark) and longitudinal fragmentation variables
      do 140 jt=1,ns
  120 i=n
      lrk=0
      ifl1=iflo(jt)
      iflj=iflf(2)
      ifli=iflf(3)
      if(iflj.eq.0) then
        px1=pxo(jt)
        py1=pyo(jt)
        pxj=0.
        pyj=0.
      else
        call luptdiC(94,pxr,pyr)
        px1=0.5*pxo(jt)+pxr
        py1=0.5*pyo(jt)+pyr
        pxj=0.5*pxo(jt)-pxr
        pyj=0.5*pyo(jt)-pyr
        if(n.gt.ip.and.k(ip+1,1)/10000.eq.6) p(ip+1,2)=0.
      endif
      if(iflf(3).ne.0) p(ip+1,4)=0.
      w=wo(jt)

c...new hadron: generate pt
  130 i=i+1
      if(i.ge.mst(30)-5-mst(31)) then
        mst(24)=mst(24)+1
        mst(25)=1
        if(mst(23).ge.1) return
      endif
      lrk=lrk+1
      k(i,1)=ip
      call luptdiC(ifl1,px2,py2)
      mqj=0
      mqi=0

c...check if j or i quark to be included, generate flavour and hadron
      if(iflj.ne.0.or.ifli.ne.0) then
        prji=par(37)**2+(px1+px2)**2+(py1+py2)**2
        call luzdisC(ifl1,iflj+ifli,prji,z)
        if(iflj.ne.0.and.(1.-z)*w.le.zj*wf) mqj=1
        if(mqj.eq.1.and.iabs(ifl1).gt.10) goto 120
        if(mqj.eq.1.and.lrk.eq.1) ifl1=iflf(1)
        if(ifli.ne.0.and.(1.-z)*w.le.zi*wf) mqi=1
        if(mqi.eq.1.and.iabs(ifl1).gt.100) goto 120
      endif
      call luifldC(ifl1,mqj*iflj,mqi*ifli,ifl2,k(i,2))
      if(k(i,2).eq.0) goto 120
      p(i,5)=ulmassC(1,k(i,2))
      pr=p(i,5)**2+(px1+px2+mqj*pxj)**2+(py1+py2+mqj*pyj)**2

c...four-momentum for hadron, optionally skip it if moving backwards;
c...position of j and i quark
      if(iflj.eq.0.and.ifli.eq.0) then
        call luzdisC(ifl1,0,pr,z)
      elseif(mst(4).eq.1.or.mst(4).eq.3) then
        gamji=(1.+par(35))/par(36)
        zbc=(pr-prji-z*gamji+prji/z)/(2.*gamji)
        z=sqrt(zbc**2+pr/gamji)-zbc
      endif
      p(i,1)=px1+px2+mqj*pxj
      p(i,2)=py1+py2+mqj*pyj
      p(i,3)=0.5*(z*w-pr/(z*w))
      p(i,4)=0.5*(z*w+pr/(z*w))
      if(mod(mst(6),10).gt.0.and.lrk.eq.1.and.max(mod(iabs(iflf(1)),
     &10),iabs(iflf(1))/10).ge.4.and.p(i,3).le.0.001) then
        if(w.ge.p(i,5)+0.5*par(22)) goto 120
        p(i,3)=0.0001
        p(i,4)=sqrt(pr)
        z=p(i,4)/w
      endif
      if(mst(5).ge.2.and.mst(6).ge.0.and.p(i,3).lt.0.) i=i-1
      if(i.eq.n+lrk.and.mqj*n.gt.ip.and.k(ip+1,1)/10000.eq.6)
     &p(ip+1,2)=i
      if(i.eq.n+lrk.and.mqi.eq.1) p(ip+1,4)=i

c...remaining flavour and momentum, go back and generate new hadron
      ifl1=-ifl2
      if(mqi.eq.1) call luifldC((-1)**int(rluC(0)+0.5),0,0,ifl1,kdump)
      if(mqj.eq.1) iflj=0
      if(mqi.eq.1) ifli=0
      px1=-px2
      py1=-py2
      w=(1.-z)*w
      if(mst(10).le.1.and.w.gt.par(21)) goto 130
      if(mst(10).ge.2.and.(iflj.ne.0.or.ifli.ne.0)) goto 130
  140 n=i

c...rotate jet to right direction
      if(mod(mst(6),5).eq.4.and.mst(1).eq.n+1) wf=wf+0.1*par(22)
      if(mod(mst(6),5).eq.4.and.mst(1).eq.n+1) goto 100
      the=ulanglC(p(ip,3),sqrt(p(ip,1)**2+p(ip,2)**2))
      phi=ulanglC(p(ip,1),p(ip,2))
      call luroboC(the,phi,0.,0.,0.)
      mst(1)=0
      mst(32)=0
      k(ip,1)=k(ip,1)+20000

      return
      end

c*********************************************************************

      subroutine lusysjC(ip)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
c        &&&&&&&&&&
      real*8 prev2
c        &&&&&&&&&&
      dimension ps(5),ifl(3),px(4),py(4),gam(3),pr(2),in(9),hm(4),hg(4),
     &lrk(2),ie(2),iflf(3),iflj(2),ifli(2),pxj(2),pyj(2),zj(2),zi(2),
     &zpos(2),pmq(3)
      double precision dp(5,5),dfour,hkc,hks,hk1,hk2,hc12,hcx1,hcx2,
     &hcxx,hcy1,hcy2,hcyx,hcyy
c...function: four-product of two vectors
      four(i,j)=p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3)
      dfour(i,j)=dp(i,4)*dp(j,4)-dp(i,1)*dp(j,1)-dp(i,2)*dp(j,2)-
     &dp(i,3)*dp(j,3)

c...begin kinematics definition: identify jets in system
      ntry=0
      np=0
      do 100 j=1,5
  100 ps(j)=0.
      i=ip-1
  110 i=i+1
      if(i.gt.min(n,mst(30)-np-5-mst(31))) then
        mst(24)=mst(24)+1
        mst(25)=2
        if(i.le.n) mst(25)=1
        if(mst(23).ge.1) return
      endif
      if(k(i,1).ge.20000.or.iabs(k(i,2)).lt.500) goto 110
      np=np+1
      k(n+np,1)=i
      k(n+np,2)=k(i,2)
      do 120 j=1,5
      p(n+np,j)=p(i,j)
  120 ps(j)=ps(j)+p(i,j)
      if(p(n+np,4)**2.lt.p(n+np,1)**2+p(n+np,2)**2+p(n+np,3)**2) then
        p(n+np,4)=sqrt(p(n+np,1)**2+p(n+np,2)**2+p(n+np,3)**2+
     &  p(n+np,5)**2)
        ps(4)=ps(4)+max(0.,p(n+np,4)-p(i,4))
      endif
      if(k(i,1).ge.10000) goto 110

c...boost to cm frame for rapidly moving system
      mbst=0
      if(ps(1)**2+ps(2)**2+ps(3)**2.gt.0.5*ps(4)**2) then
        mbst=1
        pebst=max(ps(4),1.0001*sqrt(ps(1)**2+ps(2)**2+ps(3)**2))
        mst(1)=n+1
        mst(2)=n+np
        call luroboC(0.,0.,-ps(1)/pebst,-ps(2)/pebst,-ps(3)/pebst)
      endif

c...search for very nearby partons that may be recombined
      nr=np
  130 if(nr.le.2) goto 180
      drmin=2.*par(59)
      do 140 i=n+1,n+nr
      if(i.eq.n+nr.and.iabs(k(n+1,2)).ne.500) goto 140
      i1=i+1-nr*(i/(n+nr))
      pap=sqrt((p(i,1)**2+p(i,2)**2+p(i,3)**2)*(p(i1,1)**2+
     &p(i1,2)**2+p(i1,3)**2))
      pvp=p(i,1)*p(i1,1)+p(i,2)*p(i1,2)+p(i,3)*p(i1,3)
      dr=4.*(pap-pvp)**2/(par(60)**2*pap+2.*(pap-pvp))
      if(dr.lt.drmin) then
        ir=i
        drmin=dr
      endif
  140 continue

c...recombine very nearby partons to avoid machine precision problems
      if(drmin.lt.par(59).and.ir.eq.n+nr) then
        do 150 j=1,4
  150   p(n+1,j)=p(n+1,j)+p(n+nr,j)
        p(n+1,5)=sqrt(max(0.,p(n+1,4)**2-p(n+1,1)**2-p(n+1,2)**2-
     &  p(n+1,3)**2))
        nr=nr-1
        goto 130
      elseif(drmin.lt.par(59)) then
        do 160 j=1,4
  160   p(ir,j)=p(ir,j)+p(ir+1,j)
        p(ir,5)=sqrt(max(0.,p(ir,4)**2-p(ir,1)**2-p(ir,2)**2-
     &  p(ir,3)**2))
        do 170 i=ir+1,n+nr-1
        k(i,2)=k(i+1,2)
        do 170 j=1,5
  170   p(i,j)=p(i+1,j)
        if(ir.eq.n+nr-1) k(ir,2)=k(n+nr,2)
        nr=nr-1
        goto 130
      endif
  180 if(n+5*nr+11.ge.mst(30)-5-mst(31)) then
        mst(24)=mst(24)+1
        mst(25)=1
        if(mst(23).ge.1) return
      endif

c...open versus closed strings, choose breakup region for latter
      if(iabs(k(n+1,2)).ne.500) then
        ns=nr-1
        nb=1
      else
        ns=nr+1
        w2sum=0.
        do 190 is=1,nr
        p(n+nr+is,1)=0.5*four(n+is,n+is+1-nr*(is/nr))
  190   w2sum=w2sum+p(n+nr+is,1)
        w2ran=rluC(0)*w2sum
        nb=0
  200   nb=nb+1
        w2sum=w2sum-p(n+nr+nb,1)
        if(w2sum.gt.w2ran.and.nb.lt.nr) goto 200
      endif

c...find longitudinal string directions (i.e. lightlike four-vectors)
      do 220 is=1,ns
      is1=n+is+nb-1-nr*((is+nb-2)/nr)
      is2=n+is+nb-nr*((is+nb-1)/nr)
      do 210 j=1,5
      dp(1,j)=p(is1,j)
      if(iabs(k(is1,2)).eq.500) dp(1,j)=0.5*dp(1,j)
      dp(2,j)=p(is2,j)
  210 if(iabs(k(is2,2)).eq.500) dp(2,j)=0.5*dp(2,j)
      dp(3,5)=dfour(1,1)
      dp(4,5)=dfour(2,2)
      hkc=dfour(1,2)
      if(dp(3,5)+2.*hkc+dp(4,5).le.0.) then
        dp(3,5)=dp(1,5)**2
        dp(4,5)=dp(2,5)**2
        dp(1,4)=dsqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2+dp(1,5)**2)
        dp(2,4)=dsqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2+dp(2,5)**2)
        hkc=dfour(1,2)
      endif
      hks=dsqrt(hkc**2-dp(3,5)*dp(4,5))
      hk1=0.5*((dp(4,5)+hkc)/hks-1.)
      hk2=0.5*((dp(3,5)+hkc)/hks-1.)
      in1=n+nr+4*is-3
      p(in1,5)=sqrt(dp(3,5)+2.*hkc+dp(4,5))
      do 220 j=1,4
      p(in1,j)=(1.+hk1)*dp(1,j)-hk2*dp(2,j)
  220 p(in1+1,j)=(1.+hk2)*dp(2,j)-hk1*dp(1,j)
      nrs=nr+4*ns+7

c...begin initialization: sum up energy, set starting positions
  230 ntry=ntry+1
      if(ntry.gt.200) then
        mst(24)=mst(24)+1
        mst(25)=3
        if(mst(23).ge.1) return
      endif
      i=n+nrs
      do 240 j=1,4
      p(i,j)=0.
      do 240 is=1,nr
  240 p(i,j)=p(i,j)+p(n+is,j)
      do 250 jt=1,2
      lrk(jt)=0
      ie(jt)=k(n+1+(jt/2)*(np-1),1)
      iflj(jt)=0
      ifli(jt)=0
      pxj(jt)=0.
      pyj(jt)=0.
      zj(jt)=0.
      zi(jt)=0.
      zpos(jt)=1.
      in(3*jt+1)=n+nr+1+4*(jt/2)*(ns-1)
      in(3*jt+2)=in(3*jt+1)+1
      in(3*jt+3)=n+nr+4*ns+2*jt-1
      do 250 in1=n+nr+2+jt,n+nr+4*ns-2+jt,4
      p(in1,1)=2-jt
      p(in1,2)=jt-1
  250 p(in1,3)=1.
      iflstr=0

      if(ns.eq.nr-1) then
c...initialize flavour and pt variables for open string
        px(1)=0.
        py(1)=0.
        if(ns.eq.1) call luptdiC(93,px(1),py(1))
        px(2)=-px(1)
        py(2)=-py(1)
        kfsum=0
        do 260 jt=1,2
        iflf(jt)=mod(k(ie(jt),2),500)
        kfsum=kfsum+isign(1,iflf(jt)*(10-iabs(iflf(jt))))
        ifl(jt)=iflf(jt)
        gam(jt)=0.

c...order and position of quarks in diquark, relative pt in diquark
        if(mod(mst(10),2).eq.1.and.iabs(iflf(jt)).ge.10) then
          ifla=iflf(jt)/10
          iflb=iflf(jt)-10*ifla
          iflj(jt)=ifla+int(rluC(0)+0.5)*(iflb-ifla)
          if(n.gt.ie(jt).and.k(ie(jt)+1,1)/10000.eq.6.and.iabs(k(ie(jt)+
     &    1,2)).ge.610) iflj(jt)=mod(k(ie(jt)+1,2)/10,10)
          ifl(jt)=ifla+iflb-iflj(jt)
          call luzdisC(0,1,0.,zj(jt))
          call luptdiC(94,pxr,pyr)
          px(jt)=0.5*px(jt)+pxr
          py(jt)=0.5*py(jt)+pyr
          pxj(jt)=px(jt)-2.*pxr
          pyj(jt)=py(jt)-2.*pyr
          if(n.gt.ie(jt).and.k(ie(jt)+1,1)/10000.eq.6) p(ie(jt)+1,1)=
     &    zj(jt)
        endif
        pmq(jt)=ulmassC(2,ifl(jt))

c...flavour and position of extra quarks in hadron jets (i quarks)
        if(n.gt.ie(jt).and.k(ie(jt)+1,1)/10000.eq.6) ifli(jt)=
     &  mod(k(ie(jt)+1,2),10)
        if(ifli(jt).ne.0) then
          kfsum=kfsum+isign(1,ifli(jt))
          call luzdisC(0,2+(90+iabs(iflf(jt)))/100,0.,zi(jt))
          if(iabs(iflf(jt)).gt.10.and.mod(mst(10),2).eq.1) zi(jt)=
     &    zi(jt)*zj(jt)
          p(ie(jt)+1,3)=zi(jt)
        endif
        call luifldC(int(1.+(2.+par(2))*rluC(0))*(-1)**int(rluC(0)+0.5),
     &  0,0,iflf(3),kdump)
  260   iflf(3)=mod(iflf(3),100)
        if(kfsum.ne.0) then
          mst(24)=mst(24)+1
          mst(25)=2
          if(mst(23).ge.1) return
        endif

      else
c...closed string: random initial breakup flavour, pt and vertex
        ifl(3)=int(1.+(2.+par(2))*rluC(0))*(-1)**int(rluC(0)+0.5)
        call luifldC(ifl(3),0,0,ifl(1),kdump)
        call luifldC(-ifl(3),0,0,ifl(2),kdump)
        if(iabs(ifl(2)).lt.100) ifl(2)=-mod(ifl(1),100)
        if(iabs(ifl(2)).gt.100) ifl(1)=-mod(ifl(2),100)
        iflstr=(iabs(ifl(1))+90)/100
        call luptdiC(ifl(1),px(1),py(1))
        px(2)=-px(1)
        py(2)=-py(1)
        pr3=min(25.,0.1*p(n+nr+1,5)**2)
  270   call luzdisC(ifl(1),0,pr3,z)
        zr=pr3/(z*p(n+nr+1,5)**2)
        if(zr.ge.1.) goto 270
        do 280 jt=1,2
        pmq(jt)=ulmassC(2,ifl(jt))
        gam(jt)=pr3*(1.-z)/z
        in1=n+nr+3+4*(jt/2)*(ns-1)
        p(in1,jt)=1.-z
        p(in1,3-jt)=jt-1
        p(in1,3)=(2-jt)*(1.-z)+(jt-1)*z
        p(in1+1,jt)=zr
        p(in1+1,3-jt)=2-jt
  280   p(in1+1,3)=(2-jt)*(1.-zr)+(jt-1)*zr
      endif

c...find initial transverse directions (i.e. spacelike four-vectors)
      do 320 jt=1,2
      if(jt.eq.1.or.ns.eq.nr-1) then
        in1=in(3*jt+1)
        in3=in(3*jt+3)
        do 290 j=1,4
        dp(1,j)=p(in1,j)
        dp(2,j)=p(in1+1,j)
        dp(3,j)=0.
  290   dp(4,j)=0.
        dp(1,4)=dsqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
        dp(2,4)=dsqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
        dp(5,1)=dp(1,1)/dp(1,4)-dp(2,1)/dp(2,4)
        dp(5,2)=dp(1,2)/dp(1,4)-dp(2,2)/dp(2,4)
        dp(5,3)=dp(1,3)/dp(1,4)-dp(2,3)/dp(2,4)
        if(dp(5,1)**2.le.dp(5,2)**2+dp(5,3)**2) dp(3,1)=1.
        if(dp(5,1)**2.gt.dp(5,2)**2+dp(5,3)**2) dp(3,3)=1.
        if(dp(5,2)**2.le.dp(5,1)**2+dp(5,3)**2) dp(4,2)=1.
        if(dp(5,2)**2.gt.dp(5,1)**2+dp(5,3)**2) dp(4,3)=1.
        hc12=dfour(1,2)
        hcx1=dfour(3,1)/hc12
        hcx2=dfour(3,2)/hc12
        hcxx=1./dsqrt(1.+2.*hcx1*hcx2*hc12)
        hcy1=dfour(4,1)/hc12
        hcy2=dfour(4,2)/hc12
        hcyx=hcxx*(hcx1*hcy2+hcx2*hcy1)*hc12
        hcyy=1./dsqrt(1.+2.*hcy1*hcy2*hc12-hcyx**2)
        do 300 j=1,4
        dp(3,j)=hcxx*(dp(3,j)-hcx2*dp(1,j)-hcx1*dp(2,j))
        p(in3,j)=dp(3,j)
  300   p(in3+1,j)=hcyy*(dp(4,j)-hcy2*dp(1,j)-hcy1*dp(2,j)-
     &  hcyx*dp(3,j))
      else
        do 310 j=1,4
        p(in3+2,j)=p(in3,j)
  310   p(in3+3,j)=p(in3+1,j)
      endif
  320 continue

c...produce new particle: side, pt
  330 i=i+1
      if(i.ge.mst(30)-5-mst(31)) then
        mst(24)=mst(24)+1
        mst(25)=1
        if(mst(23).ge.1) return
      endif
      jt=1.5+rluC(0)
      if(iabs(ifl(3-jt)).gt.100) jt=3-jt
      lrk(jt)=lrk(jt)+1
      jr=3-jt
      js=3-2*jt
      k(i,1)=ie(jt)
      call luptdiC(ifl(jt),px(3),py(3))
      px(4)=px(3)
      py(4)=py(3)
      mqj=0
      mqi=0

c...check if j or i quark to be included, generate flavour and hadron
      if(iflj(jt).ne.0.or.ifli(jt).ne.0) then
        prji=par(37)**2+(px(jt)+px(3))**2+(py(jt)+py(3))**2
        call luzdisC(ifl(jt),iflj(jt)+ifli(jt),prji,z)
        if(iflj(jt).ne.0.and.(1.-z)*zpos(jt).le.zj(jt)) mqj=1
        if(mqj.eq.1.and.iabs(ifl(jt)).gt.10) goto 230
        if(mqj.eq.1.and.lrk(jt).eq.1) ifl(jt)=iflf(jt)
        if(mqj.eq.1.and.lrk(jt).eq.1) pmq(jt)=ulmassC(2,ifl(jt))
        if(ifli(jt).ne.0.and.(1.-z)*zpos(jt).le.zi(jt)) mqi=1
        if(mqi.eq.1.and.iabs(ifl(jt)).gt.100) goto 230
        if(mqj.eq.1) px(jt)=px(jt)+pxj(jt)
        if(mqj.eq.1) py(jt)=py(jt)+pyj(jt)
        if(mqj*n.gt.ie(jt).and.k(ie(jt)+1,1)/10000.eq.6) p(ie(jt)+1,2)=
     &  i-nrs
        if(mqi.eq.1) p(ie(jt)+1,4)=i-nrs
      endif
      call luifldC(ifl(jt),mqj*iflj(jt),mqi*ifli(jt),ifl(3),k(i,2))
      if(k(i,2).eq.0) goto 230
      pmq(3)=ulmassC(2,ifl(3))

c...final hadrons for small invariant mass, particle mass
      wmin=par(22+mst(4))+pmq(1)+pmq(2)+par(26)*pmq(3)
      if(iflj(jt).ne.0.and.mqj.eq.0) wmin=wmin+ulmassC(2,iflj(jt))
      if(iflj(jr).ne.0) wmin=wmin+ulmassC(2,iflj(jr))
      if(iabs(ifl(jt)).gt.100) wmin=wmin+par(26)*(ulmassC(2,mst(33))-
     &pmq(3))
      wrem2=four(n+nrs,n+nrs)
      if(wrem2.lt.0.10) goto 230
      if(wrem2.lt.max(wmin*(1.+(2.*rluC(0)-1.)*par(27)),
     &par(22)+pmq(1)+pmq(2))**2) goto 460
      p(i,5)=ulmassC(1,k(i,2))
      pr(jt)=p(i,5)**2+(px(jt)+px(3))**2+(py(jt)+py(3))**2

c...choose z (gives gamma), shift z for heavy flavours, i and j quarks
      if(iflj(jt).eq.0.and.ifli(jt).eq.0) then
        call luzdisC(ifl(jt),0,pr(jt),z)
        if(max(mod(iabs(ifl(1)),10),mod(iabs(ifl(2)),10)).ge.4) then
          pr(jr)=(pmq(jr)+pmq(3))**2+(px(jr)-px(3))**2+(py(jr)-py(3))**2
          pw12=sqrt(max(0.,(wrem2-pr(1)-pr(2))**2-4.*pr(1)*pr(2)))
          z=(wrem2+pr(jt)-pr(jr)+pw12*(2.*z-1.))/(2.*wrem2)
          pr(jr)=(pmq(jr)+par(22+mst(4)))**2+(px(jr)-px(3))**2+
     &    (py(jr)-py(3))**2
          if((1.-z)*(wrem2-pr(jt)/z).lt.pr(jr)) goto 460
        endif
      elseif(mst(4).eq.1.or.mst(4).eq.3) then
        gamji=(1.+par(35))/par(36)
        zbc=(pr(jt)-prji-z*gamji+prji/z)/(2.*gamji)
        z=sqrt(zbc**2+pr(jt)/gamji)-zbc
      endif
      gam(3)=(1.-z)*(gam(jt)+pr(jt)/z)
      do 340 j=1,3
  340 in(j)=in(3*jt+j)

c...stepping within or from 'low' region easy
      if(in(1)+1.eq.in(2).and.z*p(in(1)+2,3)*p(in(2)+2,3)*
     &p(in(1),5)**2.ge.pr(jt)) then
        p(in(jt)+2,4)=z*p(in(jt)+2,3)
        p(in(jr)+2,4)=pr(jt)/(p(in(jt)+2,4)*p(in(1),5)**2)
        do 350 j=1,4
  350   p(i,j)=(px(jt)+px(3))*p(in(3),j)+(py(jt)+py(3))*p(in(3)+1,j)
        goto 420
      elseif(in(1)+1.eq.in(2)) then
        p(in(jr)+2,4)=p(in(jr)+2,3)
        p(in(jr)+2,jt)=1.
        in(jr)=in(jr)+4*js
        if(js*in(jr).gt.js*in(4*jr)) goto 230
        if(four(in(1),in(2)).le.1e-2) then
          p(in(jt)+2,4)=p(in(jt)+2,3)
          p(in(jt)+2,jt)=0.
          in(jt)=in(jt)+4*js
        endif
      endif

c...find new transverse directions (i.e. spacelike four-vectors)
  360 if(js*in(1).gt.js*in(3*jr+1).or.js*in(2).gt.js*in(3*jr+2).or.
     &in(1).gt.in(2)) goto 230
      if(in(1).ne.in(3*jt+1).or.in(2).ne.in(3*jt+2)) then
        do 370 j=1,4
        dp(1,j)=p(in(1),j)
        dp(2,j)=p(in(2),j)
        dp(3,j)=0.
  370   dp(4,j)=0.
        dp(1,4)=dsqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
        dp(2,4)=dsqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
        hc12=dfour(1,2)
        if(hc12.le.1e-2) then
          p(in(jt)+2,4)=p(in(jt)+2,3)
          p(in(jt)+2,jt)=0.
          in(jt)=in(jt)+4*js
          goto 360
        endif
        in(3)=n+nr+4*ns+5
        dp(5,1)=dp(1,1)/dp(1,4)-dp(2,1)/dp(2,4)
        dp(5,2)=dp(1,2)/dp(1,4)-dp(2,2)/dp(2,4)
        dp(5,3)=dp(1,3)/dp(1,4)-dp(2,3)/dp(2,4)
        if(dp(5,1)**2.le.dp(5,2)**2+dp(5,3)**2) dp(3,1)=1.
        if(dp(5,1)**2.gt.dp(5,2)**2+dp(5,3)**2) dp(3,3)=1.
        if(dp(5,2)**2.le.dp(5,1)**2+dp(5,3)**2) dp(4,2)=1.
        if(dp(5,2)**2.gt.dp(5,1)**2+dp(5,3)**2) dp(4,3)=1.
        hcx1=dfour(3,1)/hc12
        hcx2=dfour(3,2)/hc12
        hcxx=1./dsqrt(1.+2.*hcx1*hcx2*hc12)
        hcy1=dfour(4,1)/hc12
        hcy2=dfour(4,2)/hc12
        hcyx=hcxx*(hcx1*hcy2+hcx2*hcy1)*hc12
        hcyy=1./dsqrt(1.+2.*hcy1*hcy2*hc12-hcyx**2)
        do 380 j=1,4
        dp(3,j)=hcxx*(dp(3,j)-hcx2*dp(1,j)-hcx1*dp(2,j))
        p(in(3),j)=dp(3,j)
  380   p(in(3)+1,j)=hcyy*(dp(4,j)-hcy2*dp(1,j)-hcy1*dp(2,j)-
     &  hcyx*dp(3,j))
c...express pt with respect to new axes if sensible
        px(3)=-(px(4)*four(in(3*jt+3),in(3))+py(4)*
     &  four(in(3*jt+3)+1,in(3)))
        py(3)=-(px(4)*four(in(3*jt+3),in(3)+1)+py(4)*
     &  four(in(3*jt+3)+1,in(3)+1))
        if(abs(px(3)**2+py(3)**2-px(4)**2-py(4)**2).gt.0.01) then
          px(3)=px(4)
          py(3)=py(4)
        endif
      endif

c...sum up known four-momentum, gives coefficients for m2 expression
      do 400 j=1,4
      hg(j)=0.
      p(i,j)=px(jt)*p(in(3*jt+3),j)+py(jt)*p(in(3*jt+3)+1,j)+
     &px(3)*p(in(3),j)+py(3)*p(in(3)+1,j)
      do 390 in1=in(3*jt+1),in(1)-4*js,4*js
  390 p(i,j)=p(i,j)+p(in1+2,3)*p(in1,j)
      do 400 in2=in(3*jt+2),in(2)-4*js,4*js
  400 p(i,j)=p(i,j)+p(in2+2,3)*p(in2,j)
      hm(1)=four(i,i)
      hm(2)=2.*four(i,in(1))
      hm(3)=2.*four(i,in(2))
      hm(4)=2.*four(in(1),in(2))

c...find coefficients for gamma expression
      do 410 in2=in(1)+1,in(2),4
      do 410 in1=in(1),in2-1,4
      hc=2.*four(in1,in2)
      hg(1)=hg(1)+p(in1+2,jt)*p(in2+2,jt)*hc
      if(in1.eq.in(1)) hg(2)=hg(2)-js*p(in2+2,jt)*hc
      if(in2.eq.in(2)) hg(3)=hg(3)+js*p(in1+2,jt)*hc
  410 if(in1.eq.in(1).and.in2.eq.in(2)) hg(4)=hg(4)-hc

c...solve mass-square, gamma equation system for energies taken
      hs1=hm(jr+1)*hg(4)-hm(4)*hg(jr+1)
      if(abs(hs1).lt.1e-4) goto 230
      hs2=hm(4)*(gam(3)-hg(1))-hm(jt+1)*hg(jr+1)-hg(4)*
     &(p(i,5)**2-hm(1))+hg(jt+1)*hm(jr+1)
      hs3=hm(jt+1)*(gam(3)-hg(1))-hg(jt+1)*(p(i,5)**2-hm(1))
      p(in(jr)+2,4)=0.5*(sqrt(max(0.,hs2**2-4.*hs1*hs3))/abs(hs1)-
     &hs2/hs1)
      if(hm(jt+1)+hm(4)*p(in(jr)+2,4).le.0.) goto 230
      p(in(jt)+2,4)=(p(i,5)**2-hm(1)-hm(jr+1)*p(in(jr)+2,4))/
     &(hm(jt+1)+hm(4)*p(in(jr)+2,4))

c...step to new region if necessary
      if(p(in(jr)+2,4).gt.p(in(jr)+2,3)) then
        p(in(jr)+2,4)=p(in(jr)+2,3)
        p(in(jr)+2,jt)=1.
        in(jr)=in(jr)+4*js
        if(js*in(jr).gt.js*in(4*jr)) goto 230
        if(four(in(1),in(2)).le.1e-2) then
          p(in(jt)+2,4)=p(in(jt)+2,3)
          p(in(jt)+2,jt)=0.
          in(jt)=in(jt)+4*js
        endif
        goto 360
      elseif(p(in(jt)+2,4).gt.p(in(jt)+2,3)) then
        p(in(jt)+2,4)=p(in(jt)+2,3)
        p(in(jt)+2,jt)=0.
        in(jt)=in(jt)+4*js
        goto 360
      endif

c...four-momentum of particle, remaining quantities, loop back
  420 do 430 j=1,4
      p(i,j)=p(i,j)+p(in(1)+2,4)*p(in(1),j)+p(in(2)+2,4)*p(in(2),j)
  430 p(n+nrs,j)=p(n+nrs,j)-p(i,j)
      if(p(i,4).le.0.) goto 230
      ifl(jt)=-ifl(3)
      pmq(jt)=pmq(3)
      if(mqi.eq.1) ifl(jt)=iflf(3)*(-1)**jt
      if(mqi.eq.1) pmq(jt)=ulmassC(2,ifl(jt))
      if(mqj.eq.1) iflj(jt)=0
      if(mqi.eq.1) ifli(jt)=0
      zpos(jt)=(1.-z)*zpos(jt)
      px(jt)=-px(3)
      py(jt)=-py(3)
      gam(jt)=gam(3)
      if(in(3).ne.in(3*jt+3)) then
        do 440 j=1,4
        p(in(3*jt+3),j)=p(in(3),j)
  440   p(in(3*jt+3)+1,j)=p(in(3)+1,j)
      endif
      do 450 jq=1,2
      in(3*jt+jq)=in(jq)
      p(in(jq)+2,3)=p(in(jq)+2,3)-p(in(jq)+2,4)
  450 p(in(jq)+2,jt)=p(in(jq)+2,jt)-js*(3-2*jq)*p(in(jq)+2,4)
      goto 330

c...final two hadrons: side information for last one
  460 if(max(iabs(ifl(jr)),iabs(ifl(3))).gt.100) goto 230
      do 470 jf=jt,jr,js
      if(jf.eq.jr) i=i+1
      if(jf.eq.jr) lrk(jf)=lrk(jf)+1
      if(jf.eq.jr) k(i,1)=ie(jf)

c...accept generated flavour when no j or i quarks remaining
      if(jf.eq.jt.and.(1-mqj)*iflj(jf).eq.0.and.(1-mqi)*ifli(jf)
     & .eq.0.and.ifli(3-jf).eq.0) then
      elseif(jf.eq.jr.and.iflj(jf).eq.0.and.ifli(jf).eq.0.and.
     &ifli(3-jf).eq.0) then
        if(min(iabs(ifl(jf)),iabs(ifl(3))).gt.10) goto 230
        call luifldC(ifl(jf),0,-ifl(3),ifldmp,k(i,2))

c...else generate new flavour including j and i quarks
      else
        if(iabs(ifl(jf)).gt.100) goto 230
        if(iflj(jf).ne.0.and.iabs(ifl(jf)).gt.10) goto 230
        if(iflj(jf).ne.0.and.lrk(jf).eq.1) ifl(jf)=iflf(jf)
        ifl3=ifli(jf)
        if(jf.eq.jr.and.ifl3.eq.0) ifl3=-ifl(3)
        if(ifli(jf).eq.0.and.ifli(3-jf).ne.0) ifl3=-iflf(3)*(-1)**jf
        if((iabs(ifl(jf)).gt.10.or.iflj(jf).ne.0).and.iabs(ifl3).gt.10)
     &  goto 230
        call luifldC(ifl(jf),iflj(jf),ifl3,ifl(3),k(i,2))
        if(k(i,2).eq.0) goto 230
        if(iflj(jf).ne.0.and.mqj*jf.ne.jt) px(jf)=px(jf)+pxj(jf)
        if(iflj(jf).ne.0.and.mqj*jf.ne.jt) py(jf)=py(jf)+pyj(jf)
        if(iflj(jf).ne.0.and.n.gt.ie(jf).and.k(ie(jf)+1,1)/10000.eq.6)
     &  p(ie(jf)+1,2)=i-nrs
        if(ifli(jf).ne.0) p(ie(jf)+1,4)=i-nrs
      endif

c...find masses and transverse momenta
      p(i,5)=ulmassC(1,k(i,2))
      if(jf.eq.jt) pr(jf)=p(i,5)**2+(px(jf)+px(3))**2+(py(jf)+py(3))**2
  470 if(jf.eq.jr) pr(jf)=p(i,5)**2+(px(jf)-px(3))**2+(py(jf)-py(3))**2

c...find common setup of four-vectors
      jq=1
      if(p(in(4)+2,3)*p(in(5)+2,3)*four(in(4),in(5)).lt.p(in(7),3)*
     &p(in(8),3)*four(in(7),in(8))) jq=2
      hc12=four(in(3*jq+1),in(3*jq+2))
      hr1=four(n+nrs,in(3*jq+2))/hc12
      hr2=four(n+nrs,in(3*jq+1))/hc12
      if(in(4).ne.in(7).or.in(5).ne.in(8)) then
        px(3-jq)=-four(n+nrs,in(3*jq+3))-px(jq)
        py(3-jq)=-four(n+nrs,in(3*jq+3)+1)-py(jq)
        pr(3-jq)=p(i+(jt+jq-3)**2-1,5)**2+(px(3-jq)+(2*jq-3)*js*
     &  px(3))**2+(py(3-jq)+(2*jq-3)*js*py(3))**2
      endif

c...solve kinematics for final two (if possible)
      wrem2=wrem2+(px(1)+px(2))**2+(py(1)+py(2))**2
      hd=(sqrt(pr(1))+sqrt(pr(2)))/sqrt(wrem2)
      if(hd.ge.1.) goto 230
      ha=wrem2+pr(jt)-pr(jr)
      if(mst(4).eq.2) prev=0.5*hd**par(27+mst(4))
c        &&&&&&&&&&&&
      if(mst(4).ne.2) then
c                strangely enough, if hd = 0.999972
c               and exponet becomes 2.139834E+08 in the
c               next equation,  absoft execution core dumps.
c               so we put real*8 prev2
c           prev=0.5*hd**(par(27+mst(4))*(pr(1)+pr(2))**2)
         prev2 = (par(27+mst(4))*(pr(1)+pr(2))**2)
         prev = 0.5*hd**prev2
      endif
c       &&&&&&&&&&
      hb=sign(sqrt(max(0.,ha**2-4.*wrem2*pr(jt))),js*(rluC(0)-prev))
      if(max(mod(iabs(ifl(1)),10),mod(iabs(ifl(2)),10)).ge.6)
     &hb=sign(sqrt(max(0.,ha**2-4.*wrem2*pr(jt))),float(js))
      do 480 j=1,4
      p(i-1,j)=(px(jt)+px(3))*p(in(3*jq+3),j)+(py(jt)+py(3))*
     &p(in(3*jq+3)+1,j)+0.5*(hr1*(ha+hb)*p(in(3*jq+1),j)+
     &hr2*(ha-hb)*p(in(3*jq+2),j))/wrem2
  480 p(i,j)=p(n+nrs,j)-p(i-1,j)

c...boost back for rapidly moving system
      if(mbst.eq.1) then
        mst(1)=n+nrs+1
        mst(2)=i
        call luroboC(0.,0.,ps(1)/pebst,ps(2)/pebst,ps(3)/pebst)
        mst(1)=0
        mst(2)=0
      endif

c...mark jets as fragmented, move up particles
      ni=i-nrs
      do 490 i=n+1,n+np
  490 k(k(i,1),1)=k(k(i,1),1)+20000
      do 500 i=n+1,ni
      k(i,1)=k(i+nrs,1)
      k(i,2)=k(i+nrs,2)
      do 500 j=1,5
  500 p(i,j)=p(i+nrs,j)

c...optional ordering along chain, i.e. in rank
      if(mst(22).ge.1) then
        if(2*ni-n.ge.mst(30)-5-mst(31)) then
          mst(24)=mst(24)+1
          mst(25)=1
          if(mst(23).ge.1) return
        endif
        do 510 i=n+1,ni
        k(i-n+ni,1)=k(i,1)
        k(i-n+ni,2)=k(i,2)
        do 510 j=1,5
  510   p(i-n+ni,j)=p(i,j)
        i1=n
        do 530 i=ni+1,2*ni-n
        if(k(i,1).ne.ie(1)) goto 530
        i1=i1+1
        k(i1,1)=k(i,1)
        k(i1,2)=k(i,2)
        do 520 j=1,5
  520   p(i1,j)=p(i,j)
  530   continue
        do 550 i=2*ni-n,ni+1,-1
        if(k(i,1).eq.ie(1)) goto 550
        i1=i1+1
        k(i1,1)=k(i,1)
        k(i1,2)=k(i,2)
        do 540 j=1,5
  540   p(i1,j)=p(i,j)
  550   continue
      endif

c...bring baryons together inside gluon loop
      if(mst(22).ge.2.and.iflstr.ne.0) then
        do 560 i=n+1,n+iflstr
        k(i-n+ni,1)=k(i,1)
        k(i-n+ni,2)=k(i,2)
        do 560 j=1,5
  560   p(i-n+ni,j)=p(i,j)
        do 570 i=n+1,ni
        k(i,1)=k(i+iflstr,1)
        k(i,2)=k(i+iflstr,2)
        do 570 j=1,5
  570   p(i,j)=p(i+iflstr,j)
      endif
      n=ni

      return
      end

c*********************************************************************

      subroutine ludecyC(ip)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludat2C/ktyp(120),pmas(120),pwid(60),kfr(80),cfr(40)
      common/ludat3C/dpar(20),idb(120),cbr(400),kdp(1600)
      dimension iflo(4),ifl1(4),pv(10,5),rord(10),ue(3),be(3)

c//////////
c      real temp
c//////////
c...functions : momentum in two-particle decays and four-product

c//////////  max must be taken to avoid negative sqrt .2001.9/27 kk
      pawt(a,b,c)=sqrt( max((a**2-(b+c)**2)*(a**2-(b-c)**2),0.))/(2.*a)
c//////////
      four(i,j)=p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3)

c...choose decay channel
      ntry=0
      nsav=n
      kfa=iabs(k(ip,2))
      kfs=isign(1,k(ip,2))
  100 rbr=rluC(0)
      if(kfa.le.100) then
        idc=idb(kfa)-1
      else
        call luiflvC(kfa,ifla,iflb,iflc,ksp)
        idc=idb(76+5*ifla+ksp)-1
      endif
  110 idc=idc+1
      if(rbr.gt.cbr(idc)) goto 110

c...start readout of decay channel: matrix element, reset counters
      mmat=iabs(kdp(4*idc-3))/1000
  120 ntry=ntry+1
      if(ntry.gt.1000) then
        mst(24)=mst(24)+1
        mst(25)=5
        if(mst(23).ge.1) return
      endif
      i=n
      np=0
      nq=0
      mbst=0
      if(mmat.ge.5.and.p(ip,4).gt.20.*p(ip,5)) mbst=1
      do 130 j=1,4
      pv(1,j)=0.
  130 if(mbst.eq.0) pv(1,j)=p(ip,j)
      if(mbst.eq.1) pv(1,4)=p(ip,5)
      pv(1,5)=p(ip,5)
      ps=0.
      psq=0.
      nm=0

      do 140 i1=4*idc-3,4*idc
c...read out decay product, convert to standard flavour code
      kp=mod(kdp(i1),1000)
      if(kp.eq.0) goto 140
      if(iabs(kp).le.100) then
        kfp=kfs*kp
        if(mod(ktyp(iabs(kp)),10).eq.0) kfp=kp
      elseif(iabs(kp).lt.590) then
        kfp=kfs*kp
        if(kp.eq.500) kfp=kp
      elseif(iabs(kp).eq.590) then
        if(ksp.le.1) kfp=kfs*(-500+iflb)
        if(ksp.eq.3) kfp=kfs*(500+10*iflc+iflb)
        if(ksp.eq.2.or.ksp.eq.4) kfp=kfs*(500+10*iflb+iflc)
      elseif(iabs(kp).eq.591) then
        call luifldC(-kfs*int(1.+(2.+par(2))*rluC(0)),0,0,kfp,kdump)
        if(pv(1,5).lt.par(22)+2.*ulmassC(2,kfp)) goto 120
        kfp=mod(kfp,100)+isign(500,kfp)
      elseif(iabs(kp).eq.592) then
        kfp=-kfp
      endif

c...add decay product to event record or to iflo list
      if(mmat.ge.6.and.mmat.le.8.and.iabs(kfp).ge.500) then
        nq=nq+1
        iflo(nq)=mod(kfp,500)
        psq=psq+ulmassC(3,iflo(nq))
      elseif(mmat.ge.12.and.np.eq.3) then
        nq=nq-1
        ps=ps-p(i,5)
        k(i,1)=ip
        call luifldC(mod(kfp,500),0,mod(k(i,2),500),ifldmp,k(i,2))
        p(i,5)=ulmassC(1,k(i,2))
        ps=ps+p(i,5)
      else
        i=i+1
        np=np+1
        if(iabs(kfp).ge.500) nq=nq+1
        k(i,1)=ip+10000*(nq-2*(nq/2))
        k(i,2)=kfp
        p(i,5)=ulmassC(1+2*(iabs(kfp)/500),kfp)
        ps=ps+p(i,5)
      endif
  140 continue

  150 if(mmat.ge.6.and.mmat.le.8) then
c...choose decay multiplicity in phase space model
        psp=ps
        cnde=dpar(11)*alog(max((pv(1,5)-ps-psq)/dpar(12),1.1))
        if(mmat.eq.8) cnde=cnde+dpar(13)
  160   ntry=ntry+1
        if(ntry.gt.1000) then
          mst(24)=mst(24)+1
          mst(25)=5
          if(mst(23).ge.1) return
        endif
c//////////
c        temp =-2.*cnde*alog(max(1e-10,rluC(0)))
c        if( temp .lt. 0.) then
c           write(0,*) '2:', temp, ' cnde', cnde
c        endif
c///////////
        gauss=sqrt(-2.*cnde*alog(max(1e-10,rluC(0))))*
     &   sin(par(72)*rluC(0))
        nd=0.5+0.5*np+0.25*nq+cnde+gauss
        if(nd.lt.np+nq/2.or.nd.lt.2.or.nd.gt.10) goto 160
        if(mmat.eq.7.and.nd.eq.2) goto 160

c...form hadrons from flavour content
        do 170 jt=1,4
  170   ifl1(jt)=iflo(jt)
        if(nd.eq.np+nq/2) goto 190
        do 180 i=n+np+1,n+nd-nq/2
        jt=1+int((nq-1)*rluC(0))
        call luifldC(ifl1(jt),0,0,ifl2,k(i,2))
  180   ifl1(jt)=-ifl2
  190   jt=2+2*(nq/4)*int(rluC(0)+0.5)
        if(min(iabs(ifl1(1)),iabs(ifl1(jt))).gt.10.or.(nq.eq.4.and.
     &  min(iabs(ifl1(3)),iabs(ifl1(6-jt))).gt.10)) goto 160
        if(max(iabs(ifl1(1)),iabs(ifl1(jt))).gt.100.or.(nq.eq.4.and.
     &  max(iabs(ifl1(3)),iabs(ifl1(6-jt))).gt.100)) goto 160
        call luifldC(ifl1(1),0,ifl1(jt),ifldmp,k(n+nd-nq/2+1,2))
        if(nq.eq.4) call luifldC(ifl1(3),0,ifl1(6-jt),ifldmp,k(n+nd,2))

c...check that sum of decay product masses not too large
        ps=psp
        do 200 i=n+np+1,n+nd
        k(i,1)=ip
        p(i,5)=ulmassC(1,k(i,2))
  200   ps=ps+p(i,5)
        if(ps+dpar(14).gt.pv(1,5)) goto 160

      elseif(mmat.eq.5.or.mmat.eq.11) then
c...rescale energy to subtract off spectator quark mass
        ps=ps-p(n+np,5)
        pqt=(p(n+np,5)+dpar(15))/pv(1,5)
        do 210 j=1,5
        p(n+np,j)=pqt*pv(1,j)
  210   pv(1,j)=(1.-pqt)*pv(1,j)
        nd=np-1

      else
c...fully specified final states, check mass broadening effects
        if(np.ge.2.and.ps+dpar(14).gt.pv(1,5)) goto 120
        nd=np
      endif

      if(nd.eq.1) then
c...kinematics of one-particle decays
        do 220 j=1,4
  220   p(n+1,j)=p(ip,j)
        goto 430
      endif

c...calculate maximum weight nd-particle decay
      pv(nd,5)=p(n+nd,5)
      if(nd.eq.2) goto 280
      wtmax=1./dpar(nd-2)
      pmax=pv(1,5)-ps+p(n+nd,5)
      pmin=0.
      do 230 il=nd-1,1,-1
      pmax=pmax+p(n+il,5)
      pmin=pmin+p(n+il+1,5)
  230 wtmax=wtmax*pawt(pmax,pmin,p(n+il,5))

c...m-generator gives weight, if rejected try again
  240 rord(1)=1.
      do 260 il1=2,nd-1
      rsav=rluC(0)
      do 250 il2=il1-1,1,-1
      if(rsav.le.rord(il2)) goto 260
  250 rord(il2+1)=rord(il2)
  260 rord(il2+1)=rsav
      rord(nd)=0.
      wt=1.
      do 270 il=nd-1,1,-1
      pv(il,5)=pv(il+1,5)+p(n+il,5)+(rord(il)-rord(il+1))*(pv(1,5)-ps)
  270 wt=wt*pawt(pv(il,5),pv(il+1,5),p(n+il,5))
      if(wt.lt.rluC(0)*wtmax) goto 240

c...perform two-particle decays in respective cm frame
  280 do 300 il=1,nd-1
      pa=pawt(pv(il,5),pv(il+1,5),p(n+il,5))
      ue(3)=2.*rluC(0)-1.
      phi=par(72)*rluC(0)
      ue(1)=sqrt(1.-ue(3)**2)*cos(phi)
      ue(2)=sqrt(1.-ue(3)**2)*sin(phi)
      do 290 j=1,3
      p(n+il,j)=pa*ue(j)
  290 pv(il+1,j)=-pa*ue(j)
      p(n+il,4)=sqrt(pa**2+p(n+il,5)**2)
  300 pv(il+1,4)=sqrt(pa**2+pv(il+1,5)**2)

c...lorentz transform decay products to lab frame
      do 310 j=1,4
  310 p(n+nd,j)=pv(nd,j)
      do 340 il=nd-1,1,-1
      do 320 j=1,3
  320 be(j)=pv(il,j)/pv(il,4)
      ga=pv(il,4)/pv(il,5)
      do 340 i=n+il,n+nd
      bep=be(1)*p(i,1)+be(2)*p(i,2)+be(3)*p(i,3)
      do 330 j=1,3
  330 p(i,j)=p(i,j)+ga*(ga*bep/(1.+ga)+p(i,4))*be(j)
  340 p(i,4)=ga*(p(i,4)+bep)

      if(mmat.eq.1) then
c...matrix elements for omega and phi decays
        wt=(p(n+1,5)*p(n+2,5)*p(n+3,5))**2-(p(n+1,5)*four(n+2,n+3))**2
     &  -(p(n+2,5)*four(n+1,n+3))**2-(p(n+3,5)*four(n+1,n+2))**2
     &  +2.*four(n+1,n+2)*four(n+1,n+3)*four(n+2,n+3)
        if(max(wt*dpar(9)/p(ip,5)**6,0.001).lt.rluC(0)) goto 240

      elseif(mmat.eq.3) then
c...matrix element for s0 -> s1 + v1 -> s1 + s2 + s3 (s scalar,
c...v vector), of form cos**2(theta02) in v1 rest frame
        if(nm.ne.2) then
          im=mod(k(ip,1),10000)
          if(im.eq.0) goto 360
          do 350 il=max(ip-2,im+1),min(ip+2,n)
  350     if(mod(k(il,1),10000).eq.im) nm=nm+1
          call luiflvC(k(im,2),iflam,iflbm,iflcm,kspm)
          if(nm.ne.2.or.kspm.ne.0) goto 360
        endif
        if((p(ip,5)**2*four(im,n+1)-four(ip,im)*four(ip,n+1))**2.le.
     &  rluC(0)*(four(ip,im)**2-(p(ip,5)*p(im,5))**2)*(four(ip,n+1)**2-
     &  (p(ip,5)*p(n+1,5))**2)) goto 280
  360   nm=0

      elseif(mmat.ge.11) then
c...matrix elements for weak decays (only semileptonic for c and b)
        if(mbst.eq.0) wt=four(ip,n+1)*four(n+2,n+3)
        if(mbst.eq.1) wt=p(ip,5)*p(n+1,4)*four(n+2,n+3)
        if(wt.lt.rluC(0)*p(ip,5)*pv(1,5)**3/dpar(10)) goto 240
      endif

      if(mmat.eq.5.or.mmat.eq.11) then
c...scale back energy and reattach spectator
        do 370 j=1,5
  370   pv(1,j)=pv(1,j)/(1.-pqt)
        nd=nd+1
      endif

c...low invariant mass for system with spectator quark gives particle,
c...not two jets, readjust momenta accordingly
      if(mmat.eq.5) then
        if(p(n+2,5)**2+p(n+3,5)**2+2.*four(n+2,n+3).ge.
     &  (par(22)+ulmassC(0,k(n+2,2))+ulmassC(0,k(n+3,2)))**2) goto 430
        k(n+2,1)=ip
        call luifldC(mod(k(n+2,2),500),0,mod(k(n+3,2),500),ifldmp,
     &  k(n+2,2))
        p(n+2,5)=ulmassC(1,k(n+2,2))
        ps=p(n+1,5)+p(n+2,5)
        pv(2,5)=p(n+2,5)
        mmat=0
        nd=2
        goto 280
      elseif(mmat.eq.11) then
        if(p(n+3,5)**2+p(n+4,5)**2+2.*four(n+3,n+4).ge.
     &  (par(22)+ulmassC(0,k(n+3,2))+ulmassC(0,k(n+4,2)))**2) goto 400
        k(n+3,1)=ip
        call luifldC(mod(k(n+3,2),500),0,mod(k(n+4,2),500),ifldmp,
     &  k(n+3,2))
        p(n+3,5)=ulmassC(1,k(n+3,2))
        do 380 j=1,3
  380   p(n+3,j)=p(n+3,j)+p(n+4,j)
        p(n+3,4)=sqrt(p(n+3,1)**2+p(n+3,2)**2+p(n+3,3)**2+p(n+3,5)**2)
        ha=p(n+1,4)**2-p(n+2,4)**2
        hb=ha-(p(n+1,5)**2-p(n+2,5)**2)
        hc=(p(n+1,1)-p(n+2,1))**2+(p(n+1,2)-p(n+2,2))**2+
     &  (p(n+1,3)-p(n+2,3))**2
        hd=(pv(1,4)-p(n+3,4))**2
        he=ha**2-2.*hd*(p(n+1,4)**2+p(n+2,4)**2)+hd**2
        hf=hd*hc-hb**2
        hg=hd*hc-ha*hb
c////////////
c        temp = hg**2+he*hf
c        if(temp .lt. 0.) then
c           write(0,*) '3:', temp, ' hg, he, hf', hg, he, hf
c        endif
c///////////
        hh=(sqrt(hg**2+he*hf)-hg)/(2.*hf)
        do 390 j=1,3
        pcor=hh*(p(n+1,j)-p(n+2,j))
        p(n+1,j)=p(n+1,j)+pcor
  390   p(n+2,j)=p(n+2,j)-pcor
        p(n+1,4)=sqrt(p(n+1,1)**2+p(n+1,2)**2+p(n+1,3)**2+p(n+1,5)**2)
        p(n+2,4)=sqrt(p(n+2,1)**2+p(n+2,2)**2+p(n+2,3)**2+p(n+2,5)**2)
        nd=nd-1
      endif

  400 if(mmat.ge.11.and.iabs(k(n+1,2)).ge.500) then
c...check invariant mass of w jets, may give one particle or start over
        pmr=sqrt(max(0.,p(n+1,5)**2+p(n+2,5)**2+2.*four(n+1,n+2)))
        if(pmr.gt.par(22)+ulmassC(0,k(n+1,2))+ulmassC(0,k(n+2,2)))
     &  goto 410
        call luifldC(mod(k(n+1,2),500),0,-isign(1,k(n+1,2)),ifldmp,kf1)
        call luifldC(mod(k(n+2,2),500),0,-isign(1,k(n+2,2)),ifldmp,kf2)
        psm=ulmassC(0,kf1)+ulmassC(0,kf2)
        if(mmat.le.12.and.pmr.gt.0.2*par(22)+psm) goto 410
        if(mmat.eq.13.and.pmr.gt.dpar(14)+psm) goto 410
        if(nd.eq.4.or.kfa.eq.11) goto 120
        k(n+1,1)=ip
        call luifldC(mod(k(n+1,2),500),0,mod(k(n+2,2),500),ifldmp,
     &  k(n+1,2))
        p(n+1,5)=ulmassC(0,k(n+1,2))
        k(n+2,2)=k(n+3,2)
        p(n+2,5)=p(n+3,5)
        ps=p(n+1,5)+p(n+2,5)
        pv(2,5)=p(n+3,5)
        mmat=0
        nd=2
        goto 280
      endif

  410 if(mmat.eq.13) then
c...phase space decay of partons from w decay
        iflo(1)=mod(k(n+1,2),500)
        iflo(2)=mod(k(n+2,2),500)
        k(n+1,1)=k(n+3,1)
        k(n+1,2)=k(n+3,2)
        do 420 j=1,5
        pv(1,j)=p(n+1,j)+p(n+2,j)
  420   p(n+1,j)=p(n+3,j)
        pv(1,5)=pmr
        n=n+1
        np=0
        nq=2
        ps=0.
        psq=ulmassC(3,iflo(1))+ulmassC(3,iflo(2))
        mmat=6
        goto 150
      endif

c...boost back for rapidly moving particle
  430 n=n+nd
      if(mbst.eq.1) then
        do 440 j=1,3
  440   be(j)=p(ip,j)/p(ip,4)
        ga=p(ip,4)/p(ip,5)
        do 460 i=nsav+1,n
        bep=be(1)*p(i,1)+be(2)*p(i,2)+be(3)*p(i,3)
        do 450 j=1,3
  450   p(i,j)=p(i,j)+ga*(ga*bep/(1.+ga)+p(i,4))*be(j)
  460   p(i,4)=ga*(p(i,4)+bep)
      endif
      k(ip,1)=k(ip,1)+20000

      return
      end

c*********************************************************************

      subroutine luifldC(ifl1,ifl2,ifl3,ifl4,kf)
      common/ludat1C/mst(40),par(80)
      common/ludat2C/ktyp(120),pmas(120),pwid(60),kfr(80),cfr(40)

c...preliminaries, optional enhancements behind heavy quark
      ifla=iabs(ifl1)
      iflb=iabs(ifl2)
      iflc=iabs(ifl3)
      par1=par(1)
      par2=par(2)
      par3=par(3)
      par4=3.*par(4)
      if(ifla.ge.4.and.ifla.lt.10.and.abs(par(16)-1.).gt.0.1) then
        par1=2.5*(0.4*par(1))**(1./par(16))
        par2=par(2)**(1./par(16))
        par3=par(3)**(1./par(16))
        par4=3.*par(4)**(1./par(16))
      endif

c...meson or baryon (from existing diquark or not) to be generated
      iflg=0
      ifli=0
      ifl4=0
      kf=0
      mb=1
      if(ifla.gt.10.or.iflc.gt.10) mb=2
      if(ifla.lt.10.and.iflb.eq.0.and.iflc.eq.0.and.
     &(1.+par1)*rluC(0).lt.1.) mb=0
      if(ifla.lt.10.and.iflb.eq.0.and.(iflc+9)/10.eq.1) mb=0

c...parameter combinations for breaking diquark
      if((ifla.gt.100.or.mb.eq.1).and.par(5).gt.0.) then
        par3m=sqrt(par(3))
        par4m=1./(3.*sqrt(par(4)))
        pardm=par(7)/(par(7)+par3m*par(6))
        pars0=par(5)*(2.+(1.+par2*par3m*par(7))*(1.+par4m))
        pars1=par(7)*pars0/(2.*par3m)+par(5)*(par(6)*(1.+par4m)+
     &  par2*par3m*par(6)*par(7))
        pars2=par(5)*2.*par(6)*par(7)*(par2*par(7)+(1.+par4m)/par3m)
        parsm=max(pars0,pars1,pars2)
        par4=par4*(1.+parsm)/(1.+parsm/(3.*par4m))
      endif

      if(mb.eq.0.or.ifla.gt.100) then
c...flavour for meson, possibly with new quark
        if(mb.eq.0) then
          if(iflc.eq.0) ifl4=isign(1+int((2.+par2)*rluC(0)),-ifl1)
          ifld=max(ifla,iflc+iabs(ifl4))
          ifle=min(ifla,iflc+iabs(ifl4))

        else
c...splitting of diquark into meson pluCs new diquark
  100     iflg=mod(ifla,10)+int(rluC(0)+0.5)*
     &    ((ifla-100)/10-mod(ifla,10))
          iflh=mod(ifla,10)+(ifla-100)/10-iflg
          if((iflg.eq.3.and.rluC(0).gt.pardm).or.
     &(iflh.eq.3.and.rluC(0).lt.pardm)) then
            ifli=iflg
            iflg=iflh
            iflh=ifli
          endif
          ifli=1+int((2.+par2*par3m*par(7))*rluC(0))
          if((iflh.ne.ifli.and.rluC(0).gt.(1.+par4m)/max(2.,1.+par4m))
     &    .or.(iflh.eq.ifli.and.rluC(0).gt.2./max(2.,1.+par4m)))
     &    goto 100
          ifld=max(iflg,ifli)
          ifle=min(iflg,ifli)
          ifl4=isign(10*min(ifli,iflh)+max(ifli,iflh)+9*int(rluC(0)+
     &    1./(1.+par4m))*iabs(ifli-iflh),-ifl1)
          mst(33)=ifli
        endif

c...form meson with spin and flavour mixing for diagonal states
        ksp=int(par(8)+rluC(0))
        if(ifld.eq.3) ksp=int(par(9)+rluC(0))
        if(ifld.ge.4) ksp=int(par(10)+rluC(0))
        if(ifld.ne.ifle) then
          kf=isign(kfr(8*ksp+ifld)+ifle,(ifl1+ifl3+ifl4)*(2*ifld-7))
          if(ifla.gt.100.and.ifli.gt.iflg) kf=-kf
        else
          rfr=rluC(0)
          if(ifld.le.3) kf=23+10*ksp+int(rfr+cfr(6*ksp+2*ifld-1))+
     &    int(rfr+cfr(6*ksp+2*ifld))
          if(ifld.eq.4) kf=26+10*ksp
          if(ifld.ge.5) kf=78+4*ksp+ifld
        endif

      else
  110   if(ifla.lt.10.and.iflb.eq.0.and.iflc.eq.0) then
c...generate diquark flavour
          mb=3
          ifld=ifla
  120     ifle=1+int((2.+par2*par3)*rluC(0))
          iflf=1+int((2.+par2*par3)*rluC(0))
          if(ifle.ge.iflf.and.par4.lt.rluC(0)) goto 120
          if(ifle.lt.iflf.and.par4*rluC(0).gt.1.) goto 120
          ifl4=isign(10*ifle+iflf,ifl1)

        elseif(ifla.lt.10.and.iflb.eq.0) then
c...take diquark flavour from input
          ifld=ifla
          ifle=iflc/10
          iflf=mod(iflc,10)

        elseif(ifla.lt.10) then
c...combine diquark flavour from input (and new-generated quark)
          ifld=iflb
          if(iflc.eq.0) ifl4=isign(1+int((2.+par2)*rluC(0)),ifl1)
          ifle=ifla+int(rluC(0)+0.5)*(iflc+iabs(ifl4)-ifla)
          iflf=ifla+iflc+iabs(ifl4)-ifle

        else
c... generate (or take  from input) quark to go with diquark
          if(iflc.eq.0) ifl4=isign(1+int((2.+par2)*rluC(0)),ifl1)
          ifld=iflc+iabs(ifl4)
          ifle=ifla/10
          iflf=mod(ifla,10)
        endif

c...su(6) factors for formation of baryon, return or try again if fails
        lfr=3+2*((2*(ifle-iflf))/(1+iabs(ifle-iflf)))
        if(ifld.ne.ifle.and.ifld.ne.iflf) lfr=lfr+1
        wt=cfr(2*lfr+11)+par(11)*cfr(2*lfr+12)
        if(mb.eq.1.and.ifle.lt.iflf) wt=wt/3.
        if(mb.eq.1.and.iflb.ne.0) wt=0.75*wt
        if(mb.eq.3.and.par(5).gt.0.) then
          wtdq=pars0
          if(max(ifle,iflf).eq.3) wtdq=pars1
          if(min(ifle,iflf).eq.3) wtdq=pars2
          if(ifle.lt.iflf) wtdq=wtdq/(3.*par4m)
          if((1.+wtdq)*rluC(0).gt.1.) ifl4=ifl4+isign(100,ifl1)
          if(ifle.ge.iflf) wt=wt*(1.+wtdq)/(1.+parsm)
          if(ifle.lt.iflf) wt=wt*(1.+wtdq)/(1.+parsm/(3.*par4m))
        endif
        if(iflb.ne.0.and.wt.lt.rluC(0)) return
        if(iflb.eq.0.and.iflc.eq.0.and.wt.lt.rluC(0)) goto 110

c...form baryon
        iflg=max(ifld,ifle,iflf)
        ifli=min(ifld,ifle,iflf)
        iflh=ifld+ifle+iflf-iflg-ifli
        ksp=2+2*int(1.-cfr(2*lfr+11)+(cfr(2*lfr+11)+par(11)*
     &  cfr(2*lfr+12))*rluC(0))

        if(ksp.eq.2.and.iflg.gt.iflh.and.iflh.gt.ifli) then
c...distinguish lambda- and sigma-like particles
          if(ifle.gt.iflf.and.ifld.ne.iflg) ksp=2+int(0.75+rluC(0))
          if(ifle.lt.iflf.and.ifld.eq.iflg) ksp=3
          if(ifle.lt.iflf.and.ifld.ne.iflg) ksp=2+int(0.25+rluC(0))
        endif

c      &&&&&&&&&&&& kk from version  > uv.6.20
        i = 16*ksp-16+iflg
        j = 16*ksp-8+iflh
        if(i .lt. 1 .or. i .gt. 80) then
           kf = 1
        elseif(j .lt. 1 .or. j .gt. 80 ) then
           kf = 1
        else
c           kf=isign(kfr(16*ksp-16+iflg)+kfr(16*ksp-8+iflh)+ifli,ifl1)
           kf=isign(kfr(i)+kfr(j)+ifli,ifl1)
        endif
c       &&&&&&&&&&&& 
      endif

      return
      end

c*********************************************************************

      function ulmassC(mmass,kf)
      common/ludat1C/mst(40),par(80)
      common/ludat2C/ktyp(120),pmas(120),pwid(60),kfr(80),cfr(40)

      kfa=iabs(kf)
      ulmassC=0.
      kty=0

      if(mmass.le.1.and.kfa.le.100) then
c...ordinary particle masses
        ulmassC=pmas(kfa)
        kty=ktyp(kfa)/10

      elseif(mmass.le.1.and.kfa.lt.500) then
c...heavy hadron masses including chromomagnetic spin-spin interactions
        call luiflvC(kfa,ifla,iflb,iflc,ksp)
        if(ksp.eq.2.and.ifla.eq.iflb) then
          ifld=ifla
          ifla=iflc
          iflc=ifld
        endif
        pma=pmas(100+ifla)
        pmb=pmas(100+iabs(iflb))
        pmc=pmas(100+iflc)
        if(ksp.le.1) ulmassC=pmas(113)+pma+pmb+pmas(115)*pmas(101)**2*
     &  cfr(25+ksp)/(pma*pmb)
        if(ksp.ge.2) ulmassC=pmas(114)+pma+pmb+pmc+pmas(116)*
     &  pmas(101)**2*(cfr(21+3*ksp)/(pma*pmb)+cfr(22+3*ksp)/(pma*pmc)+
     &  cfr(23+3*ksp)/(pmb*pmc))
        kty=ktyp(100+ifla)/10

      else
c...quark and diquark masses: constituent and current algebra-like
        kfa=mod(kfa,100)
        if(kfa.ge.1.and.kfa.le.10) then
          ulmassC=pmas(100+kfa)
          if(mmass.eq.3) ulmassC=ulmassC-pmas(111)
          kty=ktyp(100+kfa)/10
        elseif(kfa.gt.10) then
          ksp=0
          if(kfa/10.ge.mod(kfa,10)) ksp=1
          pma=pmas(100+kfa/10)
          pmb=pmas(100+mod(kfa,10))
          ulmassC=pma+pmb
          if(mmass.eq.3) ulmassC=ulmassC-pmas(112)+
     &     pmas(116)*pmas(101)**2*cfr(25+ksp)/(pma*pmb)
          kty=ktyp(100+max(kfa/10,mod(kfa,10)))/10
        endif
      endif

c...optional mass broadening (truncated breit-wigner shape)
      if(mst(8).eq.1.and.mmass.ge.1.and.kty.ge.1) ulmassC=ulmassC+0.5*
     &pwid(2*kty-1)*tan((2.*rluC(0)-1.)*atan(2.*pwid(2*kty)/
     &pwid(2*kty-1)))

      return
      end

c*********************************************************************

      subroutine luptdiC(ifl,px,py)
      common/ludat1C/mst(40),par(80)

      ifla=iabs(ifl)
c...generate pt according to gaussian, some cases with different widths
      pt=par(12)*sqrt(-alog(max(1e-10,rluC(0))))
      if(mst(32).eq.1) pt=par(13)*pt
      if(ifla.ge.4.and.ifla.lt.10.and.abs(par(16)-1.).gt.0.1) pt=
     &sqrt(par(16))*pt
      if(ifla.eq.93.and.mst(11).eq.1) pt=par(14)*pt
      if(ifla.eq.93.and.mst(11).ne.1) pt=0.
      if(ifla.eq.94) pt=par(15)*pt/par(12)
      phi=par(72)*rluC(0)
      px=pt*cos(phi)
      py=pt*sin(phi)

      return
      end

c*********************************************************************

      subroutine luzdisC(ifl1,ifl3,pr,z)
      common/ludat1C/mst(40),par(80)

      ifla=max(mod(iabs(ifl1),100)/10,mod(iabs(ifl1),10))
      if(ifla.ne.0.and.(mst(4).eq.1.or.(mst(4).eq.3.and.ifla.le.3)))
     &then
c...symmetric scaling function: position of maximum, divide interval
        fa=par(31)
        fb=par(32)*pr
        if(ifl3.ne.0) fa=par(35)
        if(ifl3.ne.0) fb=par(36)*pr
        if(mst(32).eq.1) fa=par(33)
        if(mst(32).eq.1) fb=par(34)*pr
        if(fa.le.0.01) zmax=min(1.,fb)
        if(fa.gt.0.01.and.abs(fa-1.)/fb.le.0.01) zmax=fb/(1.+fb)+
     &  (1.-fa)*fb**2/(1.+fb)**3
        if(fa.gt.0.01.and.abs(fa-1.)/fb.gt.0.01) zmax=0.5*(1.+fb-
     &  sqrt((1.-fb)**2+4.*fa*fb))/(1.-fa)
        if(zmax.lt.0.1) zdiv=2.75*zmax
        if(zmax.gt.0.85) zdiv=zmax-0.6/fb**2+(fa/fb)*alog((0.01+fa)/fb)
c...choice of z, preweighted for peaks at low or high z
  100   z=rluC(0)
        idiv=1
        fpre=1.
        if(zmax.lt.0.1) then
          if(1..lt.rluC(0)*(1.-alog(zdiv))) idiv=2
          if(idiv.eq.1) z=zdiv*z
          if(idiv.eq.2) z=zdiv**z
          if(idiv.eq.2) fpre=zdiv/z
        elseif(zmax.gt.0.85) then
          if(1..lt.rluC(0)*(fb*(1.-zdiv)+1.)) idiv=2
          if(idiv.eq.1) z=zdiv+alog(z)/fb
          if(idiv.eq.1) fpre=exp(fb*(z-zdiv))
          if(idiv.eq.2) z=zdiv+z*(1.-zdiv)
        endif
c...weighting according to correct formula
        if(z.le.fb/(50.+fb).or.z.ge.1.) goto 100
        fval=(zmax/z)*exp(fb*(1./zmax-1./z))
        if(fa.gt.0.01) fval=((1.-z)/(1.-zmax))**fa*fval
        if(fval.lt.rluC(0)*fpre) goto 100

      elseif(ifl1.ne.0) then
c...generate z according to field-feynman, slac, (1-z)**c or z**c
        fc=par(40+ifla)
        if(mst(32).eq.1) fc=par(49)
        if(ifl3.ne.0) fc=par(50)
  110   z=rluC(0)
        if(fc.ge.0..and.fc.le.1.) then
          if(fc.gt.rluC(0)) z=1.-z**(1./3.)
        elseif(fc.gt.-1.) then
          if(-4.*fc*z*(1.-z)**2.lt.rluC(0)*((1.-z)**2-fc*z)**2) goto 110
        else
          if(fc.gt.0.) z=1.-z**(1./fc)
          if(fc.lt.0.) z=z**(-1./fc)
        endif

      else
c...position of j or i quark
  120   z=rluC(0)**(1./(1.+max(par(49+2*ifl3),par(50+2*ifl3))))
        if((1.-z)**min(par(49+2*ifl3),par(50+2*ifl3)).lt.rluC(0))
     &  goto 120
        if(par(50+2*ifl3).gt.par(49+2*ifl3)) z=1.-z
      endif

      return
      end

c*********************************************************************

      subroutine luiflvC(kf,ifla,iflb,iflc,ksp)
      common/ludat2C/ktyp(120),pmas(120),pwid(60),kfr(80),cfr(40)

      kfa=iabs(kf)
      kfs=isign(1,kf)
      ifla=0
      iflb=0
      iflc=0

c...reconstruct spin for hadron
      ksp=-1
      if((kfa.ge.17.and.kfa.le.26).or.kfa.eq.37.or.kfa.eq.38.or.
     &(kfa.ge.83.and.kfa.le.86).or.kfa.ge.101) ksp=0
      if((kfa.ge.27.and.kfa.le.36).or.(kfa.ge.87.and.kfa.le.90).or.
     &kfa.ge.123) ksp=1
      if((kfa.ge.41.and.kfa.le.56).or.kfa.ge.145) ksp=2
      if((kfa.ge.57.and.kfa.le.60).or.kfa.ge.241) ksp=3
      if((kfa.ge.61.and.kfa.le.80).or.kfa.ge.293) ksp=4
      if(kfa.ge.393) ksp=-1

c...reconstruct flavour content for meson
      if((kfa.ge.23.and.kfa.le.26).or.(kfa.ge.33.and.kfa.le.36).or.
     &(kfa.ge.83.and.kfa.le.90)) then
        if(kfa.le.40) ifla=kfa-22-10*ksp
        if(kfa.ge.80) ifla=kfa-78-4*ksp
        iflb=-ifla
      elseif(kfa.eq.37.or.kfa.eq.38) then
        ifla=isign(3,(-1)**int(rluC(0)+0.5))
        iflb=isign(2,-ifla)
      elseif(ksp.eq.0.or.ksp.eq.1) then
  100   ifla=ifla+1
        if(ifla.lt.8.and.kfr(8*ksp+ifla+1).lt.kfa) goto 100
        iflb=-(kfa-kfr(8*ksp+ifla))
        if(ifla.le.3) iflb=-iflb
        if(ifla.le.3) ifla=-ifla

c...reconstruct flavour content for baryon
      elseif(ksp.ge.2.and.ksp.le.4) then
  110   ifla=ifla+1
        if(ifla.lt.8.and.kfr(16*ksp+ifla-15).lt.kfa) goto 110
  120   iflb=iflb+1
        if(iflb.lt.8.and.kfr(16*ksp+iflb-7).lt.kfa-kfr(16*ksp+ifla-16))
     &  goto 120
        iflc=kfa-kfr(16*ksp+ifla-16)-kfr(16*ksp+iflb-8)
      endif

      ifla=kfs*ifla
      iflb=kfs*iflb
      iflc=kfs*iflc

      return
      end

c*********************************************************************

      function luchgeC(kf)
      common/ludat2C/ktyp(120),pmas(120),pwid(60),kfr(80),cfr(40)

      kfa=iabs(kf)
      luchgeC=0

c...calculate 3*charge for particles and partons
      if(kfa.le.100) then
        kty=mod(ktyp(kfa),10)
        if(kty.ge.1) luchgeC=3*kty-6

      elseif(kfa.lt.500) then
        call luiflvC(kfa,ifla,iflb,iflc,ksp)
        luchgeC=3*mod(ktyp(100+ifla),10)-16+
     &  (3*mod(ktyp(100+iabs(iflb)),10)-16)*isign(1,iflb)
        if(iflc.ne.0) luchgeC=luchgeC+3*mod(ktyp(100+iflc),10)-16

      elseif(kfa.le.600) then
        if(mod(kfa,10).ne.0) luchgeC=3*mod(ktyp(100+mod(kfa,10)),10)-16
        if(kfa.gt.510) luchgeC=luchgeC+3*mod(ktyp(50+kfa/10),10)-16

      elseif(kfa.le.700) then
        if(mod(kfa,10).ne.0) luchgeC=3*mod(ktyp(100+mod(kfa,10)),10)-16
      endif

      luchgeC=luchgeC*isign(1,kf)

      return
      end

c*********************************************************************

      subroutine lunameC(kf,chau)
      common/ludat2C/ktyp(120),pmas(120),pwid(60),kfr(80),cfr(40)
      common/ludat4C/chag(50),chaf(100)
      character chau*8,chag*4,chaf*4

      chau=chag(1)//chag(1)
      kfa=iabs(kf)
      kfs=isign(1,kf)

      if(kfa.eq.0) then
c...particle names: ordinary and heavy hadrons
      elseif(kfa.le.100) then
        chau=chaf(kfa)//chag(27+kfs*mod(ktyp(kfa),10))
      elseif(kfa.lt.500) then
        call luiflvC(kfa,ifla,iflb,iflc,ksp)
        if(iflc.eq.0) chau=chag(10+ifla)(1:1)//chag(10-iflb)(1:2)//
     &  chag(35-ksp)(1:1)//chag(27+kfs*(luchgeC(kfa)/3+2))
        if(iflc.ne.0) chau=chag(10+ifla)(1:1)//chag(10+iflb)(1:1)//
     &  chag(10+iflc)(1:1)//chag(30+ksp)(1:1)//chag(27+kfs*
     &  (luchgeC(kfa)/3+2))

c...jet names: gluon, quarks and diquarks; also spectator and phasespace
      elseif(kfa.lt.590) then
        ifla=max(kfa/10-50,kfa-10*(kfa/10))
        iflb=min(kfa/10-50,kfa-10*(kfa/10))
        if(iflb.eq.0) chau=chag(10+kfs*ifla)//chag(22)
        ksp=32
        if(kfa/10-50.lt.ifla) ksp=33
        if(iflb.ne.0) chau=chag(10+ifla)(1:1)//chag(10+kfs*iflb)(1:2)//
     &  chag(ksp)(1:1)//chag(22)
      elseif(kfa.le.600) then
        chau=chag(kfa-571)//chag(22)

      elseif(kfa.lt.700) then
c...hadron jets: j and i quark
        ifla=isign(kfa/10-60,kf)
        iflb=isign(mod(kfa,10),kf)
        if(ifla.ne.0) chau(1:4)=chag(10+ifla)(1:2)//chag(37)(1:2)
        if(iflb.ne.0) chau(5:8)=chag(10+iflb)(1:2)//chag(38)(1:2)
      endif

      return
      end

c*********************************************************************

      function ulanglC(x,y)
      common/ludat1C/mst(40),par(80)

      ulanglC=0.
c...reconstruct the angle from x and y coordinates
      r=sqrt(x**2+y**2)
      if(r.lt.1e-20) return
      if(abs(x)/r.lt.0.8) then
        ulanglC=sign(acos(x/r),y)
      else
        ulanglC=asin(y/r)
        if(x.lt.0..and.ulanglC.ge.0.) then
          ulanglC=par(71)-ulanglC
        elseif(x.lt.0.) then
          ulanglC=-par(71)-ulanglC
        endif
      endif

      return
      end

c*********************************************************************

      block data ludataC
      common/ludat1C/mst(40),par(80)
      common/ludat2C/ktyp(120),pmas(120),pwid(60),kfr(80),cfr(40)
      common/ludat3C/dpar(20),idb(120),cbr(400),kdp(1600)
      common/ludat4C/chag(50),chaf(100)
      character chag*4,chaf*4

c...ludat1C, containing status codes and most parameters
      data mst/
     1    0,    0,    0,    1,    1,    0,    2,    0,    0,    1,
     2    0,    1,   10,    0,    0,    0,    0,    0,    1,    6,
     3    1,    0,    1,    0,    0,    0,    0,    0,    0, 2000,
     4    0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      data par/
     1 0.10, 0.30, 0.40, 0.05, 0.50, 0.50, 0.50, 0.50, 0.60, 0.75,
     2  1.0, 0.35,  1.0,  1.0,   0.,  1.0,  1.0,  2.0,   0.,   0.,
     3 0.10,  1.0,  0.8,  1.5,  0.8,  2.0,  0.2,  2.5,  0.6,  2.5,
     4  0.5,  0.9,  0.5,  0.9,  0.5,  0.9,  1.0,   0.,   0.,   0.,
     5 0.77, 0.77, 0.77,   0.,   0.,   0.,   0.,   0.,  1.0, 0.77,
     6  1.0,  1.0,   0.,  1.0,   0.,   0.,   0.,   0., 0.09, 0.01,
     7   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     8  3.1415927,  6.2831854,   0.,0.001,   0., 5*0./

c...ludat2C, with particle data and flavour treatment
      data ktyp/0,10,23,3*0,1,2,1,2,1,2,1,2,2*0,2*3,2*2,2*3,4*0,33,43,
     &52,2,2*3,60,70,80,5*0,3,2,3,2,1,2,1,4,3,2,3,2*2,4,2*3,2,2*3,2,94,
     &103,112,121,133,142,151,162,171,1,4,3,2,3,2*2,4,2*3,4,10*0,2,3,
     &8*0,6,2*5,6,5,6,5,6,12*0/
      data pmas/0.,94.,83.,15.,2*0.,.00051,0.,.1057,0.,1.7842,0.,60.,
     &3*0.,.1396,.4937,.4977,1.8646,1.8693,1.9705,.135,.5488,.9576,
     &2.981,.7714,.8921,.8965,2.0072,2.0101,2.11,.7717,.7826,1.0195,
     &3.0969,2*.4977,2*0.,.9383,.9396,1.1894,1.1925,1.1973,1.3149,
     &1.3213,3*2.44,2*2.55,2.74,2*3.63,3.81,1.1156,2.2812,2*2.46,1.23,
     &1.231,1.232,1.233,1.3828,1.3837,1.3872,1.5318,1.535,1.6724,3*2.5,
     &2*2.63,2.8,2*3.69,3.85,4.9,2*0.,9.4,77.99001,118.,397.,9.46,78.,
c                                    !101
     &118.,397.,5000.,300.,900.,7*0.,2*.325,.5,1.6,5.,40.,60.,200.,
     &2*0.,.2,.1,0.,.11,.16,.048,4*0./
      data pwid/2.8,20.,2.8,20.,.148,.4,.05,.2,.052,.2,.153,.4,.01,.1,
     &.004,.015,.115,.14,.115,.14,.115,.14,.115,.14,.036,.035,.036,
     &.035,.039,.04,.009,.05,.01,.05,26*0./
      data kfr/0,16,17,19,100,104,109,115,0,26,27,29,122,126,131,137,
     &0,40,42,47,144,158,178,205,0,1,3,6,10,15,21,28,0,0,56,57,240,
     &246,256,271,0,0,1,3,6,10,15,21,60,61,64,70,292,307,328,356,
     &0,1,3,6,10,15,21,28,16*0/
      data cfr/0.5,0.25,0.5,0.25,1.,0.5,0.5,0.,0.5,0.,1.,1.,0.75,0.,
     &0.5,0.,0.,1.,0.1667,0.3333,0.0833,0.6667,0.1667,0.3333,-3.,1.,
     &-2.,-2.,1.,0.,0.,-3.,1.,1.,1.,5*0./

c...ludat3C, devoted to particle decay parameters and data
      data dpar/2.,5.,15.,60.,250.,1500.,1.2e4,1.2e5,150.,16.,4.5,0.7,
     &0.,0.003,0.5,5*0./
      data idb/0,1,13,19,6*0,23,0,30,0,36,3*0,37,39,70,88,92,94,100,
     &106,107,108,110,112,114,117,118,119,122,126,129,5*0,131,133,134,
     &135,136,137,138,139,140,141,142,145,148,151,154,156,161,164,168,
     &169,171,173,174,177,180,183,185,187,190,191,192,193,194,195,196,
     &197,198,199,2*0,202,203,204,205,207,215,227,240,10*0,242,249,
     &3*250,5*257,5*264,5*271/
      data (cbr(j),j=1,200)/.03,.088,.118,.176,.206,.264,.375,.518,
     &.661,.772,.915,1.,.08,.16,.24,.5,.76,1.,.82,.91,.99,1.,.175,.35,
     &.451,.669,.676,.693,1.,.1,.2,.3,.6,.9,2*1.,.5,1.,.07,.14,.194,
     &.265,.364,.464,.504,.527,.528,.544,.559,.574,.579,.584,.604,.634,
     &.649,.655,.666,.676,.686,.688,.69,.692,.694,.696,.699,.702,.707,
     &.92,1.,.182,.364,.405,.465,.525,.605,.607,.61,.615,.622,.632,
     &.647,.659,.664,.67,.676,.681,1.,.1,.2,.7,1.,.988,1.,.389,.708,
     &.945,.994,.999,1.,.427,.652,.952,.979,.998,3*1.,.667,1.,.667,1.,
     &.515,1.,.49,.83,3*1.,.896,.983,1.,.495,.838,.987,1.,.069,.138,1.,
     &.686,1.,.516,10*1.,.15,.3,1.,.15,.3,1.,.15,.3,1.,.15,.3,1.,.642,
     &1.,.045,.09,.6,.8,1.,.15,.3,1.,.05,.1,.6,2*1.,.67,1.,.33,2*1.,
     &.88,.94,1.,.88,.94,1.,.88,.94,1.,.33,1.,.67,1.,.678,.914,10*1.,
     &.15,.3/
      data (cbr(j),j=201,400)/4*1.,.5,1.,.03,.06,.09,.135,.15,.165,.21,
     &1.,.02,.04,.06,.08,.1,.12,.17,.22,.27,.32,.37,1.,.02,.04,.06,.08,
     &.1,.12,.17,.22,.27,.32,.37,.42,1.,.5,1.,.112,.224,.274,.77,.85,
     &.99,2*1.,.112,.224,.274,.77,.85,.99,1.,.12,.24,.34,.67,.93,.97,
     &1.,.11,.22,.33,.64,.94,.97,2*1.,129*0./
      data (kdp(j),j=1,320)/7,-7,2*0,8,-8,2*0,9,-9,2*0,10,-10,2*0,11,
     &-11,2*0,12,-12,2*0,501,-501,2*0,502,-502,2*0,503,-503,2*0,504,
     &-504,2*0,505,-505,2*0,506,-506,2*0,-7,8,2*0,-9,10,2*0,-11,12,2*0,
     &501,-502,2*0,504,-503,2*0,506,-505,2*0,505,-505,2*0,11,-11,2*0,
     &504,-504,2*0,503,-503,2*0,-12008,7,12,0,-12010,9,12,0,12,-17,2*0,
     &12,-27,2*0,12,-18,2*0,12,-28,2*0,-13501,502,12,0,-12008,7,14,0,
     &-12010,9,14,0,-12012,11,14,0,-12501,502,14,0,-12504,503,14,0,
     &-12506,505,14,0,6591,592,2*0,37,3*0,38,3*0,-12007,8,503,-501,
     &-12009,10,503,-501,-18,17,2*0,-28,17,2*0,-18,27,2*0,-28,27,2*0,
     &-19,23,2*0,-29,23,2*0,-19,33,2*0,-29,33,2*0,-19,24,2*0,-29,24,
     &2*0,-19,25,2*0,-29,25,2*0,-19,34,2*0,-29,34,2*0,-19,35,2*0,-18,
     &18,2*0,-28,18,2*0,-18,28,2*0,-28,28,2*0,-19,19,2*0,-29,19,2*0,
     &-19,29,2*0,-29,29,2*0,17,-17,2*0,27,-17,2*0,17,-27,2*0,27,-27,
     &2*0,7501,-502,503,-501,7503,-502,2*0,-12007,8,503,-502,-12009,10,
     &503,-502,-19,17,2*0,-29,17,2*0,-19,27,2*0,-29,27,2*0,24,17,2*0,
     &24,27,2*0,25,17,2*0,25,27,2*0,35,17,2*0/
      data (kdp(j),j=321,640)/35,27,2*0,18,-19,2*0,18,-29,2*0,28,-19,
     &2*0,28,-29,2*0,2*17,-17,0,7501,-502,503,-502,-12007,8,503,-503,
     &-12009,10,503,-503,6501,-502,503,-503,6501,-502,2*0,2*1,2*0,1,7,
     &-7,0,2*1,2*0,3*23,0,17,-17,23,0,1,17,-17,0,1,7,-7,0,17,-17,7,-7,
     &17,-17,24,0,2*23,24,0,1,33,2*0,1,34,2*0,2*1,2*0,3*23,0,8591,592,
     &2*0,3017,23,2*0,3019,17,2*0,3018,23,2*0,3018,-17,2*0,3019,23,2*0,
     &3020,23,2*0,20,1,2*0,3020,17,2*0,3021,23,2*0,21,1,2*0,22,1,2*0,
     &3017,-17,2*0,1017,-17,23,0,1,23,2*0,3017,-17,2*0,3018,-18,2*0,
     &3037,38,2*0,1017,-17,23,0,1,24,2*0,7,-7,2*0,9,-9,2*0,8591,592,
     &2*0,17,-17,2*0,2*23,2*0,41,23,2*0,42,17,2*0,57,1,2*0,42,-17,2*0,
     &57,23,2*0,57,-17,2*0,58,17,2*0,58,23,2*0,58,-17,2*0,59,1,2*0,60,
     &1,2*0,-12007,8,503,533,-12009,10,503,533,6501,-502,503,533,
     &-12007,8,503,541,-12009,10,503,541,6501,-502,503,541,-12007,8,
     &503,542,-12009,10,503,542,6501,-502,503,542,-12007,8,503,543,
     &-12009,10,503,543,6501,-502,503,543,41,-17,2*0,42,23,2*0,-12007,
     &8,503,512,-12009,10,503,512,6501,-502,503,512,6531,501,2*0,6513,
     &501,2*0/
      data (kdp(j),j=641,960)/-12007,8,503,513,-12009,10,503,513,6501,
     &-502,503,513,-12007,8,503,523,-12009,10,503,523,6501,-502,503,
     &523,6533,501,2*0,41,17,2*0,41,23,2*0,42,17,2*0,41,-17,2*0,42,23,
     &2*0,42,-17,2*0,57,17,2*0,43,23,2*0,44,17,2*0,57,23,2*0,43,-17,
     &2*0,45,17,2*0,57,-17,2*0,44,-17,2*0,45,23,2*0,46,23,2*0,47,17,
     &2*0,46,-17,2*0,47,23,2*0,57,-18,2*0,46,-17,2*0,47,23,2*0,58,17,
     &2*0,58,23,2*0,58,-17,2*0,59,1,2*0,60,1,2*0,53,1,2*0,54,1,2*0,55,
     &1,2*0,56,1,2*0,-12007,8,503,544,-12009,10,503,544,6501,-502,503,
     &544,2*500,2*0,2*500,2*0,2*500,2*0,5003,507,-508,0,-5003,-507,508,
     &0,7,-7,2*0,9,-9,2*0,11,-11,2*0,501,-501,2*0,502,-502,2*0,503,
     &-503,2*0,504,-504,2*0,2*500,2*0,7,-7,2*0,8,-8,2*0,9,-9,2*0,10,
     &-10,2*0,11,-11,2*0,12,-12,2*0,501,-501,2*0,502,-502,2*0,503,-503,
     &2*0,504,-504,2*0,505,-505,2*0,2*500,2*0,7,-7,2*0,8,-8,2*0,9,-9,
     &2*0,10,-10,2*0,11,-11,2*0,12,-12,2*0,501,-501,2*0,502,-502,2*0,
     &503,-503,2*0,504,-504,2*0,505,-505,2*0,506,-506,2*0,2*500,2*0,
     &5003,507,-508,0/
      data (kdp(j),j=961,1280)/-5003,-507,508,0,-12008,7,504,590,
     &-12010,9,504,590,-12012,11,504,590,-13501,502,504,590,-13501,504,
     &502,590,-13504,503,504,590,-13504,504,503,590,6001,505,590,0,
     &-12008,7,504,590,-12010,9,504,590,-12012,11,504,590,-13501,502,
     &504,590,-13501,504,502,590,-13504,503,504,590,-13504,504,503,590,
     &-11007,8,505,590,-11009,10,505,590,-11011,12,505,590,-11502,501,
     &505,590,-11503,504,505,590,-11502,505,501,590,-11503,505,504,590,
     &-11008,7,506,590,-11010,9,506,590,-11012,11,506,590,-11501,502,
     &506,590,-11504,503,506,590,-11501,506,502,590,-11504,506,503,590,
     &5003,507,590,197*0/
      data (kdp(j),j=1281,1600)/320*0/

c...ludat4C, containing character strings
      data chag/' ','ha','la','ta','ba','ca','sa','da','ua','g','u','d',
     &'s','c','b','t','l','h','spec','qra','qbra','jet','b--','b-','b',
     &'b+',' ','-','0','+','++','1','0','*',' ','dfbv','jq','iq',
     &' ',' ','stab','unst',8*' '/
      data chaf/'gamm','z0','w','higg','ga/z',' ','e','nue','mu',
     &'numu','tau','nuta','chi','nuch','phas',' ','pi',2*'k',2*'d','f',
     &'pi0','eta','eta''','etac','rho',2*'k*',2*'d*','f*','rho0',
     &'omeg','phi','jpsi','k0s','k0l',2*' ','p','n',3*'sig',2*'xi',
     &3*'sic','csu1','csd1','css1','ccu1','ccd1','ccs1','lam','lamc',
     &'csu0','csd0',4*'delt',3*'sig*',2*'xi*','ome*',3*'sic*','csu*',
     &'csd*','css*','ccu*','ccd*','ccs*','ccc*',2*' ','etab','etat',
     &'etal','etah','upsi','phit','phil','phih','r','higg','z''0',
     &7*' '/

      end

c*********************************************************************
c***  jetset version 6.3, e+e- part  *********************************

      subroutine lueevtC(ifl,ecm)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)

c...printout or resetting of statistics
      if(ifl.ge.10) call luestaC(ifl,0,0,0.,0,0.,0)
      if(ifl.ge.10) return

c...initialize total cross section
      if(mst(19).ge.1) call lulistC(-1)
      if(mste(9).gt.0.and.(mste(9).ge.2.or.abs(ecm-pare(51)).ge.
     &pare(17).or.10*mste(2)+ifl.ne.mste(36))) call luxtotC(ifl,ecm,
     &xtot)
      if(mste(9).ge.3) mste(9)=1

c...add initial e+e- to event record (documentation only)
      ntry=0
  100 ntry=ntry+1
      nc=0
      if(mste(30).ge.2) then
        nc=nc+2
        call lupartC(nc-1,7,0.5*ecm,0.,0.)
        k(nc-1,1)=40000
        call lupartC(nc,-7,0.5*ecm,par(71),0.)
        k(nc,1)=40000
      endif

c...radiative photon (in initial state)
      mk=0
      ecmc=ecm
      if(mste(7).ge.1.and.mste(9).ge.1) call luradkC(ecm,mk,pak,
     &thek,phik,alpk)
      if(mk.eq.1) ecmc=sqrt(ecm*(ecm-2.*pak))
      if(mste(30).ge.1.and.mk.eq.1) then
        nc=nc+1
        call lupartC(nc,1,pak,thek,phik)
        k(nc,1)=min(mste(30)/2,1)
      endif

c...virtual exchange boson (gamma or z0)
      if(mste(30).ge.3) then
        nc=nc+1
        kf=1
        if(mste(2).eq.2) kf=5
        mst(9)=1
        p(nc,5)=ecmc
        call lupartC(nc,kf,ecmc,0.,0.)
        k(nc,1)=50001
        mst(9)=0
      endif

c...choice of flavour and jet configuration
      call luxiflC(ifl,ecm,ecmc,iflc)
      if(iflc.eq.0) goto 100
      call luxjetC(ecmc,njet,cut)
      ifln=0
      if(njet.eq.4) call lux4jtC(njet,cut,iflc,ecmc,ifln,x1,x2,x4,
     &x12,x14)
      if(njet.eq.3) call lux3jtC(njet,cut,iflc,ecmc,x1,x3)
      if(njet.eq.2) mste(35)=1

c...fill jet configuration and origin
      if(njet.eq.2) call lu2jetC(nc+1,iflc,-iflc,ecmc)
      if(njet.eq.3) call lu3jetC(nc+1,iflc,-iflc,ecmc,x1,x3)
      if(njet.eq.4) call lu4jetC(nc+1,iflc,-ifln,ifln,-iflc,ecmc,
     &x1,x2,x4,x12,x14)
      do 110 ip=nc+1,n
  110 k(ip,1)=k(ip,1)+min(mste(30)/2,1)+(mste(30)/3)*(nc-1)

c...angular orientation according to matrix element
      if(mste(6).eq.1) then
        call luxdifC(nc,njet,iflc,ecmc,chi,the,phi)
        mst(1)=nc+1
        call luroboC(0.,chi,0.,0.,0.)
        call luroboC(the,phi,0.,0.,0.)
        mst(1)=0
      endif

c...rotation and boost from radiative photon
      if(mk.eq.1) then
        bek=-pak/(ecm-pak)
        mst(1)=nc+1-mste(30)/3
        call luroboC(0.,-phik,0.,0.,0.)
        call luroboC(alpk,0.,bek*sin(thek),0.,bek*cos(thek))
        call luroboC(0.,phik,0.,0.,0.)
        mst(1)=0
      endif
      nc=n

      if(mste(1).ge.3) then
c...prepare event record for parton shower evolution
        k(n+1,1)=k(n,1)
        k(n+1,2)=k(n,2)
        do 120 j=1,5
  120   p(n+1,j)=p(n,j)
        n=n+2
        do 130 i=n-2,n,2
        k(i,1)=70000+i-1
        k(i,2)=1000+mod(k(i-1,1),1000)
        do 130 j=1,5
  130   p(i,j)=0.
        p(n-2,1)=n-1
        p(n,2)=n-3

c...generate parton shower, rearrange along strings and check
        mste(13)=mste(13)+1
        call lushowC(n-3,n-1,ecmc)
        mste(13)=mste(13)-1
        njet=0
        ifln=-2
        do 140 i=1,n
        if(k(i,1).lt.20000.and.iabs(k(i,2)).ge.500) njet=njet+1
  140   if(k(i,1).lt.20000.and.iabs(k(i,2)).ge.501) ifln=ifln+1
        mst12s=mst(12)
        if(mste(5).eq.-1) mst(12)=0
        mst21s=mst(21)
        if(mste(30).le.3) mst(21)=1
        if(mste(30).ge.4) mst(21)=0
        if(mste(5).ge.0) mst(26)=0
        call luprepC
        mst(12)=mst12s
        mst(21)=mst21s
        if(mste(5).ge.0.and.mst(26).ne.0) goto 100
        nc=n+1
  150   nc=nc-1
        if(iabs(k(nc,2)).lt.500.and.mod(k(nc,1),10000).gt.0) then
          if(k(mod(k(nc,1),10000),1)/20000.eq.1) goto 150
        endif
      endif

c...event generation and statistics
      if(mste(5).eq.1) call luexecC
      if(njet.eq.4.and.ifln.ne.0) njet=-4
      if(mste(31).ne.0) call luestaC(iflc,njet,nc,ecm,mk,ecmc,ntry)

      return
      end

c*********************************************************************

      subroutine luxtotC(ifl,ecm,xtot)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)

c...status, alpha-strong, calculate z0 width for mste(2)=3
      pare(51)=ecm
      mste(36)=10*mste(2)+ifl
      alspi=ulalpsC(ecm**2)/par(71)
      rqcd=1.
      if(iabs(mste(1)).eq.1) rqcd=1.+alspi
      if(iabs(mste(1)).ge.2) rqcd=1.+alspi+pare(65)*alspi**2
      if(mste(2).ge.3) then
        rva=3.*(3.+(4.*pare(5)-1.)**2)+6.*rqcd*(2.+(1.-8.*pare(5)/
     &  3.)**2+(4.*pare(5)/3.-1.)**2)
        do 100 iflc=5,6
        vq=1.
        if(mod(mste(3),2).eq.1) vq=sqrt(max(0.,1.-(2.*ulmassC(2,iflc)/
     &  ecm)**2))
        if(iflc.eq.5) vf=4.*pare(5)/3.-1.
        if(iflc.eq.6) vf=1.-8.*pare(5)/3.
  100   rva=rva+3.*rqcd*(0.5*vq*(3.-vq**2)*vf**2+vq**3)
        pare(7)=pare(4)*pare(6)*rva/(48.*pare(5)*(1.-pare(5)))
      endif

c...calculate propagator and related constants for qfd case
      poll=1.-pare(11)*pare(12)
      if(mste(2).ge.2) then
        sff=1./(16.*pare(5)*(1.-pare(5)))
        sfw=ecm**4/((ecm**2-pare(6)**2)**2+(pare(6)*pare(7))**2)
        sfi=sfw*(1.-(pare(6)/ecm)**2)
        ve=4.*pare(5)-1.
        sf1i=sff*(ve*poll+pare(12)-pare(11))
        sf1w=sff**2*((ve**2+1.)*poll+2.*ve*(pare(12)-pare(11)))
        hf1i=sfi*sf1i
        hf1w=sfw*sf1w
      endif

c...loop over different flavours: charge, velocity
      rtot=0.
      rqq=0.
      rqv=0.
      rva=0.
      do 110 iflc=1,max(mste(4),ifl)
      if(ifl.gt.0.and.iflc.ne.ifl) goto 110
      pmq=ulmassC(2,iflc)
      if(ecm.lt.2.*pmq+pare(10)) goto 110
      qf=2./3.
      if(iflc.eq.2.or.iflc.eq.3.or.iflc.eq.5.or.iflc.eq.7) qf=-1./3.
      vq=1.
      if(mod(mste(3),2).eq.1) vq=sqrt(1.-(2.*pmq/ecm)**2)

c...calculate r and sum of charges for qed or qfd case
      rqq=rqq+3.*qf**2*poll
      if(mste(2).le.1) then
        rtot=rtot+3.*0.5*vq*(3.-vq**2)*qf**2*poll
      else
        vf=sign(1.,qf)-4.*qf*pare(5)
        rqv=rqv-6.*qf*vf*sf1i
        rva=rva+3.*(vf**2+1.)*sf1w
        rtot=rtot+3.*(0.5*vq*(3.-vq**2)*(qf**2*poll-2.*qf*vf*hf1i+
     &  vf**2*hf1w)+vq**3*hf1w)
      endif
  110 continue
      rsum=rqq
      if(mste(2).ge.2) rsum=rqq+sfi*rqv+sfw*rva

c...calculate cross section including qcd corrections
      pare(41)=rqq
      pare(42)=rtot
      pare(43)=rtot*rqcd
      pare(44)=pare(43)
      pare(45)=pare(41)*86.8/ecm**2
      pare(46)=pare(42)*86.8/ecm**2
      pare(47)=pare(43)*86.8/ecm**2
      pare(48)=pare(47)
      pare(57)=rsum*rqcd
      pare(58)=0.
      pare(59)=0.
      xtot=pare(47)
      if(mste(7).le.0) return

c...virtual cross section, soft and hard radiative x-sect in qed case
      xkl=pare(15)
      xku=min(pare(16),1.-(2.*pare(10)/ecm)**2)
      ale=2.*alog(ecm/ulmassC(0,7))-1.
      sigv=ale/3.+2.*alog(ecm**2/(ulmassC(0,9)*ulmassC(0,11)))/3.-4./3.+
     &1.526*alog(ecm**2/0.932)
      if(mste(2).le.1) then
        sigv=1.5*ale-0.5+par(71)**2/3.+2.*sigv
        sigs=ale*(2.*alog(xkl)-alog(1.-xkl)-xkl)
        sigh=ale*(2.*alog(xku/xkl)-alog((1.-xku)/(1.-xkl))-(xku-xkl))

c...ditto in qfd case, total x-sect and fraction hard photon events
      else
        szm=1.-(pare(6)/ecm)**2
        szw=pare(6)*pare(7)/ecm**2
        pare(61)=-rqq/rsum
        pare(62)=-(rqq+rqv+rva)/rsum
        pare(63)=(rqv*(1.-0.5*szm-sfi)+rva*(1.5-szm-sfw))/rsum
        pare(64)=(rqv*szw**2*(1.-2.*sfw)+rva*(2.*sfi+szw**2-4.+3.*szm-
     &  szm**2))/(szw*rsum)
        sigv=1.5*ale-0.5+par(71)**2/3.+((2.*rqq+sfi*rqv)/rsum)*sigv+
     &  (szw*sfw*rqv/rsum)*par(71)*20./9.
        sigs=ale*(2.*alog(xkl)+pare(61)*alog(1.-xkl)+pare(62)*xkl+
     &  pare(63)*alog(((xkl-szm)**2+szw**2)/(szm**2+szw**2))+
     &  pare(64)*(atan((xkl-szm)/szw)-atan(-szm/szw)))
        sigh=ale*(2.*alog(xku/xkl)+pare(61)*alog((1.-xku)/(1.-xkl))+
     &  pare(62)*(xku-xkl)+pare(63)*alog(((xku-szm)**2+szw**2)/
     &  ((xkl-szm)**2+szw**2))+pare(64)*(atan((xku-szm)/szw)-
     &  atan((xkl-szm)/szw)))
      endif
      pare(60)=sigh/(par(71)/pare(4)+sigv+sigs+sigh)
      pare(57)=rsum*(1.+(pare(4)/par(71))*(sigv+sigs+sigh))*rqcd
      pare(44)=pare(57)
      pare(48)=pare(44)*86.8/ecm**2
      xtot=pare(48)

      return
      end

c*********************************************************************

      subroutine luradkC(ecm,mk,pak,thek,phik,alpk)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)

c...function: cumulative hard photon spectrum in qfd case
      fxk(xk)=2.*alog(xk)+pare(61)*alog(1.-xk)+pare(62)*xk+pare(63)*
     &alog((xk-1.+(pare(6)/ecm)**2)**2+(pare(6)*pare(7)/ecm**2)**2)+
     &pare(64)*atan((ecm**2*(xk-1.)+pare(6)**2)/(pare(6)*pare(7)))

c...determine whether radiative photon or not
      mk=0
      pak=0.
      if(pare(60).lt.rluC(0)) return
      mk=1

c...find photon momentum in qed case
      xkl=pare(15)
      xku=min(pare(16),1.-(2.*pare(10)/ecm)**2)
      if(mste(2).le.1) then
  100   xk=1./(1.+(1./xkl-1.)*((1./xku-1.)/(1./xkl-1.))**rluC(0))
        if(1.+(1.-xk)**2.lt.2.*rluC(0)) goto 100

c...ditto in qfd case, by numerical inversion of integrated spectrum
      else
        fxkl=fxk(xkl)
        fxku=fxk(xku)
        dfxk=1e-4*(fxku-fxkl)
        fxkr=fxkl+rluC(0)*(fxku-fxkl)
        nxk=0
  110   nxk=nxk+1
        xk=0.5*(xkl+xku)
        fxkv=fxk(xk)
        if(fxkv.gt.fxkr) then
          xku=xk
          fxku=fxkv
        else
          xkl=xk
          fxkl=fxkv
        endif
        if(nxk.lt.15.and.fxku-fxkl.gt.dfxk) goto 110
        xk=xkl+(xku-xkl)*(fxkr-fxkl)/(fxku-fxkl)
      endif
      pak=0.5*ecm*xk

c...photon polar and azimuthal angle
      pme=2.*(ulmassC(0,7)/ecm)**2
  120 cthm=pme*(2./pme)**rluC(0)
      if(1.-(xk**2*cthm*(1.-0.5*cthm)+2.*(1.-xk)*pme/max(pme,
     &cthm*(1.-0.5*cthm)))/(1.+(1.-xk)**2).lt.rluC(0)) goto 120
      cthe=1.-cthm
      if(rluC(0).gt.0.5) cthe=-cthe
      sthe=sqrt(max(0.,(cthm-pme)*(2.-cthm)))
      thek=ulanglC(cthe,sthe)
      phik=par(72)*rluC(0)

c...rotation angle for hadronic system
      sgn=1.
      if(0.5*(2.-xk*(1.-cthe))**2/((2.-xk)**2+(xk*cthe)**2).gt.
     &rluC(0)) sgn=-1.
      alpk=asin(sgn*sthe*(xk-sgn*(2.*sqrt(1.-xk)-2.+xk)*cthe)/
     &(2.-xk*(1.-sgn*cthe)))

      return
      end

c*********************************************************************

      subroutine luxiflC(ifl,ecm,ecmc,iflc)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)

c...calculate maximum weight in qed or qfd case
      if(mste(2).le.1) then
        rfmax=4./9.
      else
        poll=1.-pare(11)*pare(12)
        sff=1./(16.*pare(5)*(1.-pare(5)))
        sfw=ecmc**4/((ecmc**2-pare(6)**2)**2+(pare(6)*pare(7))**2)
        sfi=sfw*(1.-(pare(6)/ecmc)**2)
        ve=4.*pare(5)-1.
        hf1i=sfi*sff*(ve*poll+pare(12)-pare(11))
        hf1w=sfw*sff**2*((ve**2+1.)*poll+2.*ve*(pare(12)-pare(11)))
        rfmax=max(4./9.*poll-4./3.*(1.-8.*pare(5)/3.)*hf1i+
     &  ((1.-8.*pare(5)/3.)**2+1.)*hf1w,1./9.*poll+2./3.*
     &  (-1.+4.*pare(5)/3.)*hf1i+((-1.+4.*pare(5)/3.)**2+1.)*hf1w)
      endif

c...choose flavour, gives charge and velocity
  100 iflc=ifl
      if(ifl.le.0) iflc=1+int(mste(4)*rluC(0))
      pmq=ulmassC(2,iflc)
      if(ecm.lt.2.*pmq+pare(10)) goto 100
      qf=2./3.
      if(iflc.eq.2.or.iflc.eq.3.or.iflc.eq.5.or.iflc.eq.7) qf=-1./3.
      vq=1.
      if(mod(mste(3),2).eq.1) vq=sqrt(max(0.,1.-(2.*pmq/ecmc)**2))

c...calculate weight in qed or qfd case
      if(mste(2).le.1) then
        rf=qf**2
        rfv=0.5*vq*(3.-vq**2)*qf**2
      else
        vf=sign(1.,qf)-4.*qf*pare(5)
        rf=qf**2*poll-2.*qf*vf*hf1i+(vf**2+1.)*hf1w
        rfv=0.5*vq*(3.-vq**2)*(qf**2*poll-2.*qf*vf*hf1i+vf**2*hf1w)+
     &  vq**3*hf1w
      endif

c...weighting or new event (radiative photon), cross section update
      if(ifl.le.0.and.rf.lt.rluC(0)*rfmax) goto 100
      pare(58)=pare(58)+1.
      if(ecmc.lt.2.*pmq+pare(10).or.rfv.lt.rluC(0)*rf) iflc=0
      if(mste(7).le.0.and.iflc.eq.0) goto 100
      if(iflc.ne.0) pare(59)=pare(59)+1.
      pare(44)=pare(57)*pare(59)/pare(58)
      pare(48)=pare(44)*86.8/ecm**2

      return
      end

c*********************************************************************

      subroutine luxjetC(ecm,njet,cut)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)

c...alpha-strong, total cross section, initial value for cut
      alspi=ulalpsC(ecm**2)/par(71)
      rqcd=1.+alspi
      if(iabs(mste(1)).ge.2) rqcd=1.+alspi+pare(65)*alspi**2
      cut=max(0.001,(pare(9)/ecm)**2,pare(8),exp(-sqrt(0.75/alspi))/2.)
      if(iabs(mste(1)).ge.2) cut=max(cut,0.2*alspi)

c...parametrization of first and second order three-jet cross section
  100 if(mste(1).eq.0.or.mste(1).ge.3.or.cut.ge.0.25) then
        pare(52)=0.
      else
        pare(52)=(2.*alspi/3.)*((3.-6.*cut+2.*alog(cut))*alog(cut/(1.-
     &  2.*cut))+(2.5+1.5*cut-6.571)*(1.-3.*cut)+5.833*(1.-3.*cut)**2-
     &  3.894*(1.-3.*cut)**3+1.342*(1.-3.*cut)**4)/rqcd
      endif
      if(iabs(mste(1)).le.1.or.mste(1).ge.3.or.cut.ge.0.25) then
        pare(53)=0.
      else
        ct=alog(1./cut-2.)
        pare(53)=alspi**2*ct**2*(2.419+0.5989*ct+0.6782*ct**2-
     &  0.2661*ct**3+0.01159*ct**4)/rqcd
      endif

c...parametrization of second order four-jet cross section
      if(iabs(mste(1)).le.1.or.mste(1).ge.3.or.cut.ge.0.125) then
        pare(54)=0.
      else
        ct=alog(1./cut-5.)
        if(cut.le.0.018) then
          xqqgg=6.349-4.330*ct+0.8304*ct**2
          xqqqq=-0.1080+0.01486*ct+0.009364*ct**2
        else
          xqqgg=-0.09773+0.2959*ct-0.2764*ct**2+0.08832*ct**3
          xqqqq=0.003661-0.004888*ct-0.001081*ct**2+0.002093*ct**3
        endif
        pare(54)=alspi**2*ct**2*(xqqgg+xqqqq)/rqcd
        pare(55)=xqqqq/(xqqgg+xqqqq)
      endif

c...if too high cross section harder cuts, else choose jet number
      if(pare(52)+pare(53)+pare(54).ge.1.) then
        cut=0.51*(2.*cut)**sqrt(1./(pare(52)+pare(53)+pare(54)))
        goto 100
      endif
      pare(50)=cut
      if(mste(1).le.0) then
        njet=min(4,2-mste(1))
      elseif(mste(1).ge.3) then
        njet=2
      else
        rnj=rluC(0)
        njet=2
        if(pare(52)+pare(53)+pare(54).gt.rnj) njet=3
        if(pare(54).gt.rnj) njet=4
      endif

      return
      end

c*********************************************************************

      subroutine lux3jtC(njet,cut,ifl,ecm,x1,x2)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)
c...dilogarithm of x for x<0.5 (x>0.5 obtained by analytic trick)
      dilog(x)=x+x**2/4.+x**3/9.+x**4/16.+x**5/25.+x**6/36.+x**7/49.

c...event type, mass effect factors and other common constants
      mste(35)=2
      pmq=ulmassC(2,ifl)
      qme=(2.*pmq/ecm)**2
      cutl=alog(cut)
      cutd=alog(1./cut-2.)
      wtmx=min(20.,37.-6.*cutd)
  100 njet=3

c...choose three-jet events in allowed region
  110 y13l=cutl+cutd*rluC(0)
      y23l=cutl+cutd*rluC(0)
      y13=exp(y13l)
      y23=exp(y23l)
      y12=1.-y13-y23
      if(y12.le.cut) goto 110
      if(y13**2+y23**2+2.*y12.le.2.*rluC(0)) goto 110

c...second order corrections
      if(mste(1).eq.2) then
        y12l=alog(y12)
        y13m=alog(1.-y13)
        y23m=alog(1.-y23)
        y12m=alog(1.-y12)
        if(y13.le.0.5) y13i=dilog(y13)
        if(y13.ge.0.5) y13i=1.644934-y13l*y13m-dilog(1.-y13)
        if(y23.le.0.5) y23i=dilog(y23)
        if(y23.ge.0.5) y23i=1.644934-y23l*y23m-dilog(1.-y23)
        if(y12.le.0.5) y12i=dilog(y12)
        if(y12.ge.0.5) y12i=1.644934-y12l*y12m-dilog(1.-y12)
        wt1=(y13**2+y23**2+2.*y12)/(y13*y23)
        wt2=(4./3.)*(-2.*(cutl-y12l)**2-3.*cutl-1.+3.289868+
     &  2.*(2.*cutl-y12l)*cut/y12)+3.*((cutl-y12l)**2-(cutl-y13l)**2-
     &  (cutl-y23l)**2-11.*cutl/6.+67./18.+1.644934-(2.*cutl-y12l)*
     &  cut/y12+(2.*cutl-y13l)*cut/y13+(2.*cutl-y23l)*cut/y23)+
     &  2.*(2.*cutl/3.-10./9.)+(4./3.)*(y12/(y12+y13)+
     &  y12/(y12+y23)+(y12+y23)/y13+(y12+y13)/y23+y13l*(4.*y12**2+
     &  2.*y12*y13+4.*y12*y23+y13*y23)/(y12+y23)**2+y23l*(4.*y12**2+
     &  2.*y12*y23+4.*y12*y13+y13*y23)/(y12+y13)**2)/wt1+
     &  3.*(y13l*y13/(y12+y23)+y23l*y23/(y12+y13))/wt1+
     &  (1./3.)*((y12**2+(y12+y13)**2)*(y12l*y23l-y12l*y12m-y23l*y23m+
     &  1.644934-y12i-y23i)/(y13*y23)+(y12**2+(y12+y23)**2)*
     &  (y12l*y13l-y12l*y12m-y13l*y13m+1.644934-y12i-y13i)/
     &  (y13*y23)+(y13**2+y23**2)/(y13*y23*(y13+y23))-
     &  2.*y12l*y12**2/(y13+y23)**2-4.*y12l*y12/(y13+y23))/wt1-
     &  3.*(y13l*y23l-y13l*y13m-y23l*y23m+1.644934-y13i-y23i)
        if(1.+pare(49)*wt2/par(72).le.(1.+pare(49)*wtmx/par(72))*
     &  rluC(0)) goto 110
        pare(56)=pare(49)*wt2/par(72)/(1.+pare(49)*wt2/par(72))
      endif

c...impose mass cuts (gives two jets), for fixed jet number new try
      x1=1.-y23
      x2=1.-y13
      x3=1.-y12
      if(4.*y23*y13*y12/x3**2.le.qme) njet=2
      if(mod(mste(3),4).ge.2.and.iabs(mste(1)).le.1.and.qme*x3+
     &0.5*qme**2+(0.5*qme+0.25*qme**2)*((1.-x2)/(1.-x1)+
     &(1.-x1)/(1.-x2)).gt.(x1**2+x2**2)*rluC(0)) njet=2
      if(mste(1).eq.-1.and.njet.eq.2) goto 100

      return
      end

c*********************************************************************

      subroutine lux4jtC(njet,cut,ifl,ecm,ifln,x1,x2,x4,x12,x14)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)
      dimension wta(4),wtb(4),wtc(4),wtd(4),wte(4)

c...common constants, choice of process (qqgg or qqqq)
      pmq=ulmassC(2,ifl)
      qme=(2.*pmq/ecm)**2
      ct=alog(1./cut-5.)
  100 njet=4
      it=1
      if(pare(55).gt.rluC(0)) it=2
      if(mste(1).le.-3) it=-mste(1)-2
      if(it.eq.1) wtmx=0.7/cut**2
      if(it.eq.2) wtmx=0.3/cut**2
      id=1

c...sample the five kinematical variables (for qqgg preweighted in y34)
  110 y134=3.*cut+(1.-6.*cut)*rluC(0)
      y234=3.*cut+(1.-6.*cut)*rluC(0)
      if(it.eq.1) y34=(1.-5.*cut)*exp(-ct*rluC(0))
      if(it.eq.2) y34=cut+(1.-6.*cut)*rluC(0)
      if(y34.le.y134+y234-1..or.y34.ge.y134*y234) goto 110
      vt=rluC(0)
      cp=cos(par(71)*rluC(0))
      y14=(y134-y34)*vt
      y13=y134-y14-y34
      vb=y34*(1.-y134-y234+y34)/((y134-y34)*(y234-y34))
      y24=0.5*(y234-y34)*(1.-4.*sqrt(max(0.,vt*(1.-vt)*vb*(1.-vb)))*
     &cp-(1.-2.*vt)*(1.-2.*vb))
      y23=y234-y34-y24
      y12=1.-y134-y23-y24
      if(min(y12,y13,y14,y23,y24).le.cut) goto 110
      y123=y12+y13+y23
      y124=y12+y14+y24

c...calculate matrix element for qqgg or qqqq process
      ic=0
      wttot=0.
  120 ic=ic+1
      if(it.eq.1) then
        wta(ic)=(y12*y34**2-y13*y24*y34+y14*y23*y34+3.*y12*y23*y34+
     &  3.*y12*y14*y34+4.*y12**2*y34-y13*y23*y24+2.*y12*y23*y24-
     &  y13*y14*y24-2.*y12*y13*y24+2.*y12**2*y24+y14*y23**2+2.*y12*
     &  y23**2+y14**2*y23+4.*y12*y14*y23+4.*y12**2*y23+2.*y12*y14**2+
     &  2.*y12*y13*y14+4.*y12**2*y14+2.*y12**2*y13+2.*y12**3)/(2.*y13*
     &  y134*y234*y24)+(y24*y34+y12*y34+y13*y24-y14*y23+y12*y13)/(y13*
     &  y134**2)+2.*y23*(1.-y13)/(y13*y134*y24)+y34/(2.*y13*y24)
        wtb(ic)=(y12*y24*y34+y12*y14*y34-y13*y24**2+y13*y14*y24+2.*y12*
     &  y14*y24)/(y13*y134*y23*y14)+y12*(1.+y34)*y124/(y134*y234*y14*
     &  y24)-(2.*y13*y24+y14**2+y13*y23+2.*y12*y13)/(y13*y134*y14)+
     &  y12*y123*y124/(2.*y13*y14*y23*y24)
        wtc(ic)=-(5.*y12*y34**2+2.*y12*y24*y34+2.*y12*y23*y34+2.*y12*
     &  y14*y34+2.*y12*y13*y34+4.*y12**2*y34-y13*y24**2+y14*y23*y24+
     &  y13*y23*y24+y13*y14*y24-y12*y14*y24-y13**2*y24-3.*y12*y13*y24-
     &  y14*y23**2-y14**2*y23+y13*y14*y23-3.*y12*y14*y23-y12*y13*y23)/
     &  (4.*y134*y234*y34**2)+(3.*y12*y34**2-3.*y13*y24*y34+3.*y12*y24*
     &  y34+3.*y14*y23*y34-y13*y24**2-y12*y23*y34+6.*y12*y14*y34+2.*y12*
     &  y13*y34-2.*y12**2*y34+y14*y23*y24-3.*y13*y23*y24-2.*y13*y14*
     &  y24+4.*y12*y14*y24+2.*y12*y13*y24+3.*y14*y23**2+2.*y14**2*y23+
     &  2.*y14**2*y12+2.*y12**2*y14+6.*y12*y14*y23-2.*y12*y13**2-
     &  2.*y12**2*y13)/(4.*y13*y134*y234*y34)
        wtc(ic)=wtc(ic)+(2.*y12*y34**2-2.*y13*y24*y34+y12*y24*y34+
     &  4.*y13*y23*y34+4.*y12*y14*y34+2.*y12*y13*y34+2.*y12**2*y34-
     &  y13*y24**2+3.*y14*y23*y24+4.*y13*y23*y24-2.*y13*y14*y24+
     &  4.*y12*y14*y24+2.*y12*y13*y24+2.*y14*y23**2+4.*y13*y23**2+
     &  2.*y13*y14*y23+2.*y12*y14*y23+4.*y12*y13*y23+2.*y12*y14**2+4.*
     &  y12**2*y13+4.*y12*y13*y14+2.*y12**2*y14)/(4.*y13*y134*y24*y34)-
     &  (y12*y34**2-2.*y14*y24*y34-2.*y13*y24*y34-y14*y23*y34+y13*y23*
     &  y34+y12*y14*y34+2.*y12*y13*y34-2.*y14**2*y24-4.*y13*y14*y24-
     &  4.*y13**2*y24-y14**2*y23-y13**2*y23+y12*y13*y14-y12*y13**2)/
     &  (2.*y13*y34*y134**2)+(y12*y34**2-4.*y14*y24*y34-2.*y13*y24*y34-
     &  2.*y14*y23*y34-4.*y13*y23*y34-4.*y12*y14*y34-4.*y12*y13*y34-
     &  2.*y13*y14*y24+2.*y13**2*y24+2.*y14**2*y23-2.*y13*y14*y23-
     &  y12*y14**2-6.*y12*y13*y14-y12*y13**2)/(4.*y34**2*y134**2)
        wttot=wttot+y34*(4.*wta(ic)-0.5*wtb(ic)+9.*wtc(ic))/18.
      else
        wtd(ic)=(y13*y23*y34+y12*y23*y34-y12**2*y34+y13*y23*y24+2.*y12*
     &  y23*y24-y14*y23**2+y12*y13*y24+y12*y14*y23+y12*y13*y14)/(y13**2*
     &  y123**2)-(y12*y34**2-y13*y24*y34+y12*y24*y34-y14*y23*y34-y12*
     &  y23*y34-y13*y24**2+y14*y23*y24-y13*y23*y24-y13**2*y24+y14*
     &  y23**2)/(y13**2*y123*y134)+(y13*y14*y12+y34*y14*y12-y34**2*y12+
     &  y13*y14*y24+2.*y34*y14*y24-y23*y14**2+y34*y13*y24+y34*y23*y14+
     &  y34*y13*y23)/(y13**2*y134**2)-(y34*y12**2-y13*y24*y12+y34*y24*
     &  y12-y23*y14*y12-y34*y14*y12-y13*y24**2+y23*y14*y24-y13*y14*y24-
     &  y13**2*y24+y23*y14**2)/(y13**2*y134*y123)
        wte(ic)=(y12*y34*(y23-y24+y14+y13)+y13*y24**2-y14*y23*y24+y13*
     &  y23*y24+y13*y14*y24+y13**2*y24-y14*y23*(y14+y23+y13))/(y13*y23*
     &  y123*y134)-y12*(y12*y34-y23*y24-y13*y24-y14*y23-y14*y13)/(y13*
     &  y23*y123**2)-(y14+y13)*(y24+y23)*y34/(y13*y23*y134*y234)+
     &  (y12*y34*(y14-y24+y23+y13)+y13*y24**2-y23*y14*y24+y13*y14*y24+
     &  y13*y23*y24+y13**2*y24-y23*y14*(y14+y23+y13))/(y13*y14*y134*
     &  y123)-y34*(y34*y12-y14*y24-y13*y24-y23*y14-y23*y13)/(y13*y14*
     &  y134**2)-(y23+y13)*(y24+y14)*y12/(y13*y14*y123*y124)
        wttot=wttot+(2.*wtd(ic)-wte(ic)/6.)/12.
      endif

c...permutations of momenta in matrix element, weighting
  130 if(ic.eq.1.or.ic.eq.3.or.id.eq.2.or.id.eq.3) then
        ysav=y13
        y13=y14
        y14=ysav
        ysav=y23
        y23=y24
        y24=ysav
        ysav=y123
        y123=y124
        y124=ysav
      endif
      if(ic.eq.2.or.ic.eq.4.or.id.eq.3.or.id.eq.4) then
        ysav=y13
        y13=y23
        y23=ysav
        ysav=y14
        y14=y24
        y24=ysav
        ysav=y134
        y134=y234
        y234=ysav
      endif
      if(ic.le.3) goto 120
      if(id.eq.1.and.wttot.lt.rluC(0)*wtmx) goto 110
      ic=5

      if(it.eq.1) then
c...qqgg events: string configuration, event type
        if(id.eq.1) then
          if(wta(2)+wta(4)+2.*(wtc(2)+wtc(4)).gt.rluC(0)*(wta(1)+wta(2)+
     &    wta(3)+wta(4)+2.*(wtc(1)+wtc(2)+wtc(3)+wtc(4)))) id=2
          if(id.eq.2) goto 130
        endif
        mste(35)=3
        if(0.5*y34*(wtc(1)+wtc(2)+wtc(3)+wtc(4)).gt.rluC(0)*wttot)
     &  mste(35)=4
        pare(56)=y34*(2.*(wta(1)+wta(2)+wta(3)+wta(4))+4.*(wtc(1)+
     &  wtc(2)+wtc(3)+wtc(4)))/(9.*wttot)
        ifln=0

c...mass cuts, kinematical variables out
        if(y12.le.cut+qme) njet=2
        if(njet.eq.2) goto 150
        q12=0.5*(1.-sqrt(1.-qme/y12))
        x1=1.-(1.-q12)*y234-q12*y134
        x4=1.-(1.-q12)*y134-q12*y234
        x2=1.-y124
        x12=(1.-q12)*y13+q12*y23
        x14=y12-0.5*qme
        if(y134*y234/((1.-x1)*(1.-x4)).le.rluC(0)) njet=2

      else
c...qqqq events: string configuration, choose new flavour
        if(id.eq.1) then
          wtr=rluC(0)*(wtd(1)+wtd(2)+wtd(3)+wtd(4))
          if(wtr.lt.wtd(2)+wtd(3)+wtd(4)) id=2
          if(wtr.lt.wtd(3)+wtd(4)) id=3
          if(wtr.lt.wtd(4)) id=4
          if(id.ge.2) goto 130
        endif
        mste(35)=5
        pare(56)=(wtd(1)+wtd(2)+wtd(3)+wtd(4))/(6.*wttot)
  140   ifln=1+int(4.*rluC(0))
        if(ifln.ne.ifl.and.0.25*pare(56).le.rluC(0)) goto 140
        if(ifln.eq.ifl.and.1.-0.75*pare(56).le.rluC(0)) goto 140
        if(ifln.gt.mste(4)) njet=2
        pmqn=ulmassC(2,ifln)
        qmen=(2.*pmqn/ecm)**2

c...mass cuts, kinematical variables out
        if(y24.le.cut+qme.or.y13.le.1.1*qmen) njet=2
        if(njet.eq.2) goto 150
        q24=0.5*(1.-sqrt(1.-qme/y24))
        q13=0.5*(1.-sqrt(1.-qmen/y13))
        x1=1.-(1.-q24)*y123-q24*y134
        x4=1.-(1.-q24)*y134-q24*y123
        x2=1.-(1.-q13)*y234-q13*y124
        x12=(1.-q24)*((1.-q13)*y14+q13*y34)+q24*((1.-q13)*y12+q13*y23)
        x14=y24-0.5*qme
        x34=(1.-q24)*((1.-q13)*y23+q13*y12)+q24*((1.-q13)*y34+q13*y14)
        if(pmq**2+pmqn**2+min(x12,x34)*ecm**2.le.(pare(10)+pmq+pmqn)**2)
     &  njet=2
        if(y123*y134/((1.-x1)*(1.-x4)).le.rluC(0)) njet=2
      endif
  150 if(mste(1).le.-2.and.njet.eq.2) goto 100

      return
      end

c*********************************************************************

      subroutine luxdifC(nc,njet,ifl,ecm,chi,the,phi)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)

c...charge, factors depending on polarization for qed case
      qf=2./3.
      if(ifl.eq.2.or.ifl.eq.3.or.ifl.eq.5.or.ifl.eq.7) qf=-1./3.
      poll=1.-pare(11)*pare(12)
      pold=pare(12)-pare(11)
      if(mste(2).le.1) then
        hf1=poll
        hf2=0.
        hf3=pare(13)**2
        hf4=0.

c...factors depending on flavour, energy and polarization for qfd case
      else
        sff=1./(16.*pare(5)*(1.-pare(5)))
        sfw=ecm**4/((ecm**2-pare(6)**2)**2+(pare(6)*pare(7))**2)
        sfi=sfw*(1.-(pare(6)/ecm)**2)
        ae=-1.
        ve=4.*pare(5)-1.
        af=sign(1.,qf)
        vf=af-4.*qf*pare(5)
        hf1=qf**2*poll-2.*qf*vf*sfi*sff*(ve*poll-ae*pold)+
     &  (vf**2+af**2)*sfw*sff**2*((ve**2+ae**2)*poll-2.*ve*ae*pold)
        hf2=-2.*qf*af*sfi*sff*(ae*poll-ve*pold)+2.*vf*af*sfw*sff**2*
     &  (2.*ve*ae*poll-(ve**2+ae**2)*pold)
        hf3=pare(13)**2*(qf**2-2.*qf*vf*sfi*sff*ve+(vf**2+af**2)*
     &  sfw*sff**2*(ve**2-ae**2))
        hf4=-pare(13)**2*2.*qf*vf*sfw*(pare(6)*pare(7)/ecm**2)*sff*ae
      endif

c...mass factor, differential cross sections for two-jet events
      sq2=sqrt(2.)
      qme=0.
      if(mste(3).ge.4.and.iabs(mste(1)).le.1.and.mste(2).le.1) qme=
     &(2.*ulmassC(2,ifl)/ecm)**2
      if(njet.eq.2) then
        sigu=4.*sqrt(1.-qme)
        sigl=2.*qme*sqrt(1.-qme)
        sigt=0.
        sigi=0.
        siga=0.
        sigp=4.

c...kinematical variables, reduce four-jet event to three-jet one
      else
        if(njet.eq.3) then
          x1=2.*p(nc+1,4)/ecm
          x2=2.*p(nc+3,4)/ecm
        else
          ecmr=p(nc+1,4)+p(nc+4,4)+sqrt((p(nc+2,1)+p(nc+3,1))**2+
     &    (p(nc+2,2)+p(nc+3,2))**2+(p(nc+2,3)+p(nc+3,3))**2)
          x1=2.*p(nc+1,4)/ecmr
          x2=2.*p(nc+4,4)/ecmr
        endif

c...differential cross sections for three-jets (or reduced four-jets)
        xq=(1.-x1)/(1.-x2)
        ct12=(x1*x2-2.*x1-2.*x2+2.+qme)/sqrt((x1**2-qme)*(x2**2-qme))
        st12=sqrt(1.-ct12**2)
        sigu=2.*x1**2+x2**2*(1.+ct12**2)-qme*(3.+ct12**2-x1-x2)-
     &  qme*x1/xq+0.5*qme*((x2**2-qme)*st12**2-2.*x2)*xq
        sigl=(x2*st12)**2-qme*(3.-ct12**2-2.5*(x1+x2)+x1*x2+qme)+
     &  0.5*qme*(x1**2-x1-qme)/xq+0.5*qme*((x2**2-qme)*ct12**2-x2)*xq
        sigt=0.5*(x2**2-qme-0.5*qme*(x2**2-qme)/xq)*st12**2
        sigi=((1.-0.5*qme*xq)*(x2**2-qme)*st12*ct12+qme*(1.-x1-x2+
     &  0.5*x1*x2+0.5*qme)*st12/ct12)/sq2
        siga=x2**2*st12/sq2
        sigp=2.*(x1**2-x2**2*ct12)
      endif

c...upper bound for differential cross section
      hf1a=abs(hf1)
      hf2a=abs(hf2)
      hf3a=abs(hf3)
      hf4a=abs(hf4)
      sigmax=(2.*hf1a+hf3a+hf4a)*abs(sigu)+2.*(hf1a+hf3a+hf4a)*abs(
     &sigl)+2.*(hf1a+2.*hf3a+2.*hf4a)*abs(sigt)+2.*sq2*(hf1a+2.*hf3a+
     &2.*hf4a)*abs(sigi)+4.*sq2*hf2a*abs(siga)+2.*hf2a*abs(sigp)

c...generate angular orientation according to differential cross section
  100 chi=par(72)*rluC(0)
      cthe=2.*rluC(0)-1.
      phi=par(72)*rluC(0)
      cchi=cos(chi)
      schi=sin(chi)
      c2chi=cos(2.*chi)
      s2chi=sin(2.*chi)
      the=acos(cthe)
      sthe=sin(the)
      c2phi=cos(2.*(phi-pare(14)))
      s2phi=sin(2.*(phi-pare(14)))
      sig=((1.+cthe**2)*hf1+sthe**2*(c2phi*hf3-s2phi*hf4))*sigu+
     &2.*(sthe**2*hf1-sthe**2*(c2phi*hf3-s2phi*hf4))*sigl+
     &2.*(sthe**2*c2chi*hf1+((1.+cthe**2)*c2chi*c2phi-2.*cthe*s2chi*
     &s2phi)*hf3-((1.+cthe**2)*c2chi*s2phi+2.*cthe*s2chi*c2phi)*hf4)*
     &sigt-2.*sq2*(2.*sthe*cthe*cchi*hf1-2.*sthe*(cthe*cchi*c2phi-
     &schi*s2phi)*hf3+2.*sthe*(cthe*cchi*s2phi+schi*c2phi)*hf4)*sigi+
     &4.*sq2*sthe*cchi*hf2*siga+2.*cthe*hf2*sigp
      if(sig.lt.sigmax*rluC(0)) goto 100

      return
      end

c*********************************************************************

      subroutine luoniaC(ifl,ecm)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)

c...printout or resetting of statistics
      if(ifl.ge.10) call luestaC(ifl,0,0,0.,0,0.,0)
      if(ifl.ge.10) return
      if(mst(19).ge.1) call lulistC(-1)

c...initial e+e- and onium state (optional)
      nc=0
      if(mste(30).ge.2) then
        nc=nc+2
        call lupartC(nc-1,7,0.5*ecm,0.,0.)
        k(nc-1,1)=40000
        call lupartC(nc,-7,0.5*ecm,par(71),0.)
        k(nc,1)=40000
      endif
      iflc=iabs(ifl)
      if(mste(30).ge.3.and.iflc.ge.5) then
        nc=nc+1
        kf=82+iabs(ifl)
        mst(9)=1
        p(nc,5)=ecm
        call lupartC(nc,kf,ecm,0.,0.)
        k(nc,1)=50001
        mst(9)=0
      endif

c...choose x1 and x2 according to matrix element
      ntry=0
  100 x1=rluC(0)
      x2=rluC(0)
      x3=2.-x1-x2
      if(x3.ge.1..or.((1.-x1)/(x2*x3))**2+((1.-x2)/(x1*x3))**2+
     &((1.-x3)/(x1*x2))**2.le.2.*rluC(0)) goto 100
      ntry=ntry+1
      njet=3
      call lu3jetC(nc+1,0,0,ecm,x1,x3)

c...photon-gluon-gluon events, small system modifications, jet origin
      qf=0.
      if(iflc.eq.1.or.iflc.eq.4.or.iflc.eq.6.or.iflc.eq.8) qf=2./3.
      if(iflc.eq.2.or.iflc.eq.3.or.iflc.eq.5.or.iflc.eq.7) qf=-1./3.
      rgam=7.2*qf**2*pare(4)/ulalpsC(ecm**2)
      mk=0
      ecmc=ecm
      if(rluC(0).gt.rgam/(1.+rgam)) then
        if(1.-max(x1,x2,x3).le.max((pare(9)/ecm)**2,pare(8))) njet=2
        if(njet.eq.2) call lu2jetC(nc+1,0,0,ecm)
      else
        mk=1
        ecmc=sqrt(1.-x1)*ecm
        if(ecmc.lt.2.*pare(10)) goto 100
        k(nc+1,1)=0
        k(nc+1,2)=1
        njet=2
        if(ecmc.lt.4.*pare(10)) njet=0
        mst(9)=1
        if(njet.eq.0) p(nc+2,5)=ecmc
        if(njet.eq.0) call lupartC(nc+2,15,0.5*(x2+x3)*ecm,par(71),1.)
        mst(9)=0
      endif
      do 110 ip=nc+1,n
  110 k(ip,1)=k(ip,1)+(mste(30)/2)+(iflc/5)*(mste(30)/3)*(nc-1)

c...differential cross sections, upper limit for cross section
      if(mste(6).eq.1) then
        sq2=sqrt(2.)
        hf1=1.-pare(11)*pare(12)
        hf3=pare(13)**2
        ct13=(x1*x3-2.*x1-2.*x3+2.)/(x1*x3)
        st13=sqrt(1.-ct13**2)
        sigl=0.5*x3**2*((1.-x2)**2+(1.-x3)**2)*st13**2
        sigu=(x1*(1.-x1))**2+(x2*(1.-x2))**2+(x3*(1.-x3))**2-sigl
        sigt=0.5*sigl
        sigi=(sigl*ct13/st13+0.5*x1*x3*(1.-x2)**2*st13)/sq2
        sigmax=(2.*hf1+hf3)*abs(sigu)+2.*(hf1+hf3)*abs(sigl)+2.*(hf1+
     &  2.*hf3)*abs(sigt)+2.*sq2*(hf1+2.*hf3)*abs(sigi)

c...angular orientation of event
  120   chi=par(72)*rluC(0)
        cthe=2.*rluC(0)-1.
        phi=par(72)*rluC(0)
        cchi=cos(chi)
        schi=sin(chi)
        c2chi=cos(2.*chi)
        s2chi=sin(2.*chi)
        the=acos(cthe)
        sthe=sin(the)
        c2phi=cos(2.*(phi-pare(14)))
        s2phi=sin(2.*(phi-pare(14)))
        sig=((1.+cthe**2)*hf1+sthe**2*c2phi*hf3)*sigu+2.*(sthe**2*hf1-
     &  sthe**2*c2phi*hf3)*sigl+2.*(sthe**2*c2chi*hf1+((1.+cthe**2)*
     &  c2chi*c2phi-2.*cthe*s2chi*s2phi)*hf3)*sigt-2.*sq2*(2.*sthe*cthe*
     &  cchi*hf1-2.*sthe*(cthe*cchi*c2phi-schi*s2phi)*hf3)*sigi
        if(sig.lt.sigmax*rluC(0)) goto 120
        mst(1)=nc+1
        call luroboC(0.,chi,0.,0.,0.)
        call luroboC(the,phi,0.,0.,0.)
        mst(1)=0
      endif

      if(mste(1).ge.3.and.njet.ge.2) then
c...prepare event record for parton shower evolution
        do 130 i=njet,2,-1
        k(nc+mk+2*i-1,1)=k(nc+mk+i,1)
        k(nc+mk+2*i-1,2)=k(nc+mk+i,2)
        do 130 j=1,5
  130   p(nc+mk+2*i-1,j)=p(nc+mk+i,j)
        do 140 i=nc+mk+2,nc+mk+2*njet,2
        k(i,1)=70000+i-1
        k(i,2)=1000+mod(k(i-1,1),1000)
        p(i,1)=i-3
        p(i,2)=i+1
        p(i,3)=0.
        p(i,4)=0.
  140   p(i,5)=0.
        p(nc+mk+2,1)=nc+mk+2*njet-1
        p(nc+mk+2*njet,2)=nc+mk+1
        n=n+njet

c...generate parton shower, rearrange along strings and check
        call lushowC(nc+mk+1,-njet,ecmc)
        njet=0
        do 150 i=1,n
  150   if(k(i,1).lt.20000.and.iabs(k(i,2)).ge.500) njet=njet+1
        mst12s=mst(12)
        if(mste(5).eq.-1) mst(12)=0
        mst21s=mst(21)
        if(mste(30).le.3) mst(21)=1
        if(mste(30).ge.4) mst(21)=0
        if(mste(5).ge.0) mst(26)=0
        call luprepC
        mst(12)=mst12s
        mst(21)=mst21s
        if(mste(5).ge.0.and.mst(26).ne.0) goto 100
        nc=n+1
  160   nc=nc-1
        if(iabs(k(nc,2)).lt.500.and.mod(k(nc,1),10000).gt.0) then
          if(k(mod(k(nc,1),10000),1)/20000.eq.1) goto 160
        endif
      else
        njet=3-mk
      endif

c...event generation and statistics
      nc=n
      if(mste(5).eq.1) call luexecC
      if(mste(31).ne.0) call luestaC(9,3-mk,nc,ecm,mk,ecmc,ntry)

      return
      end

c***********************************************************************

      subroutine luestaC(ifl,njet,nc,ecm,mk,ecmc,ntry)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludat3C/dpar(20),idb(120),cbr(400),kdp(1600)
      common/ludat4C/chag(50),chaf(100)
      common/ludateC/mste(40),pare(80)
      character chag*4,chaf*4,chap*8,chad*4
      dimension lfl(10,10),lkf(400,4),ecms(2)
      data lfl/100*0/,lkf/1600*0/,ecms/2*0./,ntrys/0/

c...fill statistics on initial state of event
      if(ifl.lt.10) then
        ecms(1)=ecms(1)+ecm
        ecms(2)=ecms(2)+ecmc
        njt=njet
        if(njet.lt.0) njt=5
        if(njet.ge.5) njt=min(njet+1,9)
        lfl(1,1)=lfl(1,1)+1
        lfl(1,njt)=lfl(1,njt)+1
        if(mk.eq.1) lfl(1,10)=lfl(1,10)+1
        if(ifl.ge.1) lfl(ifl+1,1)=lfl(ifl+1,1)+1
        if(ifl.ge.1) lfl(ifl+1,njt)=lfl(ifl+1,njt)+1
        if(ifl.ge.1.and.mk.eq.1) lfl(ifl+1,10)=lfl(ifl+1,10)+1
        ntrys=ntrys+ntry

c...fill statistics on final state of event
        if(iabs(mste(31)).le.1) return
        do 100 i=1,n
        kfa=iabs(k(i,2))
        if(kfa.ge.400.or.k(i,1).ge.40000) goto 100
        if(mod(k(i,1),10000).le.nc) lkf(400,1)=lkf(400,1)+1
        if(k(i,1).lt.20000) lkf(400,2)=lkf(400,2)+1
        if(k(i,1).lt.20000.and.luchgeC(kfa).ne.0)
     &      lkf(400,3)=lkf(400,3)+1
        kfs=2-isign(1,k(i,2))
        if(mod(k(i,1),10000).gt.nc) kfs=kfs+1
        lkf(kfa,kfs)=lkf(kfa,kfs)+1
  100   continue

c...write e+e- parameter table and initial state statistics
      elseif(ifl.eq.10) then
        if(mste(31).ge.0) write(mst(20),1000) (l,mste(l),mste(l+20),
     &  pare(l),pare(l+20),pare(l+40),pare(l+60),l=1,20)
        if(mste(31).eq.0.or.lfl(1,1).eq.0) return
        fac=1./float(lfl(1,1))
        if(mste(31).eq.-2) goto 120
        write(mst(20),1100) (fac*ecms(l),l=1,2),(lfl(1,j),j=1,10),
     &  (fac*lfl(1,j),j=2,10)
        do 110 l=2,10
  110   if(lfl(l,1).ne.0) write(mst(20),1200) l-1,chag(9+l-9*(l/10)),
     &  (lfl(l,j),j=1,10),fac*lfl(l,1),(lfl(l,j)/float(lfl(l,1)),j=2,10)
        write(mst(20),1300) (ntrys-lfl(1,1))/float(ntrys)

c...write final state statistics
  120   if(iabs(mste(31)).le.1) return
        write(mst(20),1400) (fac*lkf(400,l),l=1,3)
        do 130 l=1,399
        lkfs=lkf(l,1)+lkf(l,2)+lkf(l,3)+lkf(l,4)
        if(lkfs.eq.0) goto 130
        call lunameC(l,chap)
        kfa=l
        if(kfa.gt.100) call luiflvC(kfa,ifla,iflb,iflc,ksp)
        if(kfa.gt.100) kfa=76+5*ifla+ksp
        chad=chag(41)
        if(idb(kfa).ge.1) chad=chag(42)
        pm=ulmassC(0,l)
        write(mst(20),1500) l,chap,chad,pm,(lkf(l,j),j=1,4),fac*lkfs
  130   continue

c...reset counters
      else
        do 140 l=1,10
        do 140 j=1,7
  140   lfl(l,j)=0
        do 150 l=1,400
        do 150 j=1,4
  150   lkf(l,j)=0
        ecms(1)=0.
        ecms(2)=0.
        ntrys=0
      endif

c...format statements for output on unit mst(20) (default is 6)
 1000 format(///20x,'e+e- parameter value table'//5x,'l mste(l) ',
     &'&(l+20)      pare(l)      &(l+20)      &(l+40)      &(l+60)'/
     &20(/1x,i5,2(1x,i7),4(1x,f12.4)))
 1100 format(//20x,'event statistics - initial state'//5x,'mean total',
     &' energy =',f10.3/5x,'mean hadronic cm energy =',f10.3//5x,
     &'number and fraction of event types'/5x,'flavour',7x,'all',7x,
     &'qq      qqg     qqgg     qqqq     5jet     6jet     7jet',6x,
     &'>=8      rad'//6x,'all',4x,10(1x,i8)/22x,9(1x,f8.3))
 1200 format(/1x,i5,3x,a4,10(1x,i8)/13x,10(1x,f8.3))
 1300 format(/5x,'fraction of original events failing cuts =',f10.4)
 1400 format(//20x,'event statistics - final state'//5x,'mean primary',
     &' multiplicity =',f8.3/5x,'mean final   multiplicity =',f8.3/
     &5x,'mean charged multiplicity =',f8.3//5x,'number and fraction ',
     &'of particles produced (directly and via decays)'/41x,'part',
     &'icles       antiparticles    total/'/5x,'kf    particle/decay',
     &5x,'mass    #prim    #seco    #prim    #seco    /event'/)
 1500 format(1x,i6,4x,a8,2x,a4,1x,f8.3,4(1x,i8),1x,f9.4)

      return
      end

c*********************************************************************

      subroutine lushowC(ip1,ip2,qmax)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)
      dimension pmth(3,0:8),ps(5),pma(4),pmsd(4),iep(4),
     &ipa(4),ifla(4),ifld(4),ifl(4),itry(4),isi(4),isl(4)

c...initialization of cutoff masses etc.
      if(mste(11).le.0.or.qmax.le.pare(22)) return
      do 100 if=0,8
      pmth(1,if)=ulmassC(2,if)
      pmth(2,if)=sqrt(pmth(1,if)**2+0.25*pare(22)**2)
  100 pmth(3,if)=pmth(2,0)+pmth(2,if)
      pt2min=max(0.5*pare(22),1.1*pare(21))**2
      alams=pare(21)**2
      alfm=alog(pt2min/alams)

c...check on phase space available for emission: pair or single parton
      m3jc=0
      if(ip1.gt.0.and.ip2.eq.0) then
        npa=1
        ipa(1)=ip1
      elseif(ip1.gt.0.and.ip2.gt.0) then
        npa=2
        ipa(1)=ip1
        ipa(2)=ip2
      elseif(ip1.gt.0.and.ip2.lt.0) then
        npa=iabs(ip2)
        do 110 i=1,npa
  110   ipa(i)=ip1+2*(i-1)
      else
        return
      endif
      irej=0
      do 120 j=1,5
  120 ps(j)=0.
      pm=0.
      do 130 i=1,npa
      ifla(i)=iabs(k(ipa(i),2))-500
      pma(i)=p(ipa(i),5)
      if(ifla(i).ge.0.and.ifla(i).le.8) pma(i)=pmth(3,ifla(i))
      pm=pm+pma(i)
      if(ifla(i).lt.0.or.ifla(i).gt.8.or.pma(i).gt.qmax) irej=
     &irej+1
      do 130 j=1,4
  130 ps(j)=ps(j)+p(ipa(i),j)
      if(irej.eq.npa) return
      ps(5)=sqrt(max(0.,ps(4)**2-ps(1)**2-ps(2)**2-ps(3)**2))
      if(npa.eq.1) ps(5)=ps(4)
      if(ps(5).le.pm+pare(22)) return
      if(npa.eq.2.and.mste(13).ge.2) then
        if(ifla(1).ge.1.and.ifla(1).le.8.and.ifla(2).ge.1.and.
     &  ifla(2).le.8) m3jc=1
      endif

c...define imagined single initiator of shower for parton system
      ns=n
      if(npa.ge.2) then
        k(n+1,1)=0
        k(n+1,2)=500
        p(n+1,1)=0.
        p(n+1,2)=0.
        p(n+1,3)=0.
        p(n+1,4)=ps(5)
        p(n+1,5)=ps(5)
        p(n+2,5)=ps(5)**2
        n=n+2
      endif

c...loop over partons that may branch
      nep=npa
      im=ns-1
      if(npa.eq.1) im=ns-3
  140 im=im+2
      if(n.gt.ns) then
        if(im.gt.n) goto 330
        iflm=iabs(k(im,2))-500
        if(iflm.lt.0.or.iflm.gt.8) goto 140
        if(p(im,5).lt.pmth(2,iflm)) goto 140
        igm=mod(k(im,1),10000)
      else
        igm=-1
      endif

c...origin and flavour of daughters
      if(igm.ge.0) then
        k(im+1,1)=n+1
        do 150 i=1,nep
  150   k(n+2*i-1,1)=im
      else
        k(n+1,1)=ipa(1)
      endif
      if(igm.le.0) then
        do 160 i=1,nep
  160   k(n+2*i-1,2)=k(ipa(i),2)
      elseif(iflm.ne.0) then
        k(n+1,2)=k(im,2)
        k(n+3,2)=500
      elseif(k(im+1,2).eq.500) then
        k(n+1,2)=500
        k(n+3,2)=500
      else
        k(n+1,2)=k(im+1,2)
        k(n+3,2)=-k(im+1,2)
      endif
      do 170 ip=1,nep
      ifld(ip)=iabs(k(n+2*ip-1,2))-500
      itry(ip)=0
      isl(ip)=0
      isi(ip)=0
  170 if(ifld(ip).ge.0.and.ifld(ip).le.8) isi(ip)=1
      islm=0

c...maximum virtuality of daughters
      if(igm.le.0) then
        do 180 i=1,npa
        imp=n+2*i-1
        if(npa.ge.3) p(imp,4)=(ps(4)*p(ipa(i),4)-ps(1)*p(ipa(i),1)-
     &  ps(2)*p(ipa(i),2)-ps(3)*p(ipa(i),3))/ps(5)
        p(imp,5)=min(qmax,ps(5))
        if(npa.ge.3) p(imp,5)=min(p(imp,5),p(imp,4))
  180   if(isi(i).eq.0) p(imp,5)=p(ipa(i),5)
      else
        if(mste(12).le.2) pem=p(im+1,2)
        if(mste(12).ge.3) pem=p(im,4)
        p(n+1,5)=min(p(im,5),p(im+1,1)*pem)
        p(n+3,5)=min(p(im,5),(1.-p(im+1,1))*pem)
      endif
      do 190 i=1,nep
      pmsd(i)=p(n+2*i-1,5)
      if(isi(i).eq.1) then
        if(p(n+2*i-1,5).le.pmth(3,ifld(i))) p(n+2*i-1,5)=pmth(1,ifld(i))
      endif
  190 p(n+2*i,5)=p(n+2*i-1,5)**2

c...choose one of the daughters for evolution
  200 inum=0
      if(nep.eq.1) inum=1
      do 210 i=1,nep
  210 if(inum.eq.0.and.isl(i).eq.1) inum=i
      do 220 i=1,nep
      if(inum.eq.0.and.itry(i).eq.0.and.isi(i).eq.1) then
        if(p(n+2*i-1,5).ge.pmth(2,ifld(i))) inum=i
      endif
  220 continue
      if(inum.eq.0) then
        rmax=0.
        do 230 i=1,nep
        if(isi(i).eq.1.and.pmsd(i).ge.pmth(2,0)) then
          rpm=p(n+2*i-1,5)/pmsd(i)
          if(rpm.gt.rmax.and.p(n+2*i-1,5).ge.pmth(2,ifld(i))) then
            rmax=rpm
            inum=i
          endif
        endif
  230   continue
      endif
      inum=max(1,inum)
      iep(1)=n+2*inum-1
      do 240 i=2,nep
      iep(i)=iep(i-1)+2
  240 if(iep(i).ge.n+2*nep) iep(i)=n+1
      do 250 i=1,nep
  250 ifl(i)=iabs(k(iep(i),2))-500
      itry(inum)=itry(inum)+1
      z=0.5
      if(ifl(1).lt.0.or.ifl(1).gt.8) goto 290
      if(p(iep(1),5).lt.pmth(2,ifl(1))) goto 290

c...calculate allowed z range and integral of altarelli-parisi z kernel
      if(nep.eq.1) then
        pmed=ps(4)
      elseif(igm.eq.0.or.mste(12).le.2) then
        pmed=p(im,5)
      else
        if(inum.eq.1) pmed=p(im+1,1)*pem
        if(inum.eq.2) pmed=(1.-p(im+1,1))*pem
      endif
      if(mod(mste(12),2).eq.1) then
        zc=pmth(2,0)/pmed
      else
        zc=0.5*(1.-sqrt(max(0.,1.-(2.*pmth(2,0)/pmed)**2)))
        if(zc.lt.1e-4) zc=(pmth(2,0)/pmed)**2
      endif
      if(zc.gt.0.49) then
        p(iep(1),5)=pmth(1,ifl(1))
        p(iep(1)+1,5)=p(iep(1),5)**2
        goto 290
      endif
      if(ifl(1).eq.0) fbr=6.*alog((1.-zc)/zc)+mste(15)*(0.5-zc)
      if(ifl(1).ne.0) fbr=(8./3.)*alog((1.-zc)/zc)

c...inner veto algorithm starts, choose m and z
  260 pms=p(iep(1)+1,5)
      if(igm.ge.0) then
        pm2=0.
        do 270 i=2,nep
        pm=p(iep(i),5)
        if(ifl(i).ge.0.and.ifl(i).le.8) pm=pmth(2,ifl(i))
  270   pm2=pm2+pm
        pms=min(pms,(p(im,5)-pm2)**2)
      endif
      b0=27./6.
      do 280 if=4,mste(15)
  280 if(pms.gt.4.*pmth(2,if)**2) b0=(33.-2.*if)/6.
      if(mste(14).le.0) p(iep(1)+1,5)=pms*rluC(0)**(par(72)/pare(3)
     &/fbr)
      if(mste(14).eq.1) p(iep(1)+1,5)=4.*alams*(0.25*pms/alams)**
     &(rluC(0)**(b0/fbr))
      if(mste(14).ge.2) p(iep(1)+1,5)=pms*rluC(0)**(alfm*b0/fbr)
      p(iep(1),5)=sqrt(p(iep(1)+1,5))
      if(p(iep(1),5).le.pmth(3,ifl(1))) then
        p(iep(1),5)=pmth(1,ifl(1))
        p(iep(1)+1,5)=p(iep(1),5)**2
        goto 290
      endif
      if(ifl(1).ne.0) then
        z=1.-(1.-zc)*(zc/(1.-zc))**rluC(0)
        if(1.+z**2.lt.2.*rluC(0)) goto 260
      elseif(mste(15)*(0.5-zc).lt.rluC(0)*fbr) then
        z=(1.-zc)*(zc/(1.-zc))**rluC(0)
        if(rluC(0).gt.0.5) z=1.-z
        if((1.-z*(1.-z))**2.lt.rluC(0)) goto 260
        k(iep(1)+1,2)=500
      else
        z=zc+(1.-2.*zc)*rluC(0)
        if(z**2+(1.-z)**2.lt.rluC(0)) goto 260
        iflb=1+int(mste(15)*rluC(0))
        pmq=4.*pmth(2,iflb)**2/p(iep(1)+1,5)
        if(pmq.ge.1.) goto 260
        pmq0=4.*pmth(2,0)**2/p(iep(1)+1,5)
        if(mod(mste(12),2).eq.0.and.(1.+0.5*pmq)*sqrt(1.-pmq).lt.
     &  rluC(0)*(1.+0.5*pmq0)*sqrt(1.-pmq0)) goto 260
        k(iep(1)+1,2)=500+iflb
      endif
      if(mste(14).ge.2.and.z*(1.-z)*p(iep(1)+1,5).lt.pt2min) goto 260
      if(mste(14).ge.2.and.alfm/alog(p(iep(1)+1,5)*z*(1.-z)/alams).lt.
     &rluC(0)) goto 260

c...check if z consistent with chosen m
      if(ifl(1).eq.0) then
        iflgd1=iabs(k(iep(1)+1,2))-500
        iflgd2=iflgd1
      else
        iflgd1=ifl(1)
        iflgd2=0
      endif
      if(nep.eq.1) then
        ped=ps(4)
      elseif(nep.ge.3) then
        ped=p(iep(1),4)
      elseif(igm.eq.0.or.mste(12).le.2) then
        ped=0.5*(p(im+1,5)+p(iep(1)+1,5)-pm2**2)/p(im,5)
      else
        if(iep(1).eq.n+1) ped=p(im+1,1)*pem
        if(iep(1).eq.n+3) ped=(1.-p(im+1,1))*pem
      endif
      if(mod(mste(12),2).eq.1) then
        pmq1=pmth(2,iflgd1)**2/p(iep(1)+1,5)
        pmq2=pmth(2,iflgd2)**2/p(iep(1)+1,5)
        dz=sqrt(max(0.,(1.-p(iep(1)+1,5)/ped**2)*((1.-pmq1-pmq2)**2-
     &  4.*pmq1*pmq2)))
        zh=1.+pmq1-pmq2
      else
        dz=sqrt(max(0.,1.-p(iep(1)+1,5)/ped**2))
        zh=1.
      endif
      zl=0.5*(zh-dz)
      zu=0.5*(zh+dz)
      if(z.lt.zl.or.z.gt.zu) goto 260
      if(ifl(1).eq.0) p(iep(1)+1,3)=alog(zu*(1.-zl)/max(1e-20,zl*
     &(1.-zu)))
      if(ifl(1).ne.0) p(iep(1)+1,3)=alog((1.-zl)/max(1e-10,1.-zu))

c...three-jet matrix element correction or angular ordering
      if(igm.eq.0.and.m3jc.eq.1) then
        x1=1.-p(iep(1)+1,5)/p(ns+2,5)
        x2=z*(1.+p(iep(1)+1,5)/p(ns+2,5))
        x3=(1.-x1)+(1.-x2)
        wshow=1.+(1.-x1)/x3*(x1/(2.-x2))**2+(1.-x2)/x3*(x2/(2.-x1))**2
        if(x1**2+x2**2.lt.rluC(0)*wshow) goto 260
      elseif(igm.gt.0.and.mste(11).ge.2) then
        zm=p(im+1,1)
        if(iep(1).eq.n+3) zm=1.-p(im+1,1)
        if(z*(1.-z)/p(iep(1)+1,5).lt.(1.-zm)/(zm*p(im+1,5))) goto 260
      endif

c...end of inner veto algorithm, check if only one leg evolved so far
  290 p(iep(1)+1,1)=z
      isl(1)=0
      isl(2)=0
      if(nep.eq.1) goto 320
      if(nep.eq.2.and.p(iep(1),5)+p(iep(2),5).ge.p(im,5)) goto 200
      do 300 i=1,nep
      if(itry(i).eq.0.and.ifld(i).ge.0.and.ifld(i).le.8) then
        if(p(n+2*i-1,5).ge.pmth(2,ifld(i))) goto 200
      endif
  300 continue

c...check if chosen multiplet m1,m2,z1,z2 is physical
      if(nep.eq.4) then
      elseif(nep.eq.3) then
        pa1s=(p(n+1,4)+p(n+1,5))*(p(n+1,4)-p(n+1,5))
        pa2s=(p(n+3,4)+p(n+3,5))*(p(n+3,4)-p(n+3,5))
        pa3s=(p(n+5,4)+p(n+5,5))*(p(n+5,4)-p(n+5,5))
        pts=0.25*(2.*pa1s*pa2s+2.*pa1s*pa3s+2.*pa2s*pa3s-
     &  pa1s**2-pa2s**2-pa3s**2)/pa1s
        if(pts.le.0.) goto 200
      elseif(igm.eq.0.or.mste(12).le.2.or.mod(mste(12),2).eq.0) then
        do 310 i1=n+1,n+3,2
        iflda=iabs(k(i1,2))-500
        if(iflda.lt.0.or.iflda.gt.8) goto 310
        if(p(i1,5).lt.pmth(2,iflda)) goto 310
        if(iflda.eq.0) then
          iflgd1=iabs(k(i1+1,2))-500
          iflgd2=iflgd1
        else
          iflgd1=iflda
          iflgd2=0
        endif
        i2=2*n+4-i1
        if(igm.eq.0.or.mste(12).le.2) then
          ped=0.5*(p(im+1,5)+p(i1+1,5)-p(i2+1,5))/p(im,5)
        else
          if(i1.eq.n+1) zm=p(im+1,1)
          if(i1.eq.n+3) zm=1.-p(im+1,1)
          pml=sqrt((p(im+1,5)-p(n+2,5)-p(n+4,5))**2-
     &    4.*p(n+2,5)*p(n+4,5))
          ped=pem*(0.5*(p(im+1,5)-pml+p(i1+1,5)-p(i2+1,5))+
     &    pml*zm)/p(im+1,5)
        endif
        if(mod(mste(12),2).eq.1) then
          pmq1=pmth(2,iflgd1)**2/p(i1+1,5)
          pmq2=pmth(2,iflgd2)**2/p(i1+1,5)
          dz=sqrt((1.-p(i1+1,5)/ped**2)*((1.-pmq1-pmq2)**2-
     &    4.*pmq1*pmq2))
          zh=1.+pmq1-pmq2
        else
          dz=sqrt(1.-p(i1+1,5)/ped**2)
          zh=1.
        endif
        zl=0.5*(zh-dz)
        zu=0.5*(zh+dz)
        if(i1.eq.n+1.and.(p(i1+1,1).lt.zl.or.p(i1+1,1).gt.zu)) isl(1)=1
        if(i1.eq.n+3.and.(p(i1+1,1).lt.zl.or.p(i1+1,1).gt.zu)) isl(2)=1
        if(iflda.eq.0) p(i1+1,4)=alog(zu*(1.-zl)/max(1e-20,zl*(1.-zu)))
        if(iflda.ne.0) p(i1+1,4)=alog((1.-zl)/max(1e-10,1.-zu))
  310   continue
        if(isl(1).eq.1.and.isl(2).eq.1.and.islm.ne.0) then
          isl(3-islm)=0
          islm=3-islm
        elseif(isl(1).eq.1.and.isl(2).eq.1) then
          dzr1=max(0.,p(n+2,3)/p(n+2,4)-1.)
          dzr2=max(0.,p(n+4,3)/p(n+4,4)-1.)
          if(dzr2.gt.rluC(0)*(dzr1+dzr2)) isl(1)=0
          if(isl(1).eq.1) isl(2)=0
          if(isl(1).eq.0) islm=1
          if(isl(2).eq.0) islm=2
        endif
        if(isl(1).eq.1.or.isl(2).eq.1) goto 200
      endif
      if(igm.gt.0.and.mod(mste(12),2).eq.1.and.(p(n+1,5).ge.
     &pmth(2,ifld(1)).or.p(n+3,5).ge.pmth(2,ifld(2)))) then
        pmq1=p(n+2,5)/p(im+1,5)
        pmq2=p(n+4,5)/p(im+1,5)
        dz=sqrt((1.-p(im+1,5)/pem**2)*((1.-pmq1-pmq2)**2-4.*pmq1*pmq2))
        zh=1.+pmq1-pmq2
        zl=0.5*(zh-dz)
        zu=0.5*(zh+dz)
        if(p(im+1,1).lt.zl.or.p(im+1,1).gt.zu) goto 200
      endif

c...accepted branch, construct four momentum
  320 if(nep.eq.1) then
        p(n+1,1)=0.
        p(n+1,2)=0.
        p(n+1,3)=sqrt(max(0.,(p(ipa(1),4)+p(n+1,5))*(p(ipa(1),4)
     &  -p(n+1,5))))
        p(n+1,4)=p(ipa(1),4)
        p(n+2,2)=p(n+1,4)
      elseif(igm.eq.0.and.nep.eq.2) then
        ped1=0.5*(p(im+1,5)+p(n+2,5)-p(n+4,5))/p(im,5)
        p(n+1,1)=0.
        p(n+1,2)=0.
        p(n+1,3)=sqrt(max(0.,(ped1+p(n+1,5))*(ped1-p(n+1,5))))
        p(n+1,4)=ped1
        p(n+3,1)=0.
        p(n+3,2)=0.
        p(n+3,3)=-p(n+1,3)
        p(n+3,4)=p(im,5)-ped1
        p(n+2,2)=p(n+1,4)
        p(n+4,2)=p(n+3,4)
      elseif(nep.eq.3) then
        p(n+1,1)=0.
        p(n+1,2)=0.
        p(n+1,3)=sqrt(max(0.,pa1s))
        p(n+3,1)=sqrt(pts)
        p(n+3,2)=0.
        p(n+3,3)=0.5*(pa3s-pa2s-pa1s)/p(n+1,3)
        p(n+5,1)=-p(n+3,1)
        p(n+5,2)=0.
        p(n+5,3)=-(p(n+3,3)+p(n+1,3))
        p(n+2,2)=p(n+1,4)
        p(n+4,2)=p(n+3,4)
        p(n+6,2)=p(n+5,4)
      elseif(nep.eq.4) then
      else
        zm=p(im+1,1)
        pzm=sqrt(max(0.,(pem+p(im,5))*(pem-p(im,5))))
        pmls=(p(im+1,5)-p(n+2,5)-p(n+4,5))**2-4.*p(n+2,5)*p(n+4,5)
        if(mod(mste(12),2).eq.1) then
          if(pzm.gt.0.) then
            pts=(pem**2*(zm*(1.-zm)*p(im+1,5)-(1.-zm)*p(n+2,5)-
     &      zm*p(n+4,5))-0.25*pmls)/pzm**2
          else
            pts=0.
          endif
          p(n+1,4)=pem*p(im+1,1)
        else
          if(pzm.gt.0.) then
            pts=pmls*(zm*(1.-zm)*pem**2/p(im+1,5)-0.25)/pzm**2
          else
            pts=0.
          endif
          p(n+1,4)=pem*(0.5*(p(im+1,5)-sqrt(pmls)+p(n+2,5)-p(n+4,5))+
     &    sqrt(pmls)*zm)/p(im+1,5)
        endif
        pt=sqrt(max(0.,pts))
        phi=par(72)*rluC(0)
        p(n+1,1)=pt*cos(phi)
        p(n+1,2)=pt*sin(phi)
        if(pzm.gt.0.) then
          p(n+1,3)=0.5*(p(n+4,5)-p(n+2,5)-p(im+1,5)+2.*pem*p(n+1,4))/
     &    pzm
        else
          p(n+1,3)=0.
        endif
        p(n+3,1)=-p(n+1,1)
        p(n+3,2)=-p(n+1,2)
        p(n+3,3)=pzm-p(n+1,3)
        p(n+3,4)=pem-p(n+1,4)
        if(mste(12).le.2) then
          p(n+2,2)=(pem*p(n+1,4)-pzm*p(n+1,3))/p(im,5)
          p(n+4,2)=(pem*p(n+3,4)-pzm*p(n+3,3))/p(im,5)
        endif
      endif

c...rotate and boost line n+1 and n+3
      if(igm.gt.0) then
        if(mste(12).le.2) then
          bex=p(igm,1)/p(igm,4)
          bey=p(igm,2)/p(igm,4)
          bez=p(igm,3)/p(igm,4)
          ga=p(igm,4)/p(igm,5)
          gabep=ga*(ga*(bex*p(im,1)+bey*p(im,2)+bez*p(im,3))/(1.+ga)-
     &    p(im,4))
        else
          bex=0.
          bey=0.
          bez=0.
          gabep=0.
        endif
        the=ulanglC(p(im,3)+gabep*bez,sqrt((p(im,1)+gabep*bex)**2+
     &  (p(im,2)+gabep*bey)**2))
        phi=ulanglC(p(im,1)+gabep*bex,p(im,2)+gabep*bey)
        mst(1)=n+1
        mst(2)=mst(1)
        call luroboC(the,phi,bex,bey,bez)
        mst(1)=n+3
        mst(2)=mst(1)
        call luroboC(the,phi,bex,bey,bez)
        mst(1)=0
        mst(2)=0
      endif

c...continue loop over partons that may branch until none left
      if(igm.ge.0) k(im,1)=k(im,1)+20000
      n=n+2*nep
      nep=2
      goto 140
  330 continue

c...reconstruct string drawing information
      do 340 i=ns+1,n-1,2
      do 340 j=1,5
  340 p(i+1,j)=0.
      if(npa.ge.2) then
        k(ns+2,1)=70000+ns+1
        k(ns+2,2)=0
        iim=2
      else
        iim=0
      endif
      do 360 i=ns+1+iim,n-1,2
      if(k(i,1).lt.20000) goto 350
      id1=k(i+1,1)
      if(k(i,2).ge.501) id1=k(i+1,1)+2
      id2=2*k(i+1,1)+2-id1
      p(i+1,3)=id1
      p(i+1,4)=id2
      p(id1+1,1)=i
      p(id1+1,2)=id2
      p(id2+1,1)=id1
      p(id2+1,2)=i
  350 k(i+1,1)=70000+i
      k(i+1,2)=k(mod(k(i,1),20000)+1,2)
      if(i.eq.ns+1+iim) k(i+1,2)=k(ipa(1)+1,2)
      if(i.eq.ns+5.and.npa.ge.2) k(i+1,2)=k(ipa(2)+1,2)
      if(i.eq.ns+7.and.npa.ge.3) k(i+1,2)=k(ipa(3)+1,2)
  360 if(i.eq.ns+9.and.npa.eq.4) k(i+1,2)=k(ipa(4)+1,2)

c...transformation from cm frame
      if(npa.ge.2) then
        bex=ps(1)/ps(4)
        bey=ps(2)/ps(4)
        bez=ps(3)/ps(4)
        ga=ps(4)/ps(5)
        gabep=ga*(ga*(bex*p(ipa(1),1)+bey*p(ipa(1),2)+bez*p(ipa(1),3))
     &  /(1.+ga)-p(ipa(1),4))
      else
        bex=0.
        bey=0.
        bez=0.
        gabep=0.
      endif
      the=ulanglC(p(ipa(1),3)+gabep*bez,sqrt((p(ipa(1),1)
     &+gabep*bex)**2+(p(ipa(1),2)+gabep*bey)**2))
      phi=ulanglC(p(ipa(1),1)+gabep*bex,p(ipa(1),2)+gabep*bey)
      mst(1)=ns+1
      if(npa.eq.3) then
        chi=ulanglC(cos(the)*cos(phi)*(p(ipa(2),1)+gabep*bex)+cos(the)*
     &  sin(phi)*(p(ipa(2),2)+gabep*bey)-sin(the)*(p(ipa(2),3)+gabep*
     &  bez),-sin(phi)*(p(ipa(2),1)+gabep*bex)+cos(phi)*(p(ipa(2),2)+
     &  gabep*bey))
        call luroboC(0.,chi,0.,0.,0.)
      endif
      call luroboC(the,phi,bex,bey,bez)
      mst(1)=0

c...delete trivial shower, else connect initiators
      if(n.eq.ns+2+2*iim) then
        n=ns
      else
        do 370 ip=1,npa
        k(ipa(ip),1)=k(ipa(ip),1)+20000
        p(ipa(ip)+1,3)=ns+iim+2*ip-1
        p(ipa(ip)+1,4)=ns+iim+2*ip-1
        p(ns+2*ip+iim,1)=ipa(ip)
  370   p(ns+2*ip+iim,2)=ipa(ip)
      endif

      return
      end

c*********************************************************************

      subroutine luspheC(sph,apl)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)
      dimension sm(3,3),sv(3,3)

      np=0
c...calculate matrix to be diagonalized
      do 100 l1=1,3
      do 100 l2=l1,3
  100 sm(l1,l2)=0.
      ps=0.
      do 120 i=1,n
      if(k(i,1).ge.20000) goto 120
      np=np+1
      pa=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
      pwt=1.
      if(abs(pare(30)-2.).gt.0.001) pwt=pa**(pare(30)-2.)
      do 110 l1=1,3
      do 110 l2=l1,3
  110 sm(l1,l2)=sm(l1,l2)+pwt*p(i,l1)*p(i,l2)
      ps=ps+pwt*pa**2
  120 continue

      if(np.le.1) then
c...very low multiplicities (0 or 1) not considered
        sph=-1.
        apl=-1.
        return
      endif
      do 130 l1=1,3
      do 130 l2=l1,3
  130 sm(l1,l2)=sm(l1,l2)/ps

c...find eigenvalues to matrix (third degree equation)
      sq=(sm(1,1)*sm(2,2)+sm(1,1)*sm(3,3)+sm(2,2)*sm(3,3)-sm(1,2)**2-
     &sm(1,3)**2-sm(2,3)**2)/3.-1./9.
      sr=-0.5*(sq+1./9.+sm(1,1)*sm(2,3)**2+sm(2,2)*sm(1,3)**2+sm(3,3)*
     &sm(1,2)**2-sm(1,1)*sm(2,2)*sm(3,3))+sm(1,2)*sm(1,3)*sm(2,3)+1./27.
      sp=cos(acos(max(min(sr/sqrt(-sq**3),1.),-1.))/3.)
      p(n+1,4)=1./3.+sqrt(-sq)*max(2.*sp,sqrt(3.*(1.-sp**2))-sp)
      p(n+3,4)=1./3.+sqrt(-sq)*min(2.*sp,-sqrt(3.*(1.-sp**2))-sp)
      p(n+2,4)=1.-p(n+1,4)-p(n+3,4)
      if(p(n+2,4).lt.1e-7) then
        sph=-1.
        apl=-1.
        return
      endif

c...find first and last eigenvector by solving equation system
      do 170 ld=1,3,2
      do 140 l1=1,3
      sv(l1,l1)=sm(l1,l1)-p(n+ld,4)
      do 140 l2=l1+1,3
      sv(l1,l2)=sm(l1,l2)
  140 sv(l2,l1)=sm(l1,l2)
      smax=0.
      do 150 l1=1,3
      do 150 l2=1,3
      if(abs(sv(l1,l2)).le.smax) goto 150
      li=l1
      lj=l2
      smax=abs(sv(l1,l2))
  150 continue
      smax=0.
      do 160 l3=li+1,li+2
      l1=l3-3*((l3-1)/3)
      rl=sv(l1,lj)/sv(li,lj)
      do 160 l2=1,3
      sv(l1,l2)=sv(l1,l2)-rl*sv(li,l2)
      if(abs(sv(l1,l2)).le.smax) goto 160
      lk=l1
      smax=abs(sv(l1,l2))
  160 continue
      lj1=lj+1-3*(lj/3)
      lj2=lj+2-3*((lj+1)/3)
      p(n+ld,lj1)=-sv(lk,lj2)
      p(n+ld,lj2)=sv(lk,lj1)
      p(n+ld,lj)=-(sv(li,lj1)*p(n+ld,lj1)+sv(li,lj2)*p(n+ld,lj2))/
     &sv(li,lj)
      pa=sqrt(p(n+ld,1)**2+p(n+ld,2)**2+p(n+ld,3)**2)
      sgn=(-1.)**int(rluC(0)+0.5)
      do 170 j=1,3
  170 p(n+ld,j)=sgn*p(n+ld,j)/pa

c...middle eigenvector orthogonal to other two, reset unused components
      sgn=(-1.)**int(rluC(0)+0.5)
      p(n+2,1)=sgn*(p(n+1,2)*p(n+3,3)-p(n+1,3)*p(n+3,2))
      p(n+2,2)=sgn*(p(n+1,3)*p(n+3,1)-p(n+1,1)*p(n+3,3))
      p(n+2,3)=sgn*(p(n+1,1)*p(n+3,2)-p(n+1,2)*p(n+3,1))
      do 180 ld=1,3
      k(n+ld,1)=ld
      k(n+ld,2)=0
  180 p(n+ld,5)=0.

      sph=1.5*(p(n+2,4)+p(n+3,4))
      apl=1.5*p(n+3,4)
      mst(3)=3

      return
      end

c*********************************************************************

      subroutine luthruC(thr,obl)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)
      dimension tdi(3),tpr(3)

      np=0
      ps=0.
      do 280 ld=1,2
      if(ld.eq.2) then
c...thrust axis along z direction for major axis search
        mst(2)=n+1
        phi=ulanglC(p(n+1,1),p(n+1,2))
        call luroboC(0.,-phi,0.,0.,0.)
        the=ulanglC(p(n+1,3),p(n+1,1))
        call luroboC(-the,0.,0.,0.,0.)
      endif

c...find and order particles with highest p (pt for major)
c...(p(i,5) is temporarily used for extra particle weight, 1 for thrust)
      if(mst(23).ge.1.and.n+mste(21)/10+15.ge.mst(30)-5-mst(31)) then
        thr=-2.
        obl=-2.
        return
      endif
      do 100 lf=n+4,n+mste(21)/10+4
  100 p(lf,4)=0.
      do 140 i=1,n
      if(k(i,1).ge.20000) goto 140
      if(ld.eq.1) then
        np=np+1
        pa=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
        p(i,5)=1.
        if(abs(pare(31)-1.).gt.0.001) p(i,5)=pa**(pare(31)-1.)
        ps=ps+p(i,5)*pa
      else
        pa=sqrt(p(i,1)**2+p(i,2)**2)
      endif
      do 110 lf=n+mste(21)/10+3,n+4,-1
      if(pa.le.p(lf,4)) goto 120
      do 110 j=1,5
  110 p(lf+1,j)=p(lf,j)
      lf=n+3
  120 do 130 j=1,3
  130 p(lf+1,j)=p(i,j)
      p(lf+1,4)=pa
      p(lf+1,5)=p(i,5)
  140 continue

      if(np.le.1) then
c...very low multiplicities (0 or 1) not considered
        thr=-1.
        obl=-1.
        return
      endif

c...find and order initial axes with highest thrust
      do 150 lg=n+mste(21)/10+5,n+mste(21)/10+15
  150 p(lg,4)=0.
      nc=2**(min(mste(21)/10,np)-1)
      do 210 lc=1,nc
      do 160 j=1,3
  160 tdi(j)=0.
      do 170 lf=1,min(mste(21)/10,np)
      sgn=p(n+lf+3,5)
      if(2**lf*((lc+2**(lf-1)-1)/2**lf).ge.lc) sgn=-sgn
      do 170 j=1,4-ld
  170 tdi(j)=tdi(j)+sgn*p(n+lf+3,j)
      tds=tdi(1)**2+tdi(2)**2+tdi(3)**2
      do 180 lg=n+mste(21)/10+min(lc,10)+4,n+mste(21)/10+5,-1
      if(tds.le.p(lg,4)) goto 190
      do 180 j=1,4
  180 p(lg+1,j)=p(lg,j)
      lg=n+mste(21)/10+4
  190 do 200 j=1,3
  200 p(lg+1,j)=tdi(j)
      p(lg+1,4)=tds
  210 continue

c...iterate direction of axis until stable maximum
      p(n+ld,4)=0.
      lg=0
  220 lg=lg+1
      thp=0.
  230 thps=thp
      do 240 j=1,3
      if(thp.le.1e-10) tdi(j)=p(n+mste(21)/10+4+lg,j)
      if(thp.gt.1e-10) tdi(j)=tpr(j)
  240 tpr(j)=0.
      do 260 i=1,n
      if(k(i,1).ge.20000) goto 260
      sgn=sign(p(i,5),tdi(1)*p(i,1)+tdi(2)*p(i,2)+tdi(3)*p(i,3))
      do 250 j=1,4-ld
  250 tpr(j)=tpr(j)+sgn*p(i,j)
  260 continue
      thp=sqrt(tpr(1)**2+tpr(2)**2+tpr(3)**2)/ps
      if(thp.ge.thps+pare(34)) goto 230

c...save good axis, try new initial axis until a number of tries agree
      if(thp.lt.p(n+ld,4)-pare(34).and.lg.lt.min(10,nc)) goto 220
      if(thp.gt.p(n+ld,4)+pare(34)) then
        lagr=0
        sgn=(-1.)**int(rluC(0)+0.5)
        do 270 j=1,3
  270   p(n+ld,j)=sgn*tpr(j)/(ps*thp)
        p(n+ld,4)=thp
      endif
      lagr=lagr+1
  280 if(lagr.lt.mod(mste(21),10).and.lg.lt.min(10,nc)) goto 220

c...find minor axis and value by orthogonality
      sgn=(-1.)**int(rluC(0)+0.5)
      p(n+3,1)=-sgn*p(n+2,2)
      p(n+3,2)=sgn*p(n+2,1)
      p(n+3,3)=0.
      thp=0.
      do 290 i=1,n
      if(k(i,1).ge.20000) goto 290
      thp=thp+p(i,5)*abs(p(n+3,1)*p(i,1)+p(n+3,2)*p(i,2))
      p(i,5)=sqrt(max(p(i,4)**2-p(i,1)**2-p(i,2)**2-p(i,3)**2,0.))
  290 continue
      p(n+3,4)=thp/ps

c...reset unused components, rotate back to original coordinate system
      do 300 ld=1,3
      k(n+ld,1)=ld
      k(n+ld,2)=0
  300 p(n+ld,5)=0.
      mst(2)=n+3
      call luroboC(the,phi,0.,0.,0.)
      mst(2)=0

      thr=p(n+1,4)
      obl=p(n+2,4)-p(n+3,4)
      mst(3)=3

      return
      end

c*********************************************************************

      subroutine luclusC(njet,tgen,dmin)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)

c...momenta and sum of momenta for particles
c...(p(i,4) is temporarily used to represent absolute momenta)
      np=0
      ps=0.
      do 100 i=1,n
      if(k(i,1).ge.20000) goto 100
      np=np+1
      p(i,4)=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
      ps=ps+p(i,4)
  100 continue

      if(np.le.2*iabs(mste(22))) then
c...very low multiplicities not considered
        njet=-1
        tgen=-1.
        dmin=-1.
        return
      endif
      nl=0

      if(mste(22).ge.0) then
c...find initial jet configuration. if too few jets, make harder cuts
        dinit=1.25*pare(32)
  110   dinit=0.8*dinit

c...sum up small momentum region, jet if enough absolute momentum
        njet=0
        na=0
        do 120 j=1,3
  120   p(n+1,j)=0.
        do 140 i=1,n
        if(k(i,1).ge.20000) goto 140
        k(i,1)=0
        if(p(i,4).gt.2.*dinit) goto 140
        na=na+1
        k(i,1)=1
        do 130 j=1,3
  130   p(n+1,j)=p(n+1,j)+p(i,j)
  140   continue
        p(n+1,4)=sqrt(p(n+1,1)**2+p(n+1,2)**2+p(n+1,3)**2)
        if(p(n+1,4).gt.2.*dinit) njet=1
        if(dinit.ge.0.2*pare(32).and.njet+np-na.lt.2*iabs(mste(22)))
     &  goto 110

c...find fastest particle, sum up jet around it. iterate until all
c...particles used up
  150   njet=njet+1
        if(mst(23).ge.1.and.n+2*njet.ge.mst(30)-5-mst(31)) then
          njet=-2
          tgen=-2.
          dmin=-2.
          return
        endif
        pmax=0.
        do 160 i=1,n
        if(k(i,1).ne.0.or.p(i,4).le.pmax) goto 160
        im=i
        pmax=p(i,4)
  160   continue
        do 170 j=1,3
  170   p(n+njet,j)=0.
        do 190 i=1,n
        if(k(i,1).ne.0) goto 190
        d2=(p(i,4)*p(im,4)-p(i,1)*p(im,1)-p(i,2)*p(im,2)-
     &  p(i,3)*p(im,3))*2.*p(i,4)*p(im,4)/(p(i,4)+p(im,4))**2
        if(d2.gt.dinit**2) goto 190
        na=na+1
        k(i,1)=njet
        do 180 j=1,3
  180   p(n+njet,j)=p(n+njet,j)+p(i,j)
  190   continue
        p(n+njet,4)=sqrt(p(n+njet,1)**2+p(n+njet,2)**2+p(n+njet,3)**2)
        if(dinit.ge.0.2*pare(32).and.njet+np-na.lt.2*iabs(mste(22)))
     &  goto 110
        if(na.lt.np) goto 150

      else
c...use given initial jet configuration
        do 200 it=n+1,n+njet
  200   p(it,4)=sqrt(p(it,1)**2+p(it,2)**2+p(it,3)**2)
      endif

c...assign all particles to nearest jet, sum up new jet momenta
  210 tsav=0.
  220 do 230 it=n+njet+1,n+2*njet
      do 230 j=1,3
  230 p(it,j)=0.
      do 270 i=1,n
      if(k(i,1).ge.20000) goto 270
      if(mste(23).eq.1) then
c...symmetric distance measure between particle and jet
        d2min=1e10
        do 240 it=n+1,n+njet
        if(p(it,4).lt.dinit) goto 240
        d2=(p(i,4)*p(it,4)-p(i,1)*p(it,1)-p(i,2)*p(it,2)-
     &  p(i,3)*p(it,3))*2.*p(i,4)*p(it,4)/(p(i,4)+p(it,4))**2
        if(d2.ge.d2min) goto 240
        im=it
        d2min=d2
  240   continue
      else
c..."multicity" distance measure between particle and jet
        pmax=-1e10
        do 250 it=n+1,n+njet
        if(p(it,4).lt.dinit) goto 250
        prod=(p(i,1)*p(it,1)+p(i,2)*p(it,2)+p(i,3)*p(it,3))/p(it,4)
        if(prod.le.pmax) goto 250
        im=it
        pmax=prod
  250   continue
      endif
      k(i,1)=im-n
      do 260 j=1,3
  260 p(im+njet,j)=p(im+njet,j)+p(i,j)
  270 continue

c...absolute value and sum of jet momenta, find two closest jets
      psjt=0.
      do 280 it=n+njet+1,n+2*njet
      p(it,4)=sqrt(p(it,1)**2+p(it,2)**2+p(it,3)**2)
  280 psjt=psjt+p(it,4)
      d2min=1e10
      do 290 it1=n+njet+1,n+2*njet-1
      do 290 it2=it1+1,n+2*njet
      d2=(p(it1,4)*p(it2,4)-p(it1,1)*p(it2,1)-p(it1,2)*p(it2,2)-
     &p(it1,3)*p(it2,3))*2.*p(it1,4)*p(it2,4)/
     &max(0.01,p(it1,4)+p(it2,4))**2
      if(d2.ge.d2min) goto 290
      im1=it1
      im2=it2
      d2min=d2
  290 continue

c...if allowed, join two closest jets and start over
      if(njet.gt.iabs(mste(22)).and.d2min.lt.pare(33)**2) then
        nr=1
        do 300 j=1,3
  300   p(n+nr,j)=p(im1,j)+p(im2,j)
        p(n+nr,4)=sqrt(p(n+nr,1)**2+p(n+nr,2)**2+p(n+nr,3)**2)
        do 320 it=n+njet+1,n+2*njet
        if(it.eq.im1.or.it.eq.im2) goto 320
        nr=nr+1
        do 310 j=1,5
  310   p(n+nr,j)=p(it,j)
  320   continue
        njet=njet-1
        goto 210

c...divide up broad jet if empty cluster in list of final ones
      elseif(njet.eq.iabs(mste(22)).and.nl.le.2) then
        do 330 it=n+1,n+njet
  330   k(it,2)=0
        do 340 i=1,n
  340   if(k(i,1).lt.20000) k(n+k(i,1),2)=k(n+k(i,1),2)+1
        im=0
        do 350 it=n+1,n+njet
  350   if(k(it,2).eq.0) im=it
        if(im.ne.0) then
          nl=nl+1
          ir=0
          d2max=0.
          do 360 i=1,n
          if(k(i,1).ge.20000) goto 360
          if(k(n+k(i,1),2).le.1.or.p(i,4).lt.dinit) goto 360
          it=n+njet+k(i,1)
          d2=(p(i,4)*p(it,4)-p(i,1)*p(it,1)-p(i,2)*p(it,2)-
     &    p(i,3)*p(it,3))*2.*p(i,4)*p(it,4)/(p(i,4)+p(it,4))**2
          if(d2.le.d2max) goto 360
          ir=i
          d2max=d2
  360     continue
          if(ir.eq.0) goto 390
          it=n+njet+k(ir,1)
          do 370 j=1,3
          p(im+njet,j)=p(ir,j)
  370     p(it,j)=p(it,j)-p(ir,j)
          p(im+njet,4)=p(ir,4)
          p(it,4)=sqrt(p(it,1)**2+p(it,2)**2+p(it,3)**2)
          do 380 it=n+1,n+njet
          do 380 j=1,5
  380     p(it,j)=p(it+njet,j)
          if(nl.le.2) goto 210
        endif
      endif

c...if generalized thrust has not yet converged, continue iteration
  390 tgen=psjt/ps
      if(tgen.gt.tsav+pare(34).and.nl.le.2) then
        tsav=tgen
        do 400 it=n+1,n+njet
        do 400 j=1,5
  400   p(it,j)=p(it+njet,j)
        goto 220
      endif

c...reorder jets after momentum, sum up jet energies and multiplicities
      do 420 it=n+1,n+njet
      pmax=0.
      do 410 ir=n+njet+1,n+2*njet
      if(p(ir,4).le.pmax) goto 410
      im=ir
      pmax=p(ir,4)
  410 continue
      k(im,1)=it-n
      p(im,4)=-1.
      k(it,1)=it-n
      k(it,2)=0
      p(it,4)=0.
      do 420 j=1,3
  420 p(it,j)=p(im,j)
      do 430 i=1,n
      if(k(i,1).ge.20000) goto 430
      k(i,1)=k(n+njet+k(i,1),1)
      p(i,4)=sqrt(p(i,5)**2+p(i,1)**2+p(i,2)**2+p(i,3)**2)
      k(n+k(i,1),2)=k(n+k(i,1),2)+1
      p(n+k(i,1),4)=p(n+k(i,1),4)+p(i,4)
  430 continue
      im=0
      do 440 it=n+1,n+njet
      if(k(it,2).eq.0) im=it
  440 p(it,5)=sqrt(max(p(it,4)**2-p(it,1)**2-p(it,2)**2-p(it,3)**2,0.))

c...values at return (negative for failure fixed number of clusters)
      dmin=sqrt(d2min)
      if(njet.eq.1) dmin=0.
      mst(3)=njet
      if(im.ne.0) then
        njet=-1
        tgen=-1.
        dmin=-1.
      endif

      return
      end

c*********************************************************************

      subroutine luorieC(mori)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      dimension ns(2),pts(2),pls(2)

c...place largest axis along z axis and second largest in xy plane
      mst(2)=n+mst(3)
      call luroboC(0.,-ulanglC(p(n+1,1),p(n+1,2)),0.,0.,0.)
      call luroboC(-ulanglC(p(n+1,3),p(n+1,1)),0.,0.,0.,0.)
      call luroboC(0.,-ulanglC(p(n+2,1),p(n+2,2)),0.,0.,0.)
      if(mori.eq.1) mst(2)=0
      if(mori.eq.1) return

      do 100 is=1,2
      ns(is)=0
      pts(is)=0.
  100 pls(is)=0.

c...rotate slim jet along +z direction
      do 110 i=1,n
      if(k(i,1).ge.20000) goto 110
      is=2.-sign(0.5,p(i,3))
      ns(is)=ns(is)+1
      pts(is)=pts(is)+sqrt(p(i,1)**2+p(i,2)**2)
  110 continue
      if(ns(1)*pts(2)**2.lt.ns(2)*pts(1)**2)
     &call luroboC(par(71),0.,0.,0.,0.)

c...rotate second largest jet into -z,+x quadrant
      do 120 i=1,n
      if(k(i,1).ge.20000.or.p(i,3).ge.0.) goto 120
      is=2.-sign(0.5,p(i,1))
      pls(is)=pls(is)-p(i,3)
  120 continue
      if(pls(2).gt.pls(1)) call luroboC(0.,par(71),0.,0.,0.)
      mst(2)=0

      return
      end

c*********************************************************************

      subroutine lucellC(njet)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)

c...loop over all particles: find cell that was hit
      nce2=2*mste(24)*mste(25)
      ptlrat=1./sinh(pare(35))**2
      nc=n
      do 110 i=1,n
      if(k(i,1).ge.20000) goto 110
      if(p(i,1)**2+p(i,2)**2.le.ptlrat*p(i,3)**2) goto 110
      pt=sqrt(p(i,1)**2+p(i,2)**2)
      eta=sign(alog((sqrt(pt**2+p(i,3)**2)+abs(p(i,3)))/pt),p(i,3))
      ieta=max(1,min(mste(24),1+int(mste(24)*0.5*(eta/pare(35)+1.))))
      phi=ulanglC(p(i,1),p(i,2))
      iphi=max(1,min(mste(25),1+int(mste(25)*0.5*(phi/par(71)+1.))))
      ietph=mste(25)*ieta+iphi

c...add to cell already hit, or book new cell
      do 100 ic=n+1,nc
      if(ietph.eq.k(ic,1)) then
        k(ic,2)=k(ic,2)+1
        p(ic,5)=p(ic,5)+pt
        goto 110
      endif
  100 continue
      if(mst(23).ge.1.and.nc.ge.mst(30)-5-mst(31)) then
        njet=-2
        return
      endif
      nc=nc+1
      k(nc,1)=ietph
      k(nc,2)=1
      p(nc,1)=(pare(35)/mste(24))*(2*ieta-1-mste(24))
      p(nc,2)=(par(71)/mste(25))*(2*iphi-1-mste(25))
      p(nc,5)=pt
  110 continue

c...smear true bin content by calorimeter resolution
      if(mste(27).ge.1) then
        do 130 ic=n+1,nc
        pei=p(ic,5)
        if(mste(27).eq.2) pei=p(ic,5)/cosh(p(ic,1))
  120   pef=pei+pare(39)*sqrt(-2.*alog(max(1e-10,rluC(0)))*pei)*
     &  cos(par(72)*rluC(0))
        if(pef.lt.0..or.pef.gt.pare(40)*pei) goto 120
        p(ic,5)=pef
  130   if(mste(27).eq.2) p(ic,5)=pef*cosh(p(ic,1))
      endif

c...find initiator cell, the one with highest pt of not yet used ones
      nj=nc
  140 etmax=0.
      do 150 ic=n+1,nc
      if(k(ic,1).eq.0.or.k(ic,1).gt.nce2) goto 150
      if(p(ic,5).le.etmax) goto 150
      icmax=ic
      eta=p(ic,1)
      phi=p(ic,2)
      etmax=p(ic,5)
  150 continue
      if(etmax.lt.pare(36)) goto 210
      if(mst(23).ge.1.and.nj.ge.mst(30)-5-mst(31)) then
        njet=-2
        return
      endif
      k(icmax,1)=k(icmax,1)+nce2
      nj=nj+1
      k(nj,1)=1
      k(nj,2)=0
      p(nj,1)=eta
      p(nj,2)=phi
      p(nj,3)=0.
      p(nj,4)=0.
      p(nj,5)=0.

c...sum up unused cells within required distance of initiator
      do 160 ic=n+1,nc
      if(k(ic,1).eq.0) goto 160
      if(abs(p(ic,1)-eta).gt.pare(38)) goto 160
      dphia=abs(p(ic,2)-phi)
      if(dphia.gt.pare(38).and.dphia.lt.par(72)-pare(38)) goto 160
      phic=p(ic,2)
      if(dphia.gt.par(71)) phic=phic+sign(par(72),phi)
      if((p(ic,1)-eta)**2+(phic-phi)**2.gt.pare(38)**2) goto 160
      k(ic,1)=-k(ic,1)
      k(nj,2)=k(nj,2)+k(ic,2)
      p(nj,3)=p(nj,3)+p(ic,5)*p(ic,1)
      p(nj,4)=p(nj,4)+p(ic,5)*phic
      p(nj,5)=p(nj,5)+p(ic,5)
  160 continue

c...reject cluster below minimum et, else accept
      if(p(nj,5).lt.pare(37)) then
        nj=nj-1
        do 170 ic=n+1,nc
  170   if(k(ic,1).lt.0) k(ic,1)=-k(ic,1)
      elseif(mste(26).le.2) then
        p(nj,3)=p(nj,3)/p(nj,5)
        p(nj,4)=p(nj,4)/p(nj,5)
        if(abs(p(nj,4)).gt.par(71)) p(nj,4)=p(nj,4)-sign(par(72),
     &  p(nj,4))
        do 180 ic=n+1,nc
  180   if(k(ic,1).lt.0) k(ic,1)=0
      else
        do 190 j=1,4
  190   p(nj,j)=0.
        do 200 ic=n+1,nc
        if(k(ic,1).ge.0) goto 200
        p(nj,1)=p(nj,1)+p(ic,5)*cos(p(ic,2))
        p(nj,2)=p(nj,2)+p(ic,5)*sin(p(ic,2))
        p(nj,3)=p(nj,3)+p(ic,5)*sinh(p(ic,1))
        p(nj,4)=p(nj,4)+p(ic,5)*cosh(p(ic,1))
        k(ic,1)=0
  200   continue
      endif
      goto 140

c...arrange clusters in falling et sequence
  210 do 230 i=1,nj-nc
      etmax=0.
      do 220 ij=nc+1,nj
      if(k(ij,1).eq.0) goto 220
      if(p(ij,5).lt.etmax) goto 220
      ijmax=ij
      etmax=p(ij,5)
  220 continue
      k(ijmax,1)=0
      k(n+i,1)=i
      k(n+i,2)=k(ijmax,2)
      do 230 j=1,5
  230 p(n+i,j)=p(ijmax,j)
      njet=nj-nc
      mst(3)=njet

c...convert to massless or massive four-vectors
      if(mste(26).eq.2) then
        do 240 i=n+1,n+njet
        eta=p(i,3)
        p(i,1)=p(i,5)*cos(p(i,4))
        p(i,2)=p(i,5)*sin(p(i,4))
        p(i,3)=p(i,5)*sinh(eta)
        p(i,4)=p(i,5)*cosh(eta)
  240   p(i,5)=0.
      elseif(mste(26).ge.3) then
        do 250 i=n+1,n+njet
  250   p(i,5)=sqrt(max(0.,p(i,4)**2-p(i,1)**2-p(i,2)**2-p(i,3)**2))
      endif

      return
      end

c*********************************************************************

      subroutine lufowoC(h10,h20,h30,h40)
      common/lujetsC/n,k(2000,2),p(2000,5)

c...momenta for particles (p(i,5) temporarily used) and h0
      np=0
      h0=0
      hd=0.
      do 100 i=1,n
      if(k(i,1).ge.20000) goto 100
      np=np+1
      p(i,5)=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
      h0=h0+p(i,5)
      hd=hd+p(i,5)**2
  100 continue
      h0=h0**2

      if(np.le.1) then
c...very low multiplicities (0 or 1) not considered
        h10=-1.
        h20=-1.
        h30=-1.
        h40=-1.
        return
      endif

c...calculate h1 - h4
      h10=0.
      h20=0.
      h30=0.
      h40=0.
      do 120 i1=1,n-1
      if(k(i1,1).ge.20000) goto 120
      do 110 i2=i1+1,n
      if(k(i2,1).ge.20000) goto 110
      cthe=(p(i1,1)*p(i2,1)+p(i1,2)*p(i2,2)+p(i1,3)*p(i2,3))/
     &(p(i1,5)*p(i2,5))
      h10=h10+p(i1,5)*p(i2,5)*cthe
      h20=h20+p(i1,5)*p(i2,5)*(1.5*cthe**2-0.5)
      h30=h30+p(i1,5)*p(i2,5)*(2.5*cthe**3-1.5*cthe)
      h40=h40+p(i1,5)*p(i2,5)*(4.375*cthe**4-3.75*cthe**2+0.375)
  110 continue
  120 continue

c...calculate h10 - h40, reset p(i,5) to mass
      h10=(hd+2.*h10)/h0
      h20=(hd+2.*h20)/h0
      h30=(hd+2.*h30)/h0
      h40=(hd+2.*h40)/h0
      do 130 i=1,n
  130 if(k(i,1).lt.20000) p(i,5)=sqrt(max(p(i,4)**2-p(i,1)**2-
     &p(i,2)**2-p(i,3)**2,0.))

      return
      end

c*********************************************************************

      function ulalpsC(q2)
      common/ludat1C/mst(40),par(80)
      common/ludateC/mste(40),pare(80)

c...the number of active flavours
      nf=3
      do 100 ifl=4,mste(4)
  100 if(q2.gt.4.*ulmassC(2,ifl)**2) nf=nf+1

c...the strong coupling constant in first and second order
      if(mste(8).eq.0) then
        ulalpsC=pare(3)
      elseif(iabs(mste(1)).le.1.or.mste(8).eq.1) then
        ulalpsC=12.*par(71)/((33.-2.*nf)*alog(q2/pare(1)**2))
      else
        algq=alog(q2/pare(2)**2)
        ulalpsC=12.*par(71)/((33-2.*nf)*algq+(6.*(153.-19.*nf)/
     &  (33.-2.*nf))*alog(algq))
      endif

      pare(49)=ulalpsC
      pare(65)=1.986-0.115*nf

      return
      end

c*********************************************************************

      block data luedatC
      common/ludateC/mste(40),pare(80)
      data mste/
     1    3,    2,    7,    5,    1,    1,    0,    2,    1,    0,
     2    2,    4,    2,    2,    5,    0,    0,    0,    0,    0,
     3   42,    1,    1,   25,   24,    1,    0,    0,    0,    1,
     4    1,    1,    0,    0,    0,    0,    0,    0,    0,    0/
      data pare/
     1  1.5,  0.5, 0.20,0.0072974,0.229,94.,2.8, 0.02,  2.0,  1.0,
     2   0.,   0.,   0.,   0., 0.01, 0.99,  0.2,   0.,   0.,   0.,
     3 0.40,  1.0,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  2.0,
     4  1.0, 0.25,  2.5,0.0001, 2.5,  1.5,  7.0,  1.0,  0.5,  2.0,
     5  40*0./
      end

c*********************************************************************
c***  jetset version 6.3, low-pt physics part  ***********************

      subroutine luloptC(kf1,kf2,pe1,pe2)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/ludathC/chr(20),khr(60)
      dimension kre(2,3)

c...fill initial hadrons, find flavour configuration of two hadron jets
      do 110 ip=1,2
      kf=(2-ip)*kf1+(ip-1)*kf2
      pe=(2-ip)*pe1+(ip-1)*pe2
      call lupartC(ip,kf,pe,(ip-1)*par(71),0.)
      k(ip,1)=40000
      kfa=iabs(kf)
      if(kfa.ge.17.and.kfa.le.19) ihr=2*kfa-34
      if(kfa.eq.41.or.kfa.eq.42) ihr=5*kfa-199
      rhr=rluC(0)
  100 ihr=ihr+1
      if(rhr.gt.chr(ihr)) goto 100
      do 110 j=1,3
  110 kre(ip,j)=khr(3*(ihr-1)+j)*isign(1,kf)

c...calculate cm energy, fill jets in cm frame
      pm1=ulmassC(1,kf1)
      per1=max(pm1,pe1)
      pz1=sqrt(per1**2-pm1**2)
      pm2=ulmassC(1,kf2)
      per2=max(pm2,pe2)
      pz2=-sqrt(per2**2-pm2**2)
      ecm=sqrt((per1+per2)**2-(pz1+pz2)**2)
      pec1=0.5*(ecm+(ulmassC(2,kre(1,1))**2-ulmassC(2,kre(2,1))**2)/ecm)
      call lu1jetC(3,kre(1,1),kre(1,2),kre(1,3),pec1,0.,0.)
      k(3,1)=10001
      call lu1jetC(5,kre(2,1),kre(2,2),kre(2,3),ecm-pec1,par(71),0.)
      k(5,1)=2

c...jet fragmentation, boost to lab frame
      call luexecC
      mst(1)=3
      call luroboC(0.,0.,0.,0.,(pz1+pz2)/(per1+per2))
      mst(1)=0

      return
      end

c*********************************************************************

      block data luhdatC
      common/ludathC/chr(20),khr(60)
c...flavour arrangement data for hadrons
      data chr/0.5,1.,0.5,1.,0.5,1.,0.3333,0.7708,0.8333,0.9792,1.,
     &0.3333,0.7708,0.8333,0.9792,1.,4*0./
      data khr/1,0,-2,-2,0,1,1,0,-3,-3,0,1,2,0,-3,-3,0,2,11,0,2,12,2,
     &1,12,1,1,21,2,1,21,1,1,22,0,1,12,1,2,12,2,2,21,1,2,21,2,2,12*0/
      end

c***  end of jetset version 6.3  *************************************
c*********************************************************************

      function rluC(idum)
c...this function is an interface to a suitable random number
c...generator. some possibilities are suggested below, but in
c...the end it is up to the user to find one that works.
c################
      real*8  u
c################

       call rndc(u)
       rluC=u
c            avoid ~1 and 0; if( rluC .eq. 1.0 ) then
c                                   is not enough in some log
c                                   which result in 0. for 0.9999..
c
       if(rluC .ge. 0.999999)  then
          rluC = 0.999999
       elseif(rluC .eq. 0.) then
          rluC = 1.e-18
       endif


c...for nd (lund) and univac
c     rluC=ranf(0)

c...for vax (other seed may be chosen if desired)
c     data iseed/65539/
c     rluC=ran(iseed)

c...for ibm (desy)
c     rluC=rn(float(idum))

c...for ibm (slac)
c     rluC=ran7(float(idum))

c...for cdc
c     rluC=ranf()

c...using cern library
c     rluC=rndmC(idum)
      return
      end
