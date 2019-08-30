*cmz :  3.14/16 13/03/89  14.48.41  by  nick van eijndhoven (cern)
*-- author :
      subroutine poisso(xav,iran)
c
c *** generation of poisson distribution ***
c *** nve 16-mar-1988 cern geneva ***
c
      dimension rndm(1)
c origin : h.fesefeldt (27-oct-1983)
c
c --- use normal distribution for <x> > 9.9 ---
      if(xav.gt.9.9) goto 2
c
      mm=ifix(5.*xav)
      iran=0
      if(mm.le.0) goto 3
      r=exp(-xav)
      call grndm(rndm,1)
      ran1=rndm(1)
      if(ran1.le.r) return
      rr=r
      do 1 i=1,mm
      iran=iran+1
      if(i.le.5) rrr=xav**i/nfac(i)
c** stirling' s formula for large numbers
      if(i.gt.5) rrr=exp(i*log(xav)-(i+0.5)*log(i*1.)+i-0.9189385)
      rr=rr+r*rrr
      if(ran1.le.rr) return
    1 continue
      return
c** normal distribution with sigma**2 = <x>
    2 call normal(ran1)
      ran1=xav+ran1*sqrt(xav)
      iran=ifix(ran1)
      if(iran.lt.0) iran=0
      return
c** for very small xav try iran=1,2,3
    3 p1=xav*exp(-xav)
      p2=xav*p1/2.
      p3=xav*p2/3.
      call grndm(rndm,1)
      ran=rndm(1)
      iran=3
      if(ran.lt.p3) return
      iran=2
      if(ran.lt.p2) return
      iran=1
      if(ran.lt.p1) return
      iran=0
      return
      end
