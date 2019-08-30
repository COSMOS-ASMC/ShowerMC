*cmz :  3.14/16 13/03/89  14.48.44  by  nick van eijndhoven (cern)
*-- author :
      subroutine dlpsor(a,n)
c
c *** double precision version of cernlib routine "flpsor" ***
c *** nve 29-mar-1988 cern geneva ***
c
c called by : phasp
c origin : h.fesefeldt (02-dec-1986)
c
      implicit double precision (a-h,o-z)
c
      dimension a(n)
      dimension lt(20),rt(20)
      integer r,rt
c
      level=1
      lt(1)=1
      rt(1)=n
   10 l=lt(level)
      r=rt(level)
      level=level-1
   20 if(r.le.l) if(level) 50,50,10
c
c   subdivide the interval l,r
c     l : lower limit of the interval (input)
c     r : upper limit of the interval (input)
c     j : upper limit of lower sub-interval (output)
c     i : lower limit of upper sub-interval (output)
c
      i=l
      j=r
      m=(l+r)/2
      x=a(m)
  220 if(a(i).ge.x) go to 230
      i=i+1
      go to 220
  230 if(a(j).le.x) go to 231
      j=j-1
      go to 230
c
  231 if(i.gt.j) go to 232
      w=a(i)
      a(i)=a(j)
      a(j)=w
      i=i+1
      j=j-1
      if(i.le.j) go to 220
c
  232 level=level+1
      if((r-i).ge.(j-l)) go to 30
      lt(level)=l
      rt(level)=j
      l=i
      go to 20
   30 lt(level)=i
      rt(level)=r
      r=j
      go to 20
   50 return
      end
