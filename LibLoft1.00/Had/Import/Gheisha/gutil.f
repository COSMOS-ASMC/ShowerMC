c        //// modified by k.k to use rndc.
      subroutine grndm(rvec,len)
      implicit none
       real*8 u
       integer *4 len,i
c       real *4 rvec(len),rando,dum
       real *4 rvec(len)
        do 1000 i=1,len
           call rndc(u)
           rvec(i)= u
1000    continue
      end
c
      subroutine granor(a,b)
c.
c.    ******************************************************************
c.    *                                                                *
c.    *       to generate 2 numbers a and b following a normal         *
c.    *       distribution (mean=0 sigma=1.)                           *
c.    *         copy of the cern library routine rannor                *
c.    *                                                                *
c.    *    ==>called by : <user>, many geant routines                  *
c.    *       author    f.carminati *********                          *
c.    *                                                                *
c.    ******************************************************************
c.
      dimension rndm(2)
*
      call grndm(rndm,2)
      y=rndm(1)
      z=rndm(2)
      x=6.283185*z
      a1=sqrt (-2.0*log(y))
      a=a1*sin (x)
      b=a1*cos (x)
      return
      end
      subroutine dummy
      entry uga013
      entry uga014
      return
      end

      subroutine vzero(array, n)
      real *4 array(n)
      
      if(n.gt.0) then
        do 10 i=1,n
           array(i) = 0
  10    continue
      endif
      return
      end

      subroutine ucopy(array1, array2, n)
      real *4 array1(n), array2(n)
      
      if(n.gt.0) then
        do 10 i=1,n
           array2(i) = array1(i)
  10    continue
      endif
      return
      end

      subroutine flpsor(a,n)
c
c ***  cernlib routine "flpsor" ***
c *** nve 29-mar-1988 cern geneva ***
c
c called by : phasp
c origin : h.fesefeldt (02-dec-1986)
c
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

