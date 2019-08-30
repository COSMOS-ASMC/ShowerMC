*cmz :  3.14/16 13/03/89  14.48.42  by  nick van eijndhoven (cern)
*-- author :
c
c *** various vector operations ***
c
c copied from f14blo.pamlib 23.4.82
c author: v.blobel (university of hamburg)
c desylib
c
c *** blank common replaced by /vecuty/ to match geant/gheisha code ***
c *** note that p(10,100) has become pv(10,200) due to this         ***
c
c un-used entries removed :
c     "pcop" "pexc" "pzer" "pwrt" "dot4" "impu" "impuli" "add3"
c     "sub3" "cross" "dot" "smul" "norz" "parper" "punit" "trap"
c
c *** all entries re-written as subroutines using only necessary ***
c *** "double precision" stmts. and all specific functions have  ***
c *** been changed to their generic equivalences                 ***
c *** nve 29-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (22-june-1984)
c
      subroutine add(k,l,m)
c
c --- pv-array ---
      common /vecuty/ pv(10,200)
c
c
      double precision a,b
c
      a=pv(4,k)+pv(4,l)
      pv(4,m)=a
      b=a*a
      do 2 i=1,3
      a=pv(i,k)+pv(i,l)
      b=b-a*a
      pv(i,m)=a
 2    continue
      pv(5,m)=sign(sqrt(abs(b)),b)
      do 3 i=9,10
      pv(i,m)=pv(i,k)+pv(i,l)
 3    continue
      return
      end
