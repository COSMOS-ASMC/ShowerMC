*cmz :  3.14/16 13/03/89  14.48.44  by  nick van eijndhoven (cern)
*-- author :
      subroutine dotes2(c,s,c2,s2,pr,i)
c
c *** double precision version of "rotes2" ***
c *** nve 29-mar-1988 cern geneva ***
c
c called by : phasp
c origin : h.fesefeldt (02-dec-1986)
c
      implicit double precision (a-h,o-z)
c
      dimension pr(50)
c
      k1 = 5*i - 4
      k2 = k1 + 1
      sa = pr(k1)
      sb = pr(k2)
      a      = sa*c - sb*s
      pr(k2) = sa*s + sb*c
      k2 = k2 + 1
      b = pr(k2)
      pr(k1) = a*c2 - b*s2
      pr(k2) = a*s2 + b*c2
      return
      end
