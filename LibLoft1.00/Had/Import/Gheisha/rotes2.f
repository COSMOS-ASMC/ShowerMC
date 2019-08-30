*cmz :  3.14/16 23/05/89  10.04.05  by  nick van eijndhoven (cern)
*-- author :
      subroutine rotes2(c,s,c2,s2,pr,i)
c
c *** nve 16-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (27-oct-1983)
c
      dimension pr(*)
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
