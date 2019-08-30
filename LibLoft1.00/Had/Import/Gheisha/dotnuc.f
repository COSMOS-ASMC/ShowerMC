*cmz :  3.14/16 10/05/90  17.25.40  by  nick van eijndhoven (cern)
*-- author :
      subroutine dotnuc(c,s,c2,s2,pr,i)
c
      implicit double precision (a-h,o-z)
      dimension pr(50)
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
