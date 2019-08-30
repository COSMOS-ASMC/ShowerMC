*cmz :  3.14/16 13/03/89  14.48.40  by  nick van eijndhoven (cern)
*-- author :
      integer function nfac(n)
c
c *** nve 16-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (27-oct-1983)
c
      nfac=1.
      m=n
      if(m.le.1) return
      if(m.gt.10) m=10
      do 1 i=2,m
    1 nfac=nfac*i
      return
      end
