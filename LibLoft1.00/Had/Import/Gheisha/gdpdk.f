*cmz :  3.14/16 29/06/89  11.36.18  by  nick van eijndhoven (cern)
*-- author :
      function gdpdk(a,b,c)
c
c *** double precision version of "pdk" ***
c *** nve 29-mar-1988 cern geneva ***
c
c called by : phasp
c origin : h.fesefeldt (02-dec-1986)
c
      implicit double precision (a-h,o-z)
c
c     gdpdk = sqrt(a*a+(b*b-c*c)**2/(a*a) - 2.0*(b*b+c*c))/2.0
c
      a2 = a*a
      b2 = b*b
      c2 = c*c
      if(a2) 21,21,61
   61 continue
      arg=a2+(b2-c2)**2/a2-2.0*(b2+c2)
      if (arg) 21,21,31
   21 gdpdk=0.0
      goto 41
   31 continue
      gdpdk = 0.5*sqrt(abs(a2 + (b2-c2)**2/a2 - 2.0*(b2+c2)))
   41 continue
      return
      end
