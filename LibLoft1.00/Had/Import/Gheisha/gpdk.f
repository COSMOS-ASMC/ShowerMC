*cmz :  3.14/16 29/06/89  11.35.32  by  nick van eijndhoven (cern)
*-- author :
      function gpdk(a,b,c)
c
c *** nve 16-mar-1988 cern geneva ***
c
c called by : phasp
c origin : h.fesefeldt (27-oct-1983)
c
c     gpdk = sqrt(a*a+(b*b-c*c)**2/(a*a) - 2.0*(b*b+c*c))/2.0
c
      a2 = a*a
      b2 = b*b
      c2 = c*c
      if(a2) 21,21,61
   61 continue
      arg=a2+(b2-c2)**2/a2-2.0*(b2+c2)
      if (arg) 21,21,31
   21 gpdk=0.0
      goto 41
   31 continue
      gpdk = 0.5*sqrt(abs(a2 + (b2-c2)**2/a2 - 2.0*(b2+c2)))
   41 continue
      return
      end
