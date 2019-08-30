*cmz :  3.14/16 10/05/90  17.25.40  by  nick van eijndhoven (cern)
*-- author :
      function dpdnuc(a,b,c)
c
c called by : phpnuc
c origin    : h.fesefeldt
c
      implicit double precision (a-h,o-z)
c
c     dpdk = sqrt(a*a+(b*b-c*c)**2/(a*a) - 2.0*(b*b+c*c))/2.0
      a2 = a*a
      b2 = b*b
      c2 = c*c
      if(a2) 21,21,61
   61 continue
      arg=a2+(b2-c2)**2/a2-2.0*(b2+c2)
      if (arg) 21,21,31
   21 dpdnuc=0.0
      goto 41
   31 continue
      dpdnuc = 0.5*sqrt(a2 + (b2-c2)**2/a2 - 2.0*(b2+c2))
   41 continue
      return
      end
