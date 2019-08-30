      subroutine epDraw_gbox(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)      ! output. (x,y,z) to describe
                               !  gbox in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      logical kdgtest
 
       type(epPos)::  o, a, b, c, ap, bp, cp, op
       integer iax,  ibx, iby, ibz, icx, icy, icz, 
     *         iapx, iapy, ibpx, ibpy, ibpz, icpx, icpy, icpz,
     *         iopx, iopy
       parameter(iax=1, ibx=2, iby=3, ibz=4, icx=5, 
     *           icy=6, icz=7, iapx=8, iapy=9, ibpx=10, 
     *           ibpy=11, ibpz=12, icpx=13, icpy=14, icpz=15,
     *           iopx=16, iopy=17)


        a%x = Volat( comp%vol + iax)
        a%y = 0.
        a%z = 0.

        b%x = Volat( comp%vol + ibx)
        b%y = Volat( comp%vol + iby)
        b%z = Volat( comp%vol + ibz)

        c%x = Volat( comp%vol + icx)
        c%y = Volat( comp%vol + icy)
        c%z = Volat( comp%vol + icz)

        ap%x = Volat( comp%vol + iapx)
        ap%y = Volat( comp%vol + iapy)
        ap%z = 0. 

        bp%x = Volat( comp%vol + ibpx)
        bp%y = Volat( comp%vol + ibpy)
        bp%z = Volat( comp%vol + ibpz)

        cp%x = Volat( comp%vol + icpx)
        cp%y = Volat( comp%vol + icpy)
        cp%z = Volat( comp%vol + icpz)

        op%x = Volat( comp%vol + iopx)
        op%y = Volat( comp%vol + iopy)
        op%z = 0.
        o%x = 0.
        o%y = 0.
        o%z = 0. 



 
!          we follow the  logic used in box drawing
      if( kdgtest(how, 1) )then
         n = n + 1
         p(n)= o
         n = n + 1
         p(n) = a
         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n) = op
         n = n + 1
         p(n) = ap

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
     
      if( kdgtest(how, 2) )then
         n = n + 1
         p(n) = o
         n = n + 1
         p(n) = a
         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n) = c
         n = n + 1
         p(n) = b
         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      if( kdgtest(how, 3) )then
         n = n + 1
         p(n) = o
         n = n + 1
         p(n) = op
         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n) = c
         n = n + 1
         p(n) = cp
         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      if( kdgtest(how, 4) )then
         n = n + 1
         p(n) = a
         n = n + 1
         p(n) = ap
         n = n + 1
         p(n)%x = gpsep
!        -----------
         n = n + 1
         p(n) = b
         n = n + 1
         p(n) = bp
         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif


      if( kdgtest(how, 5) )then
         n = n + 1
         p(n) = op
         n = n + 1
         p(n) = ap
         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n) = cp
         n = n + 1
         p(n) = bp
         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif


      if( kdgtest(how, 6) )then
         n = n + 1
         p(n) = c
         n = n + 1
         p(n) = b
         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n) = cp
         n = n + 1
         p(n) = bp
         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
      end
