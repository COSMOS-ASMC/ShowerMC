!        ***********************************************************
!        *
!        * csampCollInA: sample # of collisions inside nucleus
!        *
!        *************** tested 88.08.03***********************k.k**
!
!    Hadron nucleus collision is decompsed into successive
!         collision of incident hadron (p, pi, etc)
!         with nucleon inside the nucleus. This program
!         obtains the # of successive collisions.
!   /usage/ call csampCollInA(proj, ia, nc)
!       proj: /ptcl/   input. projectile ptcl
!         nc: output.  # of collistions sampled
!   /method/
!       Using wood-saxon density of nucleus, simplified
!         glauber calculation is done by using cwoodsaxon_den etc.
!         its results for 
!         A**(1/3) = 4**(1/3) to 208**(1/3) with step (tatal width)/15
!         and for elementary cross sections log10(15mb) to log10(300mb)
!         step (total width)/15  is tabulated. (cumProb.h)
!      Because of the descrete nature of the table, the average number
!      could be 10 % smaller than what would be obtained by
!      exact table.  The distributions for very high energy and 
!      very high mass targeti are not accurate; they will not be
!      used in actual case.
! 
        subroutine csampCollInA(pj, ia,  nc)
        implicit none

#include  "Zcode.h"
#include  "Zptcl.h"
        integer i
        type(ptcl):: pj  !input projectile particle
        integer ia       ! target mass no.
        integer nc       ! output number of collisions.

        real*8 a1, a2, a3, da
        real*8  a, xs
!
        integer  idxa, idxxs

        real*8 u
        real*8  xs1, xs2, dxs
        integer mm, nn, kk
        parameter( mm = 14, nn = 16,  kk= 16)
!                 cumProb,  xsec     A

        real*4 cumProb(mm, nn, kk)




!        parameter (a1 = 4.0**0.3333333, a2 = 208.**0.333333333333,
        parameter (a1 = 1.5874011, a2 = 5.9249921,
     *  da = (a2-a1)/(kk-1) )

!         parameter ( xs1 =log10(15.), xs2 =log10(300.) )
        parameter ( xs1 = 1.176091259, xs2 = 2.477121255, 
     *  dxs =( xs2 - xs1)/15.)

#include  "cumProb.h"
!             get cross-section for proton target.        
!        call cxpXsec(pj, xs)
        call cinelx(pj, 1.d0, 1.d0,  xs)
        xs = log10(xs)
!
         a = ia
         a3 = a**0.3333333333
         idxa = (a3- a1)/da + 1
         idxxs = (xs - xs1)/dxs + 1
         if( (a3 - idxa * da - a1) .gt. (idxa*da + da + a1 - a3) ) then
            idxa = idxa  + 1
         endif
         idxa =max(1, min(idxa, kk))

         if( (xs - idxxs * dxs - xs1) .gt.
     *        (idxxs * dxs + dxs + xs1 - xs) ) then
            idxxs = idxxs +1
         endif
         idxxs =max(1, min(idxxs, nn))
         call rndc(u)
         do   i=1,  mm
            if(u .le. cumProb(i, idxxs, idxa) ) then
               nc=i
               goto 100
            endif
         enddo
  100    continue

        end

