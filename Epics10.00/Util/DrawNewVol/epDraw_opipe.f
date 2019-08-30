      subroutine epDraw_opipe(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !   a opipe in local coordnate.
                              ! (x,y,z)= gpsep is a separator
                              ! to be converted to a blank line
                              ! dimension of p must be >+ (nvccl+2)*2
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  


       integer iir, ior, ih, isa, iea,
     *         ircossa, irsinsa, ircosea, irsinea
       parameter( iir = 9, ior=1, ih = 2,  isa=3, iea=4,
     *            ircossa=5, irsinsa=6, ircosea=7, irsinea=8  )


       real*8 r, h, sa, ea
       real*8 t1, t2
      integer n1, n2,  locoru, locorl, lociru, locirl
      logical kdgtest
      logical isinside
      logical isinside2, scut, ecut
      real*8 x
      isinside(x) = mod(thetamin-thetamax+360.d0, 360.d0) .ge.
     *               mod(x-thetamax+360.d0, 360.d0)
      isinside2(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)

       r = Volat( comp%vol + ior)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)
       scut = .false.
       ecut = .false.

       n = 0
       if(.not. isinside(sa) .and. isinside(ea) ) then
          t1 = thetamax
          t2 = ea
          ecut = .true.
       elseif(isinside(sa) .and. isinside(ea)) then
          if(isinside2(thetamin)) then
!              there is two regions
             scut = .true.
             t1 =  sa
             t2 = mod(thetamin,360.d0)
             if(t2 .lt. t1) t2 = t2+360.d0
             call epdrawCcl(r,  0.d0, t1, t2, p(n+1), n1)
             locorl = n + 1
             n = n + n1
             call epdrawCcl(r,  h, t1, t2, p(n+1), n2)
             locoru = n + 1
             n =n + n2
             n = n + 1
             p(n)%x = gpsep

             r = Volat( comp%vol + iir)
             call epdrawCcl(r,  0.d0, t1, t2, p(n+1), n1)
             locirl = n+1
             n = n + n1
             call epdrawCcl(r,  h, t1, t2, p(n+1), n2)
             lociru = n+1
             n =n + n2
             n = n + 1
             p(n)%x = gpsep

             if(kdgtest(howcyl, 1)) then
!                call epdrawPipeEdg(p(lociru), p(locoru), n1, p(n+1), n2)
                call epdrawPipeEdg(p(lociru), p(locoru), 
     *          n1-1, p(n+1), n2)
                n = n + n2
                n = n + 1
                p(n)%x = gpsep
             endif
             if(kdgtest(howcyl, 2)) then
!                call epdrawPipeEdg(p(locirl), p(locorl), n1, p(n+1), n2)
                call epdrawPipeEdg(p(locirl), p(locorl),
     *          n1-1, p(n+1), n2)
                n = n + n2
                n = n + 1
                p(n)%x = gpsep
             endif

             t1 =mod(thetamax, 360.d0)
             t2 = ea
             if(t2 .lt. t1) t2 = t2+360.d0
             r = Volat( comp%vol + ior)
             call epdrawCcl(r,  0.d0, t1, t2, p(n+1), n1)
             locorl = n+1
             n = n + n1
             call epdrawCcl(r,  h, t1, t2, p(n+1), n2)
             locoru= n+1
             n =n + n2
             n = n + 1
             p(n)%x = gpsep

             r = Volat( comp%vol + iir)
             call epdrawCcl(r,  0.d0, t1, t2, p(n+1), n1)
             locirl = n + 1
             n = n + n1
             call epdrawCcl(r,  h, t1, t2, p(n+1), n2)
             lociru= n+1
             n =n+  n2 
             n = n + 1
             p(n)%x = gpsep
             ecut = .true.

             if(kdgtest(howcyl, 1)) then
!                call epdrawPipeEdg(p(lociru), p(locoru), n1, p(n+1), n2)
                call epdrawPipeEdg(p(lociru), p(locoru),
     *          n1-1, p(n+1), n2)
                n = n + n2
                n = n + 1
                p(n)%x = gpsep
             endif
             if(kdgtest(howcyl, 2)) then
!                call epdrawPipeEdg(p(locirl), p(locorl), n1, p(n+1), n2)
                call epdrawPipeEdg(p(locirl), p(locorl),
     *          n1-1, p(n+1), n2)
                n = n + n2
                n = n + 1
                p(n)%x = gpsep
             endif
             goto 100
          else
             t1 = sa
             t2 = ea
             scut = .true.
             ecut = .true.
          endif
       elseif(isinside2(thetamax) .and. isinside2(thetamin)) then
          t1 =mod(thetamax, 360.d0)
          t2 =mod(thetamin, 360.d0)
       elseif(isinside2(thetamin) .and. isinside(sa) ) then
          t1 = sa
          t2 = mod(thetamin, 360.d0)
          scut = .true.
       else
          return  !  no part to be drawn
       endif

       if(t2 .lt. t1) t2=t2+360.d0
       call epdrawCcl(r,  0.d0, t1, t2, p(n+1), n1)
       locorl = n+1
       n = n+n1
       call epdrawCcl(r,  h, t1, t2, p(n+1), n2)
       locoru = n +1
       n = n +  n2
       n = n + 1
       p(n)%x = gpsep

       r = Volat( comp%vol + iir)
       call epdrawCcl(r,  0.d0, t1, t2, p(n+1), n1)
       locirl = n+1
       n = n +n1
       call epdrawCcl(r,  h, t1, t2, p(n+1), n2)
       lociru = n + 1
       n = n +  n2
       n = n + 1
       p(n)%x = gpsep

       if(kdgtest(howcyl, 1)) then
!         call epdrawPipeEdg(p(lociru), p(locoru), n1, p(n+1), n2)
         call epdrawPipeEdg(p(lociru), p(locoru), n1-1, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
       endif
       if(kdgtest(howcyl, 2)) then
!         call epdrawPipeEdg(p(locirl), p(locorl), n1, p(n+1), n2)
         call epdrawPipeEdg(p(locirl), p(locorl), n1-1, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
       endif
 100   continue
       if(scut) then
!           draw startng cut  square
         n = n + 1
         p(n)%x = Volat( comp%vol + ircossa)*
     *      Volat( comp%vol + iir)/Volat( comp%vol + ior)
         p(n)%y = Volat( comp%vol + irsinsa)*
     *      Volat( comp%vol + iir)/Volat( comp%vol + ior)
         p(n)%z = 0.

         n = n + 1
         p(n)%x = Volat( comp%vol + ircossa)
         p(n)%y = Volat( comp%vol + irsinsa)
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = Volat( comp%vol + ircossa)*
     *       Volat( comp%vol + iir)/Volat( comp%vol + ior)
         p(n)%y = Volat( comp%vol + irsinsa)*
     *       Volat( comp%vol + iir)/Volat( comp%vol + ior)
         p(n)%z = h

         n = n + 1
         p(n)%x = Volat( comp%vol + ircossa)
         p(n)%y = Volat( comp%vol + irsinsa)
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
      if(ecut) then
!           draw ending cut  square
         n = n + 1
         p(n)%x = Volat( comp%vol + ircosea)
         p(n)%y = Volat( comp%vol + irsinea)
         p(n)%z = 0.

         n = n + 1
         p(n)%x = Volat( comp%vol + ircosea)*
     *     Volat( comp%vol + iir)/Volat( comp%vol + ior)
         p(n)%y = Volat( comp%vol + irsinea)*
     *     Volat( comp%vol + iir)/Volat( comp%vol + ior)
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = Volat( comp%vol + ircosea)
         p(n)%y = Volat( comp%vol + irsinea)
         p(n)%z = h

         n = n + 1
         p(n)%x = Volat( comp%vol + ircosea)*
     *    Volat( comp%vol + iir)/Volat( comp%vol + ior)
         p(n)%y = Volat( comp%vol + irsinea)*
     *    Volat( comp%vol + iir)/Volat( comp%vol + ior)
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
      end
