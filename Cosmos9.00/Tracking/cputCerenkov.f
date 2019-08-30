!         ready made output for cerenkov 
!         threshold dependence on the air density  is
!         taken into account.
!
       subroutine cputCerenkov
       implicit none
#include "Zmanagerp.h"
#include "Ztrack.h"
! #include "Zmagfield.h"       
#include "Ztrackp.h"
#include "Ztrackv.h"
#include "Zcode.h"
#include "Zmass.h"
#include "Zprimary.h"
#include "Zprimaryv.h"
#include "Zheavyp.h"
#include "Zincidentv.h"
!
      real*8 crnth(knuc)
!
!  
       real*8  thg,  thpi, thk, thmu, thnuc
       parameter (thg=40.-1., thpi=thg*maspic, thk=thg*maskc,
     1  thmu=thg*masmu, thnuc=thg*masp)
!                threshold energy in GeV at sea level ???
!
        type(coord)::f
        type(coord)::t
        integer chrg, code
!
        real*8 h1, den1, e1, the
!        real*8 h2, den2, e2
        real*8 cvh2den
        logical first/.true./
        real*8 denAtSea
        save first, denAtSea
        integer ka, itb, it,  utrace
        data crnth(kelec)/20.e-6/, crnth(kpion)/thpi/,
     1  crnth(kkaon)/thk/, crnth(kmuon)/thmu/,
     2  crnth(knuc)/thnuc/
        
        if(first) then
           denAtSea = cvh2den(0.d0)
           first = .false.
        endif

        utrace = Trace
        if(utrace .gt. 160) utrace = utrace - 100
!          get air density.
        h1 = TrackBefMove%pos%height
        den1 = cvh2den(h1)
!        h2 = MovedTrack.pos.height
!        den2 = cvh2den(h2)
        ka = TrackBefMove%p%code
        e1 = TrackBefMove%p%fm%p(4)
!        e2 = MovedTrack.p.fm.p(4)
!
!        if( ka .gt. knuc) then
!        the = thg * TrackBefMove%p%mass
        the = 58.4* TrackBefMove%p%mass
!        else
!           the = crnth(ka)
!        endif
!
        if( e1 .gt. the*denAtSea/den1 ) then
!               get transformed coord.
!  utrace:      61-70:  --> 1ry
!               71-80:  --> 1ry but z is depth above 
!               81-90:  --> det 
!               91-100 -->  det but z is depth above
!
           call ccoordForTr(utrace-60, f, t)
           itb=TrackBefMove%t*100  ! centimeter/beta
           it=(MovedTrack%t - TrackBefMove%t)*100  !  length/beta (cm)
           chrg = TrackBefMove%p%charge

           if(Trace .gt. 160) then
              call chookCeren(ka, chrg, e1, itb, it, f, t)
           else
              if(mod(utrace, 2) .eq. 0 ) then
                 write(TraceDev) ka, chrg, sngl(e1), itb, it, 
     *           sngl(f%r(1)), sngl(f%r(2)), sngl(f%r(3)),
     *           sngl(t%r(1)), sngl(t%r(2)), sngl(t%r(3))
              else
                 write(TraceDev, *) ka, chrg, sngl(e1), itb, it, 
     *           sngl(f%r(1)), sngl(f%r(2)), sngl(f%r(3)),
     *           sngl(t%r(1)), sngl(t%r(2)), sngl(t%r(3))
              endif
           endif 
       endif
       return
!        *************  start of 1 shower for cerenkov trace(called from
!                      ciniTracking
!
       entry cputCerenkovS
       utrace = Trace
       if(utrace .gt. 160) utrace = utrace -100

       if(Trace .gt. 160) then
          call chookCerenS(EventNo,  Prim, AngleAtObsCopy)
       else
          code = Prim%particle%code
          if(mod(utrace, 2) .eq. 0) then
             write(TraceDev) EventNo,
     *       code, 
     *       Prim%particle%fm%p(4), 
     *       AngleAtObsCopy%r(1), AngleAtObsCopy%r(2),
     *       AngleAtObsCopy%r(3)
          else
             write(TraceDev, *) EventNo,
     *      code, 
     *      Prim%particle%fm%p(4), 
     *      AngleAtObsCopy%r(1), AngleAtObsCopy%r(2),
     *      AngleAtObsCopy%r(3)
          endif
       endif
       return
!        ************** end of  1 shower
         entry  cputCerenkovE
!        ***************
         utrace = Trace
         if(utrace .gt. 160) utrace = utrace - 100
         ka = 0
         chrg = 0
         e1 = 0.
         itb =0
         it =0
         f%r(1) = 0.
         f%r(2) = 0.
         f%r(3) =0
         t%r(1) = 0.
         t%r(2) = 0.
         t%r(3) = 0.
         if(Trace .gt. 160) then
            call chookCerenE(ka, chrg, e1, itb, it, f, t)
         else
!
            if(mod(utrace, 2)  .eq. 0) then
               write(TraceDev) ka, chrg, sngl(e1), itb, it,
     *         sngl(f%r(1)), sngl(f%r(2)), sngl(f%r(3)),
     *         sngl(t%r(1)), sngl(t%r(2)), sngl(t%r(3))
            else 
               write(TraceDev, *) ka, chrg, sngl(e1), itb, it,
     *         sngl(f%r(1)), sngl(f%r(2)), sngl(f%r(3)),
     *         sngl(t%r(1)), sngl(t%r(2)), sngl(t%r(3))
            endif
         endif
         end
