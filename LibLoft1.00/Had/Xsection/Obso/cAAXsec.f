!
      subroutine cAAXsec(pjmassn, tgmassn, ek,  xs)
! this is now obsolute
!         pjmassn:  input. integer. incident heavy mass number.
!         tgmassn:  input. real*8. average target mass number.
!              ek:  input. real*8. incident heavy kinetic energy in GeV.
!              xs:   output. real*8. collision x-section in mb.
!
!     formula: xs = pi(r1 + r2 - d)**2
!        r = r0 A**(1/3)
!        d = 1.189 r0 exp(-0.055min(A1,A2))
!       r0 = (1.29f to 1.41f)
!  Therefore
!      xs = pi r0^2 (A1^0.333 + A2^0.333 - 1.189exp(-0.055min(A1,A2)))**2
!  This is by Frier et al in ICRC Paris conf. (from Uchu Hosha Sen Edited
!      by Nishimura, p.170)
!
!  At high energies (sqrt(s) > 80 GeV for E/A), we include energy dependence
!  of cross-section as follows.
!    Let the pp cross-section increases as E**delta, then
!    AB crosssection is well fitted by the dependence of E**alfa with
!    
!     alfa = 2.5* delta/(p + t + p*t/2)
!
!   where  p = A**(1/3) and t = B**(1/3)
!
!  
!       

      implicit none
#include "Zxsectionp.h"
      integer pjmassn
      real*8  tgmassn, ek, xs
!      character*70 msg
!      ------------------------old ------(main frame age)
!      if(abs(tgmassn-14.5) .lt. 5.) then
!         xs =  45.2 * (pjmassn**.333 + 2.03)**2
!      else
!         write(msg, *) ' update cAAXsection so that target can be ',
!     *               ' non air'
!         call cerrorMsg(msg, 0)
!      endif    
!        
!      ---------------------------------------------------------
!         difference for He-Air(=14.5)) collisions ( r0=1.29 )
!      old         new 
! He   591 mb      492   mb
! Fe  1549        1710 
!
       real*8 p, t
       real*8 einc/500./
       write(0,*) ' use cAAxsec2 instead of obsolete cAAxsec; sorry' 
       stop

       if( ek .le. 0. ) then
          xs = 0.
       else
          p = pjmassn**0.3333
          t = tgmassn**0.3333 
          xs = 52.2 *( p + t - 
     *       1.189 * exp(- 0.055*min( dble(pjmassn), tgmassn)))**2
          if(ek/pjmassn .gt. einc) then
             xs = xs * (ek/einc/pjmassn)**(2.5* Deltpp/(p+t+p*t/2.))
          endif 
       endif
      end

