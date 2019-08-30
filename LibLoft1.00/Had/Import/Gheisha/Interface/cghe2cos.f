      subroutine cghe2cos(ghecode, code, subcode, charge)
      implicit none
!         convert gheisha particle code into cosmos codes 
!  **  actually, the input code is Giant code. **
!     Giant code is first converted to gheisha code, and
!     it is converted to cosmos code.

      integer ghecode ! input.  gheisha particle code
      integer code    ! output. cosmos ptcl code
      integer subcode ! output. cosmos ptcl subcode
      integer charge  ! output. cosmos ptcl charge
#include "Zcode.h"

!
!      for particls to be neglected in cosmos, code will be  0
!     (Omega, tau lepton, DS) no neutrionos are expected  come
!      out from ghesha.
!
      integer ghe2cosCode(48), ghe2cosSubC(48), ghe2cosChg(48)

      character*80 msg

!

      	data ghe2cosCode/
     * kphoton, kelec, kelec, krare, kmuon, kmuon, kpion,
     * kpion, kpion, kkaon, kkaon, kkaon, knuc, knuc,
     * knuc, kkaon, keta, klambda, ksigma, ksigma,
     * ksigma, kgzai, kgzai, kbomega, knuc, klambda,
     * ksigma, ksigma, ksigma, kgzai, kgzai, kbomega, krare,
     * krare, kdmes, kdmes, kdmes, kdmes, krare, krare, 
     * klambdac,
     * krare, krare, krare, krare, ktriton, kalfa, krare/

       	data ghe2cosSubC/
     * 0, antip, regptcl, krare, 0, 0, 0, 0, 0, k0l, 0,
     * 0, regptcl, regptcl, antip, k0s, 0, regptcl, 
     * regptcl, regptcl, regptcl, regptcl, regptcl,
     * regptcl, antip, antip, antip, antip, antip, antip,
     * antip, antip, krare, krare, 0, 0, regptcl, antip, 
     * krare, 
     * krare, regptcl, krare, krare, krare, 0, 0, 0, krare/

         data ghe2cosChg/
     * 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 
     * 0, 0, 0, 1, 0, -1, 0, -1, -1, 0, 0, -1, 0, 1, 0, 
     * 1, 1, 1, -1, 1, -1, 0, 0, 1, -1, 1, 1, -1, 0, 1, 
     * 1, 2, 0/


       if(ghecode .ge. 1 .and. ghecode .le. 48) then
          if(ghecode .eq. 45 ) then
!              deuteron
             code = kgnuc
             subcode = 2
             charge = 1
          else
             code = ghe2cosCode(ghecode)
             subcode = ghe2cosSubC(ghecode)
             charge = ghe2cosChg(ghecode)
          endif
       elseif(ghecode .eq. 200) then
!            psudo ptcl by Gheisha, I don't know this. neglect
          code = krare
       else
          write(msg, *)'ghesha (giant) code=', ghecode, 'invalid'
          call cerrorMsg(msg, 0)
       endif
       end




