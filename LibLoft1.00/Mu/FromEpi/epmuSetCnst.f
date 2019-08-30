!     *********************************
      subroutine epmuSetCnst(zin,  ain)
!     *********************************
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmuint.h"
#include "Zmass.h"

       real*8 zin  ! charge of the target atom (sigle atom)
       real*8 ain  ! mass number of the target (//)

       real*8 MubyMe
       parameter (MubyMe = masmu/masele)

       A =  ain
       A3 = A**(1d0/3.d0)
       Gzai = 0.25d0   !  standard value but may increase at high Emu
       Z = zin
       Z2 = Z*Z
       if(Z  .gt. 10.) then
          Zp2 = Z*(Z + 1.2)
       else
          Zp2 = Z*(Z +  1.)
       endif
!           z(z+1.2)*alpha*(2r0 Me/Mu)**2
       D =  Zp2 * ar02 *(2* masele/masmu)**2  !  in mb
       LogZ = log(Z)
       Z3 = Z**(1.d0/3.d0)
       Shadow = A**(-0.1d0)
       PointLike=0.22d0
       Ak = 189.d0
       Akm =   Ak * MubyMe/Z3
       Akm2  = Akm * sqrt(exp(1.d0))
       end




