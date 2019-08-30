#include "ZsaveStruc.h"
!       **************************************************************
!       *
!       * cbst0:  boost a partilce into the rest system of another
!       *         particle
!       *
!       **************************************************************
!
! /usage/  call cbst0(init, gb, p,  po)
!
!        Suppose two particles given in the system K. 
!        The Lorentz factor of particle 1 is given in K.
!        4 momentum of particle 2 is also given in K.  This boosts
!        particle 2 into  the rest frame of particle 1.

!  init: integer.  Input.  If gb is the same as the previous call to this
!        subroutine, give a value other than 1.  If gb is different
!        from the previous call, use 1. 
!    gb: type fmom. Input. Lorentz factor of particle 1 in K.
!    p:  type ptcl. Input. particle 2
!   po:  type ptcl. Output. particle seen at the rest system
!        of particle 1.
!
       subroutine cbst0(init, gb, p, po)
         implicit none

#include  "Zptcl.h"
         type(fmom):: gb
         type(ptcl):: p, po
         integer init
!
         type(fmom):: ig
#ifdef  USESAVE
         save ig
#endif
!

         if(init .eq. 1) then
            ig%p(1) = -gb%p(1)
            ig%p(2) = -gb%p(2)
            ig%p(3) = -gb%p(3)
            ig%p(4) =  gb%p(4)
         endif
         po = p
         call clorep(init, ig, p, po)
       end


