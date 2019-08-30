#include "ZsaveStruc.h"
!       **************************************************************
!       *
!       * cibst1:  boost a partilce into the moving system of another
!       *         particle.
!       *
!       **************************************************************
!
! /usage/  call cibst1(init, p1, p2,  po)
!
!        Suppose two particles p1, p2.
!        p2 is given in the rest system of p1 which is given in a system
!        K.
!        This boosts particle 2 into  the K system.
!        (x,y,z) axis are assumed to be parallel to that of the K.
!        If the z axis at the rest system of p1 is the direction of
!        p1, this must not be used. (Such case is muon decay where
!        polarization exits to the z-direction). use cibst1Pol
!       
!  init: integer.  Input.  If p1 is the same as the previous call to this
!        subroutine, give a value other than 1.  If p1 is different
!        from the previous call, use 1. 
!    p1: type ptcl. Input. particle 1
!    p2: type ptcl. Input. particle 2
!    po: type ptcl. Output. particle seen at K
!        po may be  the same one as p2
!
       subroutine cibst1(init, p1, p2, po)
         implicit none
!----         include '../Zptcl.h'
#include  "Zptcl.h"
         type(ptcl):: p1, p2, po
         integer init
!
         type(fmom):: g
#ifdef  USESAVE         
         save g
#endif
!
         if(init .eq. 1) then
            call cgetlf(p1, g)
         endif   
         po = p2

         call clorep(init, g, p2, po)

       end
#include "ZsaveStruc.h"
!       **************************************************************
!       *
!       * cibstPol:  boost a partilce into the moving system of another
!       *         particle.
!       *
!       **************************************************************
!
! /usage/  call cibstPol(init, p1, p2,  po)
!
!        Suppose two particles p1, p2.
!        p2 is given in the rest system of p1 which is given in a system
!        K.
!        This boosts particle 2 into  the K system.
!        z axis at the rest system of p1 is assumed to be the direction
!        of p1 in K.  This may be used for muon decay where polaization
!        exits.
!
!  init: integer.  Input.  If p1 is the same as the previous call to this
!        subroutine, give a value other than 1.  If p1 is different
!        from the previous call, use 1. 
!    p1: type ptcl. Input. particle 1
!    p2: type ptcl. Input. particle 2
!    po: type ptcl. Output. particle seen at K
!        po may be  the same one as p2
!
       subroutine cibstPol(init, p1, p2, po)
         implicit none
!----         include '../Zptcl.h'
#include  "Zptcl.h"
         type(ptcl):: p1, p2, po
         integer init
!
         type(fmom):: g
#ifdef  USESAVE         
         save g
#endif
!
         if(init .eq. 1) then
            call cgetlf(p1, g)
         endif   
         po = p2
         call cloreb(init, g, p2, po)
       end


