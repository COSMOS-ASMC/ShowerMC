      subroutine cgeqm(p1, p2, q, icon)
!         get equivalent mass and 4 momentum of two particles
!      p1: type ptcl. Input.  particle 1
!      p2: type ptcl. Input.  particle 2
!      q:  type ptcl. Output. equivalent particle 4 momentum
!                                    and mass are given
!   icon: integer.    Output.  if m^2 <0, icon=1, else icon=0
!                              if icon =1, m becomes the m^2.
!
         implicit none
!----         include '../Zptcl.h'
#include  "Zptcl.h"
         type(ptcl):: p1, p2
         type(ptcl):: q
         integer icon
!
!            equivalent mass
        q%mass = p1%mass**2 + p2%mass**2 + 2* (p1%fm%p(4) * p2%fm%p(4)
     *   - (p1%fm%p(1) * p2%fm%p(1) + p1%fm%p(2) * p2%fm%p(2) + 
     *      p1%fm%p(3) * p2%fm%p(3)))
        if(q%mass .ge. 0.d0) then
             q%mass = sqrt(q%mass)
             icon = 0
         else
             icon = 1
         endif
!    
         q%fm%p(1) = p1%fm%p(1) + p2%fm%p(1)
         q%fm%p(2) = p1%fm%p(2) + p2%fm%p(2)
         q%fm%p(3) = p1%fm%p(3) + p2%fm%p(3)
         q%fm%p(4) = p1%fm%p(4) + p2%fm%p(4)
        end         
      subroutine cgeqm2(p1, p2, q, icon)
      !   if p1/p2 is  nucleus, /n energy is used
      !   to get q
#include  "Zptcl.h"
#include  "Zcode.h"
         type(ptcl):: p1, p2
         type(ptcl):: q
         integer icon

         type(ptcl):: pp1, pp2
         if( p1%code == kgnuc ) then
            call cmkptc(6, -1, 1, pp1)
            pp1%fm%p(1:4) = p1%fm%p(1:4)/p1%subcode
         else
            pp1 = p1
         endif


         if( p2%code == kgnuc ) then
            call cmkptc(6,-1,1, pp2)
            pp2%fm%p(1:4) = p2%fm%p(1:4)/p2%subcode
         else
            pp2 = p2
         endif

         call cgeqm(pp1, pp2, q, icon)
         end

      

