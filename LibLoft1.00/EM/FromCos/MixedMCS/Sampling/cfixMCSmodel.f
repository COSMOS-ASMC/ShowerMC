      subroutine cfixMCSmodel( aPtcl )
      use modMCScontrol
      implicit none
!
!        fix how to treat MCS  for this e-/e+
!
#include "Zptcl.h"
#include "Zcode.h"

      type(ptcl):: aPtcl ! input. e+/e-

      integer i
      integer :: icon

      if( aPtcl%code /= kelec) then
         write(0,*) 'cfixMSmodel is called non e-/e+'
         write(0,*) ' the code,subcode,charge= ', 
     *    aPtcl%code, aPtcl%subcode, aPtcl%charge
         stop 
      endif
      KEGeV = aPtcl%fm%p(4)- aPtcl%mass
      KEeV = KEGeV*1.0d9 
      if(KEGeV > MaxErgMCS) then
         ActiveMCS = 'Mol'
      else
         do i = 1, NoOfHowMCS
            if( KEGeV <  MCSErg(i)) then
               if(KEGeV > MinErgMCS) then
                  ActiveMCS = HowMCSList(i)
               else
                  ActiveMCS = 'Mol'
               endif
!!!!!!!!!!!!!
!!               write(0,*) 'E=',KEeV,' model=',ActiveMCS
!!!!!!!!!!!
               return           !!!!!!!!!!
            endif
         enddo
         call cerrorMsg('MCSModel and energy range  strange', 1)
         do i = 1,  NoOfHowMCS
            write(0,*) HowMCSList(i), MCSErg(i)
         enddo
         stop 44444
      endif
      end
