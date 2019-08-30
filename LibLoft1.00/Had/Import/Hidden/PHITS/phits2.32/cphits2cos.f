      subroutine cphits2cos(ityp, ktyp, code, subcode, charge)
!    this is probably not used. (but next sub is used)
      implicit none
#include "Zcode.h"
      integer,intent(in)::ityp   !  phits code
      integer,intent(in)::ktyp   !  KF code (but A>1, diff)
      integer,intent(out)::code  !  cosmos  code
      integer,intent(out)::subcode  !  cosmos  subcode      
      integer,intent(out)::charge   !  cosmos  charge
      if( ityp >= 15 ) then
         code = 9   ! = kgnuc
         charge =  ktyp/1000000
         subcode =ktyp- charge*1000000
         if( subcode <= 1 ) then  
!               proton may  come out as  1000001
            code = knuc
            subcode = regptcl
         endif   
      else
         call ckf2cos(ktyp, code, subcode, charge)
      endif
*               ityp :  description
*                 1  :  proton 
*                 2  :  neutron                                                
*                 3  :  pion (+)                                               
*                 4  :  pion (0)                                               
*                 5  :  pion (-)                                               
*                 6  :  muon (+)                                               
*                 7  :  muon (-)                                               
*                 8  :  kaon (+)                                               
*                 9  :  kaon (0)                                               
*                10  :  kaon (-)
*                11  :  other ( neu etc; by kf code)
*                12  :  electron                                               
*                13  :  positron
*                14  :  photon
*                15  :  deuteron
*                16  :  triton                                                 
*                17  :  3He                                                    
*                18  :  Alpha                                                  
*                19  :  Neucleus
*                20  :  All
      end
      subroutine ccos2phits(code, subcode, charge,ityp, ktyp)
      implicit none
      integer,intent(in)::code  !  cosmos  code
      integer,intent(in)::subcode  !  cosmos  subcode      
      integer,intent(in)::charge   !  cosmos  charge
      integer,intent(out)::ityp   !  phits code
      integer,intent(out)::ktyp   !  KF code
      integer kftp
      if( code /= 9 ) then
         call ccos2kf(code, subcode, charge, ktyp)
         ityp = kftp(ktyp)
      else
         ktyp = charge*1000000 + subcode
         if( charge == 1 ) then
            if( subcode == 2 ) then
               ityp = 15  ! deuteron
            else
               ityp = 16  ! triton
            endif
         elseif(charge == 2) then
            if(subcode == 3 ) then
               ityp =  17   !  3He
            else
               ityp = 18   ! 4He
            endif
         else
            ityp = 19
         endif
      endif 

      end
