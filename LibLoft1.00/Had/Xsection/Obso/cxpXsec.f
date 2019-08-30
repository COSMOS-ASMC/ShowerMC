!        ***********************************************************
!        *
!        * cxpXsec:  get xsection for proton  target.
!        *
!
        subroutine cxpXsec(pj, xs)
        implicit none

#include  "Zcode.h"
#include  "Zptcl.h"

        type(ptcl):: pj
        real*8  xs, ke
        integer k
!

        character*100 msg



         k=pj.code
         ke = pj.fm.p(4) - pj.mass
         if(ke .le. 0.) then
            xs = 0.
            return   !*******
         endif
!

         if(k .eq. kpion .or. k .eq. keta) then
            if(pj.charge .ge. 0) then
               call cpiPluspXsec(ke, xs)
            else
               call cpiMinuspXsec(ke, xs)
            endif
         elseif(k .eq. kkaon ) then
            if(pj.charge .ge. 0.) then
               call ckPluspXsec(ke, xs)
            else
               call ckMinuspXsec(ke, xs)
            endif
         elseif(k .eq. knuc) then
            if(pj.charge .ge. 0.) then
               call cppXsec(ke, xs)
            else
               call cpbarpXsec(ke, xs)
            endif
         elseif(k .eq. krho  .or. k .eq. komega  
     *         .or. k .eq. kphi  .or. k .eq. kdmes) then
            if(pj.charge .ge. 0) then
               call cpiPluspXsec(ke, xs)
            else
               call cpiMinuspXsec(ke, xs)
            endif
         elseif( k .ge. klambda .and. k .le. kbomega) then
            call cppXsec(ke, xs)
         else
               write(msg,*) ' undef.ptcl for cxpXsec=', k,
     *         ' K.E=', ke
               call cerrorMsg(msg, 0)
         endif
        end
