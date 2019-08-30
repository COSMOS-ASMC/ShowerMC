      subroutine cconsvChg(outc, a, ntp, icon)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"

      integer outc ! input. outc + sum of charge of a(i)...a(ntp) = 0
      integer ntp  ! input. number of paticles produced except for proj and tgt
      type(ptcl):: a(ntp)  ! input/output.  partcles (except for proj and tgt)
      integer icon ! output. 0--> o%k. 1--> n%g

!
!      Total charge of 'a' is adjuested to be the same as
!      rproj.charge + rtgt.charge
!      
      integer sum, i,  ncc, nnn
      real*8 u
      sum = 0.
      do i = 1, ntp
         sum = sum + a(i)%charge
      enddo
      sum = sum + outc
      ncc = 0

      do while( sum .ne. 0 .and. ncc .lt. 10)
         nnn = 0
         do while (abs(sum) .ge. 2  .and. nnn .lt. 10)
!         very rare to have nnn > 10
!
!            charge 'sum' must be made to be 0
!
            call rndc(u)
            if(sum .gt. 0) then
               if(u .lt. 0.3333d0) then
                  call cchgopposit(1, a, ntp,icon)
                  if(icon .eq. 0) then
                     sum = sum - 2
                  endif
!               elseif(u .lt. 0.6666d0) then
!                  call cchg0(1, a, ntp, icon)
!                  if(icon .eq. 0) then
!                     sum = sum - 1
!                  endif
               else
!                  call rndc(u)
                  call cchg0c(-1, a, ntp, icon)
                  if(icon .eq. 0) then
                     sum = sum - 1
                     call cchg0c(-1, a, ntp, icon)
                     if(icon .eq. 0) then
                        sum= sum-1
                     endif
                  endif
!                  elseif(u .lt. 0.1) then
!                     call cchg0(1, a, ntp, icon)
!                     if(icon .eq. 0) then
!                        sum = sum -1
!                     endif
!                  else
!                     call cchgopposit(1, a, ntp, icon)
!                     if(icon  .eq. 0) then
!                        sum = sum - 2
!                     endif
!                  endif
               endif
            else
               if(u .lt. 0.3333d0) then
                  call cchgopposit(-1, a, ntp, icon)
                  if(icon .eq. 0) then
                     sum = sum + 2
                  endif
!               elseif(u .lt. 0.666d0) then
!               elseif(u .lt. 0.95d0) then
!                  call cchg0(-1, a, ntp, icon)
!                  if(icon .eq. 0) then
!                     sum = sum + 1
!                  endif
               else
                  call rndc(u)
                  call cchg0c(1, a, ntp, icon)
                  if(icon .eq. 0) then
                     sum = sum + 1
                     call cchg0c(1, a, ntp, icon)
                     if(icon .eq. 0) then
                        sum = sum + 1
                     endif
                  endif
!                  elseif(u .lt. 0.8d0) then
!                     call cchg0(-1, a, ntp, icon)
!                     if(icon .eq. 0) then
!                        sum = sum + 1
!                     endif
!                  else
!                     call cchgopposit(-1, a, ntp, icon)
!                     if(icon .eq. 0) then
!                        sum = sum + 2
!                     endif
!                  endif
               endif
            endif
            nnn = nnn + 1 
         enddo
         if(sum .eq. 1) then
            call cchg0(1, a, ntp, icon)
            if(icon .eq. 0) then
               sum = 0 
            else
               call cchg0c(-1, a, ntp, icon)
               if(icon .eq. 0) then
                  sum = 0
               endif
            endif
         elseif(sum .eq. -1) then
            call cchg0(-1, a, ntp, icon)
            if(icon .eq. 0) then
               sum = 0
            else
               call cchg0c(1, a, ntp, icon)
               if(icon .eq. 0) then
                  sum = 0
               endif
            endif
         endif
         if(sum .eq. 0) then
            icon = 0
         else
            ncc = ncc + 1
            icon = 1
         endif
      enddo
      end
!     
      subroutine cchgopposit(chg, a, ntp, icon)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
!         make one plus or minums charged ptcl into opposit charge
      integer chg  ! input.  1  or -1
      integer ntp
      type(ptcl):: a(ntp)
      integer icon  ! output.  0-> o%k
      real*8 u
      integer i, ncc
      logical found
      
      found = .false.
      ncc = 0
      do while ( .not. found .and.  ncc .lt. ntp*3)
         call rndc(u)
         i = ntp*u + 1
         found = a(i)%charge .eq. chg .and. 
     *       (a(i)%code .eq. kpion .or. a(i)%code .eq. kkaon)
         ncc = ncc + 1
      enddo
!
!      make this to opposit sign
!
      if(found) then
         icon = 0
         a(i)%charge =  -chg
      else
         icon = 1
      endif
      end

      subroutine cchg0(chg, a, ntp, icon)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"

!         make one plus or minums charge ptcl into 0
      integer chg  ! input.  1  or -1
      integer ntp
      type(ptcl):: a(ntp)
      integer icon

      real*8 u
      integer i, ncc
      logical found
      integer code, subcode

      found = .false.
      ncc = 0
      do while ( .not. found .and. ncc .lt. ntp*3)
         call rndc(u)
         i = ntp*u + 1
         found = a(i)%charge .eq. chg .and. 
     *        (a(i)%code .eq. kpion .or. a(i)%code .eq. kkaon)
         ncc = ncc + 1
      enddo
!
!      make this to 0
!
      if(found) then 
         icon = 0
         a(i)%charge =  0
         code = a(i)%code
         if(code .eq. kpion) then
            subcode = 0
         else
            call rndc(u)
            if(u .lt. 0.5) then
               call rndc(u)
               if(u .lt. 0.5) then
                  subcode = k0s
               else
                  subcode = -k0s
               endif
            else
               call rndc(u)
               if(u .lt. 0.5) then
                  subcode = k0l
               else
                  subcode = -k0l
               endif
            endif
         endif
         call cmkptc(code, subcode, 0, a(i))
         call cpm2e(a(i), a(i))
      else
         icon = 1
      endif
      end
!    
      subroutine cchg0c(chg, a, ntp, icon)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"

      integer chg
      integer ntp
      integer icon
      type(ptcl):: a(ntp)
      real*8 u
      integer i, ncc
      logical found
!       change neutral into charged
      integer code, subcode

      found = .false.
      ncc = 0

      do while ( .not. found .and. ncc .lt. ntp*3)
         call rndc(u)
         i = ntp*u + 1
         found = a(i)%charge .eq.  0 .and. 
     *        (a(i)%code .eq. kpion .or. a(i)%code .eq. kkaon)
         ncc = ncc + 1
      enddo
      if( found ) then
         a(i)%charge =  chg
         icon = 0
         code = a(i)%code
!           remake paticle; since mass is diff.
         call cmkptc(code, subcode, chg, a(i))
         call cpm2e(a(i), a(i))
      else
         icon = 1
      endif
      end
      

