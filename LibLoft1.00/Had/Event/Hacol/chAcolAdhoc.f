!         hadoron-A collision of my adhoc model.
      subroutine chAcolAdhoc(pj, ia, iz, a, ntp)
      implicit none
!          This is a Adhoc model of multiple production.
!          h-A collision is decomposed into multiple collision
!          of the leading particle, and for each collision
!          chncol is used with some reduction of the leading
!          particle energy for the 2nd, 3rd,... collisions.
!          Since the incident energy is high, we neglect
!          the multiple collisions at energy < 5 GeV.
!

#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnp.h"
!
      type(ptcl):: pj   ! input  projectile ptcl
      integer ia    ! input. mass no. of target
      integer iz    ! input. charge no. of target
      type(ptcl):: a(*)  !  output. produced ptcls
      integer ntp   ! number of produced ptcls
!
      integer ncoll, tgtchg, i, n
      type(ptcl):: aPtcl,  tgt
!      real*8  eminSucCol/5./
      real*8  eminSucCol/3./
      integer icon, fails
!
!          Fermi momentum can be neglected ( E > Elund= 500 GeV)
      if(SucInt .eq. 0) then
         call csampCollInA(pj, ia, ncoll)
      elseif(SucInt .eq. 1) then
         call csmpColInA2(pj, ia, ncoll)
      else
         call cerrorMsg('SucInt has an invalid value', 0)
      endif
!            fix target nucleon charge
         call cfxTgtChg(ia, iz, tgtchg)
!                make target nucleon
         call cmkptc(knuc, regptcl, tgtchg, tgt)
!              give 4 momentum
         tgt%fm%p(1) = 0.
         tgt%fm%p(2) = 0.
         tgt%fm%p(3) = 0.
         tgt%fm%p(4) = tgt%mass
!         
         fails = 0
         icon = 1
         do while( fails .lt. 10 .and. icon .ne. 0 )
            call chncol(pj, tgt, a, ntp, icon)
            if(icon .ne. 0) then
               fails = fails + 1
            endif
         enddo

         if(icon .ne. 0) then
!            generation failed.
            a(1) = pj
            ntp = 1
            goto 100
         endif

         call cslpx2(.true.)      ! specify that this is 2nd,3rd .. col.
         do i = 2, ncoll
!                extract leading ptcl
            aPtcl = a(ntp)
            if(aPtcl%fm%p(4) .gt. eminSucCol) then
                  ntp = ntp - 1
                  call cfxTgtChg(ia, iz, tgtchg)
                  tgt%charge = tgtchg
                  call chncol(aPtcl, tgt, a(ntp+1), n, icon)
                  if(icon .ne. 0) then
                     ntp = ntp + 1
                     a(ntp) = aPtcl
                     goto 50
                  endif
                  ntp = ntp + n
             endif
         enddo 
 50      continue 
         call cslpx2(.false.)  ! reset nucleus condition 
 100     continue
      end
