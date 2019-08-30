!             hadron A collision by New Lund (Fritiof v7.02)
        subroutine chANewLund(pj, ia, iz,  a, ntp)
        implicit none

#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnv.h"
!
      type(ptcl):: pj   ! input  projectile ptcl
      integer ia    ! input. nucleon no. of target
      integer iz    ! input. charge no. of target
      type(ptcl):: a(*)  !  output. produced ptcls
      integer ntp   ! number of produced ptcls
!
      type(ptcl)::   tgt, temp
      integer icon, kf, i, code, subcode, charge, maxi
      real*8 maxe
!      type(fmom):: gb  ! (g*beta, gc)
      character*100 msg
!////////////
!       logical deb
!       common /cccdeb/deb
!      integer seed(2)
!      debug = .false.
!////////////
!

!             make target nucleon simply to form cms
         call cmkptc(knuc, regptcl, 1,  tgt)
!              give 4 momentum
         tgt.fm.p(1) = 0.
         tgt.fm.p(2) = 0.
         tgt.fm.p(3) = 0.
         tgt.fm.p(4) = tgt.mass
!
!         make projectile so that its direction is z
!
         Pjlab.fm.p(1) = 0.
         Pjlab.fm.p(2) = 0.
         Pjlab.fm.p(3) = sqrt(pj.fm.p(4)**2-pj.mass**2)
         Pjlab.fm.p(4) = pj.fm.p(4)
         Pjlab.mass = pj.mass
!         Pjlab.code = pj.code
!         Pjlab.subcode = pj.subcode
!         Pjlab.charge = pj.charge

!          get cms equivlent mass and 4 momentum
       call cgeqm(Pjlab, tgt, Cmsp, icon)
       if(icon .ne. 0) then
          write(msg, *) 
     *    ' cms cannot be formed in chANewLund; proj and ',
     *    'target are '
          call cerrorMsg(msg, 1)
          call cprptc(Pjlab, 1)
          call cprptc(tgt, 1)
          stop 9999
       endif
!          get (g*beta, gc) of cms
!       call cgetlf(Cmsp,  gb)
!          boost pj into cms.
!       call cbst0(1, gb, Pjlab, Pjcms)
!
       code = pj.code         ! need substitution due to integer*2
       subcode = pj.subcode
       charge = pj.charge
!          convert particle code to kf code
       call ccos2kf(code, subcode, charge, kf)
!          generate ptcls in cms and set them in a
!/////////////
!        debug = pj.fm.p(4) .gt. 30.e3
!        if(debug) then
!           write(*,*) '-----',
!     *      pj.fm.p(1), pj.fm.p(2), pj.fm.p(3), pj.fm.p(4)
!           write(*,*) 'ia, iz=', ia,iz,code,subcode,charge, kf
!           call rnd1s(seed)
!           write(*,*) ' seed=', seed
!        endif
!       if( deb ) then
!          write(*,*) ' bef cfrevent '
!       endif
!///////////

       call cfrevent(kf, charge, ia, iz, Cmsp.mass, a, ntp)

!////////////////
!       if( deb ) then
!          write(*,*) ' after cfrevent '
!       endif
!////////

!          boost back to lab.
!          find max energy partcle and store it in the last
!          part. (to save stack area, later)
       maxe = -1. 
       do   i=1, ntp
          call cibst1(i, Cmsp, a(i), a(i))
          if(a(i).fm.p(4) .gt. maxe) then
             maxe = a(i).fm.p(4)
             maxi = i
          endif
! &&&&&&&&&&&&&&&&777
!          write(*,*) ' code=', a(i).code, ' sub=', a(i).subcode, 
!     *     ' chg=', a(i).charge, ' KE=', a(i).fm.p(4) - a(i).mass
!  &&&&&&&&&&&&&&&&
       enddo
       temp = a(maxi)
       a(maxi) = a(ntp)
       a(ntp) = temp
!         in  above, momentum  is measured from pj.fm so 
!         we need rotate it.
       call crot3mom(pj, a, ntp)
       end

