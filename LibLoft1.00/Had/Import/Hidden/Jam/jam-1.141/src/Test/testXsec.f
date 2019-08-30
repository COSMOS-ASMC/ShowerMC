#include "../jamdat.f"
#include "mydummy.f"
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "jamextrn.inc"
      integer:: msel, icltyp, kf1, kf2, mchanel, mabsrb, ijet, icon
      real(8)::srt, pr, em1, em2, sig, sigel, sigin(30)
      record /ptcl/tg, pj
      real(8):: Eklab, s, prsq, plab
      integer code, subc, charge, ia, iz
c///////
      call mydummy   ! to confirm block data
c///////////
      msel = 1

      write(0,*) 'Enter: projectile code,subc, charg '
      read(*,*)  code, subc,  charge
      if(code == 5 .and. charge == 0 ) then
         subc = k0l
      endif
      call cmkptc(code, subc, charge, pj)
      subc = pj.subcode
      call ccos2kf(code, subc, charge, kf1)
      em1 = pj.mass


      if(code == 6 ) then
         if( subc  == antip ) then
            icltyp = 4     ! pbar nbar
         else
            icltyp = 1          ! pp pn,
         endif
      else
         icltyp = 2    ! pip, pin, kp, kn
      endif
      write(0,*) 'proj=', code, subc, charge
      write(0,*) 'Enter: target code, charge (p or n)'
      subc = -1
      read(*,*)  code,  charge
      call cmkptc(code, subc, charge, tg)
      subc = tg.subcode
      call ccos2kf(code, subc, charge, kf2)
      em2 =  tg.mass
      ia = 1
      iz = tg.charge
      pj.fm.p(1:2) = 0.
      pj.fm.p(4) = 2.
      pj.fm.p(3) = sqrt(pj.fm.p(4)**2 - pj.mass**2)
      call cjamini(pj, tg)
      Eklab=0.2
      do while (Eklab < 10000.) 
         plab = sqrt((Eklab+pj.mass)**2 - pj.mass**2)
         s = em1**2 + em2**2 + 2*em2*(Eklab + em1)
         prsq= (s-(em1+em2)**2)*(s-(em1-em2)**2)/(4*s) 
         pr = sqrt(prsq)
         srt = sqrt(s)
         call jamcross(msel,icltyp,srt,pr,kf1,kf2,em1,em2,
     $       sig,sigel,sigin,mchanel,mabsrb,ijet,icon)
         write(*,'(1p,4g12.3,i3)') Eklab, sig, sigel, plab, icon
         Eklab = Eklab * 10.**0.1
      enddo
      end
