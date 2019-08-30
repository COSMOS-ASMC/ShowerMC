      subroutine cutptcl(proj)
!///////////// when the calling place of crot3mom after particle generation
!            is moved to somewhere after the call to chookNEPI, 
!            we must change call cutptcl in  chookNEPI (proj can be simpliy
!            MovedTrack.p)
      use SoftenPiK
      implicit none
!          special purpose routine to see the importance of
!          energetic particles. 
! 
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zincidentv.h"
!
!       
      type(ptcl):: proj
      integer:: nout
      real(8):: E0
      integer:: i
      integer:: stackpos
      type(ptcl):: Tglab, Cmsp, Pjcms
      integer icon
      real(8):: y, eta
!//////////
      type(ptcl):: temp
!////////////////      
      
      modified = .false.

      if( IntInfArray(ProcessNo).process == "coll") then 
         E0lab = proj.fm.p(4)         
         if(useXinCMS) then
!               make proton (only for making CMS)
            call cmkptc(knuc, regptcl, 1, Tglab)
            Tglab.fm.p(1:3)=0.
            Tglab.fm.p(4) = Tglab.mass
!              get cms equivlent mass and 4 momentum
!            call cgeqm(MovedTrack.p, Tglab, Cmsp, icon)
            call cgeqm(proj, Tglab, Cmsp, icon)
            if(icon /=  0) then
               write(0,*) ' cms cannot be formed in cutptcl'
               stop
            endif
!             boost Pwork into CMS
            do i =1, Nproduced
               call cbst1(i, Cmsp, Pwork(i), Pwork(i))
            enddo
!                 projectile in cms
            call cbst1(2, Cmsp, proj, Pjcms)
!                      only forward X is manipulated; fwbw is not used 
            call csoftenFE(Pjcms,  1,  Pwork, Nproduced, nout)  !!!!Pwork(i) bad
            Nproduced = nout
!                 re-boost into Lab
            do i =1, Nproduced
               call cibst1(i, Cmsp, Pwork(i), Pwork(i))
            enddo
         else
!                in lab. do softening
!            call csoftenPiK(MovedTrack.p, Pwork,  Nproduced, nout)
            call csoftenPiK(proj, Pwork,  Nproduced, nout)
!///////////
!            write(0,*) ' E0=',MovedTrack.p.fm.p(4)
!            write(0,*) ' Npro=',Nproduced, ' nout=',nout
!////////////
            Nproduced = nout
         endif

         if(special) then
            do i = 1, Nproduced
               if( Pwork(i).code == kpion .or.
     *              Pwork(i).code == kkaon .or.
     *              Pwork(i).code == keta ) then
                  call cyeta(PWork(i), y, eta)
                  write(*,'("xd ",2i3,  1p,3g13.4)') 
     *            Pwork(i).code, Pwork(i).charge,
     *            Pwork(i).fm.p(4)/E0lab, y, eta
               endif
            enddo
         endif
      endif

      end
