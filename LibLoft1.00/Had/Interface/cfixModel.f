      subroutine cfixModel( aPtcl )
      implicit none
!
!        fix hadronic interaction model (including 
!        photoproduction of hadrons, and possibly
!        muon nucear interaction).
!        and set ActiveMdl
!
#include "Zptcl.h"
#include "Zevhnp.h"
#include "Zevhnv.h"
#include "Zcode.h"

      type(ptcl)::aPtcl ! input. particle for hadronic interaction

      real*8 kepn 
      integer i
      integer :: icon

      kepn = aPtcl%fm%p(4)- aPtcl%mass
      if( aPtcl%code .eq. kgnuc ) then
         kepn =  kepn / aPtcl%subcode
      endif

      do i = 1, NoOfMdls
         if( kepn  .lt. InteErg(i) )  then
            ActiveMdl = ModelList(i)
            goto 10
         endif
      enddo

      call cerrorMsg('IntModel list shown below is strange', 1)
      call cerrorMsg(IntModel, 1)
      call cpdpmjetinp
      write(0,*) ' current K%E/n=',kepn,
     *  ' code=',aPtcl%code, ' subcode=',
     *  aPtcl%subcode, ' charge=',aPtcl%charge
      write(0,*) '  NoOfMdls=', NoOfMdls
      do i = 1,  NoOfMdls
         write(0,*) ' i=', i, ' InteErg(i)=',
     *   InteErg(i)
      enddo
      stop 12345
 10   continue

      do i = 1, NoOfMdls2
         if( kepn  .lt. InteErg2(i)) then
            ActiveMdl2 = ModelList2(i) ! used only for mfp calc
            goto 15
         endif
      enddo

      call cerrorMsg('XsecModel list shown below is strange', 1)
      call cerrorMsg(XsecModel, 1)
      call cpdpmjetinp
      write(0,*) ' current K%E/n=',kepn,
     *  ' code=',aPtcl%code, ' subcode=',
     *  aPtcl%subcode, ' charge=',aPtcl%charge
      write(0,*) '  NoOfMdls2=', NoOfMdls2
      do i = 1,  NoOfMdls2
         write(0,*) ' i=', i, ' InteErg2(i)=',
     *   InteErg2(i)
      enddo
!!      stop 123456
      stop 23456    !! jaxa length should be <=5
 15   continue

      if( aPtcl%code == knuc .and.  aPtcl%subcode == antip ) then
         if(ActiveMdl == "phits" ) then
          !  ' for p-bar or n-bar, phits cannot be used'
            call epRescue(icon)
            if( icon /= 0 ) then
               call epRescueErr('phits cannot be used for p-bar', aPtcl)
               stop
            endif
         endif
      elseif(aPtcl%code == kkaon ) then
         if(ActiveMdl == "phits" ) then
            call epRescue(icon)
            if( icon /= 0 ) then
               call epRescueErr('phits cannot be used for Kaon', aPtcl)
               stop
            endif
         endif
      endif
      ! If photon or muon, they will be 
         !  managed in photon interaction routine
      ! other hadronic particles,say, lambda0, will be decayed
      ! if the model cannot cope with it.  (cseeColPossible).
 20   continue      
      end
      subroutine epRescue(icon)
      implicit none
!
!        fix hadronic interaction model (including 
!        photoproduction of hadrons, and possibly
!        muon nucear interaction).
!        and set ActiveMdl
!
#include "Zptcl.h"
#include "Zevhnp.h"
#include "Zevhnv.h"
#include "Zcode.h"
      integer,intent(out):: icon  ! 0.  rescuee could be set
                           !    1.    no resucee could be fount
!             see if jam or dpmjet3 can be used, if not error
      integer::i

      icon  = 0
      do i = 1, NoOfMdls
         if( ModelList(i) == "dpmjet3") then
            ActiveMdl= "dpmjet3"
            goto 20
         elseif(ModelList(i) == "jam") then
            ActiveMdl= "jam"
            goto 20
         elseif( ModelList(i) == "nucrin") then
            ActiveMdl= "nucrin"
            goto 20
         endif
      enddo
      icon = 1
 20   continue
      end

      subroutine epRescueErr(msg, aPtcl)
      implicit none
#include "Zptcl.h"
      character(*),intent(in):: msg
      type(ptcl)::aPtcl !  input. problematic  ptcl

      write(0,'(a)')  trim(msg),
     * ' so that dpmjet3 or jam  (or nucrin)  must be given in ',
     * ' the IntModel list as a rescuee: ',
     * " If you don't want to use such a rescuee as a normal ",
     * ' interaction model, you may specify it as follows:',
     * " IntModel='"//'"abc" E1 "dpmjet3" E2 "xyz"',
     * ' and set E1=E2, then, model abc is used below E1=E2. ',
     * ' When energy is >= E1, model xyz is used. So dpmjet3 will',
     * ' not be used at any energy. It will be used only as ',
     * ' rescuee.',
     * ' ',
     * 'particle info: '
      write(0,'(a, g13.4)') ' ptcl code=',aPtcl%code,
     *               '  charge  =',aPtcl%charge,
     *               '  subcode =',aPtcl%subcode,
     *              '  Ek      =', aPtcl%fm%p(4)-aPtcl%mass, ' GeV'
      end
!       
      subroutine cqActiveMdl(model)
      implicit none
!
!        inquire the current active model
!
#include "Zptcl.h"
#include "Zevhnv.h"

      character*16 model
      model = ActiveMdl
      end
      subroutine cqActiveMdl2(model)
      implicit none
#include "Zptcl.h"
#include "Zevhnv.h"

      character*16 model
      model = ActiveMdl2
      end
