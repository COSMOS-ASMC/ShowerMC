      subroutine cseeColPossible(pj,  icon)
      implicit none

#include "Zptcl.h"
!  #include "Zkfcode.h"
! #include "Zmanagerp.h"
! #include "Zevhnp.h"
#include "Zevhnv.h"
#include "Zcode.h"
! #include "Zmass.h"
      type(ptcl):: pj  ! input.  hadronic projectile particle
                        !  which is  going to make the collision
                ! fixed by epfixProc or cfixProc 
       ! **** Must not be called for other interactions ***** 
      integer,intent(out)::icon  ! 0. The interaction is
                    ! acceptable. 
                    !-1 not acceptable. This happend when
                    ! the interacion is collision and 
                    ! the current interacion model cannot
                    ! manage it.  The user must force the
                    ! projectile to decay. (pj is Simga etc)

      icon = 0
      select case(ActiveMdl)
      case('phits')
         if(  pj%code == kpion .and. pj%charge == 0) then
            icon = -1
!         elseif( pj.code == kkaon .and. pj.charge == 0 ) then
         elseif( pj%code == kkaon ) then ! from v 7.633  
                           ! at cfixModel, phits might be
                           ! already not assigned
            icon  = -1
         elseif( pj%code >=  klambda .or. pj%code == kdmes ) then
            icon = -1
         endif
      case('jam') 
         if(  pj%code == kpion .and. pj%charge == 0) then
            icon = -1   ! revived v9.153

!         elseif( pj.code == kkaon .and. pj.charge == 0 ) then
!            original jam should be ok
!            icon  = -1
!         elseif( pj.code >=  klambda .or. pj.code == kdmes ) then
!            cjam.f manage these as n or p
!            icon = -1
!
         endif
         
!
!      data RegMdls/'phits', 'jam', 'dpmjet3', 'fritiof7.02',
!     *      'fritiof1.6',
!     *     'gheisha', 'nucrin', 'ad-hoc', 'incdpm3', 'qgsjet2'/
!
      case('nucrin' ) 
!              seems ok ??
!         if(  pj.code == kpion .and. pj.charge == 0) then
!            icon = -1
!         elseif( pj.code == kkaon .and. pj.charge == 0 ) then
!            icon  = -1
         if( pj%code >=  klambda .or. pj%code == kdmes ) then
            icon = -1
         endif

      case('epos')    
!         if(  pj.code == kdmes .or. pj.code == 'kgzai') then
         if(  pj%code == kdmes ) then
            icon = -1
         endif
      case('sibyll') 
         if( pj%code >  klambda .or. pj%code == kdmes ) then
            icon = -1
         endif
      case('fritiof1.6')         
!           seems ok??
!         if(  pj.code == kpion .and. pj.charge == 0) then
!            icon = -1
!         elseif( pj.code == kkaon .and. pj.charge == 0 ) then
!            icon  = -1
         if( pj%code >=  klambda .or. pj%code == kdmes ) then
            icon = -1
         endif
      case('dpmjet3')
         if( pj%code == keta ) then
            icon = -1
         elseif( pj%code == kgzai .and. pj%fm%p(4)<= 9.d0)  then
            icon = -1
         endif
      end select
!           managed inside, mostly
!      case('dpmjet3')
!      case('qgsjet2')
!      case('gheisha')
!      case('fritiof7.02')
!      case('incdpm3')
!      case('ad-hoc')

      end
