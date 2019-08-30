!     *************
      subroutine  epMCShard
!     *************
      use modMCScontrol
      implicit none
!
!
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
!
      real*8 E, prob,  tbrem, tknock, tanihi, t, dt, dl, dx,
     *   syncmfp
      real(8)::u
      logical:: posok, numok
      type(epPos):: posw
      integer:: i
       
      E = cTrack%p%fm%p(4)
      if( doNewMCS ) then
         if( (Media(MediaNo)%name /= 'sp' ) .and.
     *       (Media(MediaNo)%name /= 'hollow')  ) then
            call cfixMCSmodel( cTrack%p ) ! see energetically OK ?

            if( ActiveMCS /= 'Mol') then

               if( MCSzCond > 0 ) then
!                see if z is in the specifed region; else -> Mol
                  posok = .false.
                  call epl2w(cTrack%cn, cTrack%pos, posw)
                  do i = 1, MCSzCond
                     if(real(MCSzRange(i)) <= posw%z ) then
                        if( posw%z <= imag(MCSzRange(i)) ) then
                           posok = .true.
                           exit
                        endif
                     endif
                  enddo
               else
!                  no check about Z
                  posok=.true.
               endif

               if( MCSnumCond > 0 ) then
                  if( MCSdebug ) then
                     call epl2w(cTrack%cn, cTrack%pos, posw)
                  endif
! see if cn is in the ranage
                  numok = .false.
                  do i = 1, MCSnumCond
                     if(MCSnumRangeMin(i) <=  cTrack%cn ) then
                        if( cTrack%cn <= MCSnumRangeMax(i) ) then
                           numok = .true.
                           exit
                        endif
                     endif
                  enddo
               else
! no check about cn
                  numok = .true.
               endif

               if( MCSandor == 'and' ) then
                  if( posok .and. numok ) then
                     if( MCSrevert ) then
                        ActiveMCS = 'Mol'
                     else
!     keep current one
                     endif
                  else
                     if( .not. MCSrevert ) then
                        ActiveMCS = 'Mol'
                     else
!                          keep current
                     endif
                  endif
                  
               else             ! 'or'
                  if( posok .or.  numok ) then
                     if(  MCSrevert ) then
                        ActiveMCS = 'Mol'
                     else
                  !     keep current
                     endif
                  else
                     if( .not. MCSrevert ) then
                        ActiveMCS = 'Mol'
                     else
!                          keep current
                     endif
                  endif
               endif
               if( MCSdebug ) then 
                  if( ActiveMCS == 'Mol') then
                     write(0,'(a, 1p,3g14.4)') 'mol: ',  posw
                  else
                     write(0,'(a, 1p,3g14.4)') 'hin: ',  posw
                  endif
               endif
            else
               ActiveMCS = 'Mol'
            endif
         else
            ActiveMCS = 'Mol'
         endif
      endif
      
      if( ActiveMCS /= 'Mol' ) then
         call cfixMixedConst(MediaNo, int(cTrack%p%charge))
      endif


      if( ActiveMCS == 'El_hin') then
             ! sample hard scattering  mfp g/cm2.
         call cgetLamh( KEeV, mfpHardgr)  ! mfp of hard cs
!  in g/cm2
         call rndc(u)
         lHardgr = -log(u)* mfpHardgr  ! g/cm2
         lHardrl =lHardgr / Media(MediaNo)%X0g    ! rl
          !  r.l
          ! next is  cm length. this might  be used later for soft
          ! mcs treatment.
         lHardcm = lHardgr* Media(MediaNo)%gtocm
     *        / Media(MediaNo)%rhoc
         call csetIntInf(lHardrl, .false., 'hcs')
         
!         if( lHardrl < t ) then
!            Move%proc="hcs"     ! hard Coulomb scattring
!            t = lHardrl
!         endif
      endif

!      Move%dt = t
!      Move%dx = Move%dt * Media(MediaNo)%X0g
!      Move%dl = Move%dx * Media(MediaNo)%gtocm /
!     *    Media(MediaNo)%rhoc
      end
