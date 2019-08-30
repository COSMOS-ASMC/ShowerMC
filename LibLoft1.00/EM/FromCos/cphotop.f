      subroutine cphotop( pj )
!     ***********
!     use modXsecMedia
      use modColInfo    ! for TargetNucloeonO ...
      use modEMcontrol
      implicit none

#include  "ZmediaLoft.h"      
#include  "Zcode.h"
#include  "Zptcl.h"
#include  "Zpwork.h"
! #include  "Zptcl.h"  ! modColInfo contans this
      

!
      type(ptcl),intent(in) :: pj  ! incident photon
      
!
      integer ngen
      character(8):: whichcode

      whichcode = " "
      Nproduced = 0
      if( HowPhotoP == 0 ) then  ! will not happen
         ngen = 1
!     Pwork(Nproduced+1) =  MovedTrack%p
         Pwork(Nproduced+1) =  pj
      elseif( HowPhotoP == 1 ) then
         whichcode = "sofia"
      elseif( HowPhotoP == 2 ) then
!     if( MovedTrack%p%fm%p(4) < 2.5 )  then
         if( pj%fm%p(4) < 2.5 )  then
            whichcode ="current"
         else
            whichcode ="sofia"
         endif
      elseif( HowPhotoP == 3 ) then
!     if( MovedTrack%p%fm%p(4) < 2.5 )  then
         if( pj%fm%p(4) < 2.5 )  then
            whichcode ="sofia"
         else
            whichcode ="current"
         endif
      elseif( HowPhotoP == 4 ) then
         whichcode = "current"
      else
         write(0,*) 'HowPhotoP =', HowPhotoP, ' invalid '
         stop
      endif

!      call cfixTarget( Media(MediaNo) )   not needed. though not harmufl
      
      if( whichcode == "sofia" ) then
         call csofia( TargetNucleonNo, TargetProtonNo, 
!     *      MovedTrack%p, Pwork(Nproduced+1), ngen)
     *      pj, Pwork(Nproduced+1), ngen)
      elseif( whichcode == "current" ) then
         call cgpHad(TargetNucleonNo, TargetProtonNo, 
!     *        MovedTrack%p, Pwork(Nproduced+1), ngen)
     *        pj,  Pwork(Nproduced+1), ngen)
      elseif( whichcode == " " ) then
         write(0,*) ' setting mistake of  whichcode=',
     *    whichcode, ' in cphotop'
         stop
      endif
      Nproduced = Nproduced + ngen
      end
