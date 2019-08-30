#include "ZepMaxdef.h"
!!           Maxmedia;  moved to ZmediaLoft.h		    
!       integer Maxmedia, BitPhotoElec, BitPhoton,
        integer BitPhotoElec, BitPhoton,
     *    BitElectron, BitPositron, BitProton,
     *    BitNeutron, BitAntiNuc, BitDecay,
     *    BitOther
!!           Moved to ZmediaLoft.h			  
!        parameter(  Maxmedia = MAX_MEDIA,
       parameter(  
     *  BitPhotoElec=1, BitPhoton=2, BitElectron=3, BitPositron=3,
     *  BitProton=4, BitNeutron=5, BitAntiNuc=6, 
     *  BitDecay=7, BitOther=8  )
!  #include "Zmedia.h"
#include "ZepTrack.h"
#include "Zmove.h"
#include "Zep3Vec.h"

!          common variables used in tracking ptcls.
       type(epTrack)::  cTrack    !  current track
       type(epTrack)::  Incident  ! to keep 1st stacked track
       type(epmove)::  Move         ! to contain info for moving a
                                  ! track  
                   
	real*8  Ethin         ! Thin sampling threshold.
	real*8  Beta  !  v/c for MovedTrack; given if TimeStrucrue=T.
        real*8  Tcoefx, Tminx  ! used to truncate path
        

        integer Cn        !  current particle cn
        integer FirstCn   !  first collision (decay etc) occurred
                   ! in this component #
!!            moved to ZmediaLoft.h		    
!!       type(epmedia):: Media(Maxmedia)
!!        integer MediaNo   !  current particle media no = Cn2media(cn)
!!        integer NoOfMedia !  No of media 
!       >>>>>>>>>>>>>>>>>>>>>>light
        integer CrossMode ! 0-->crossing self boundary, 1-->crossing matreshka
                          ! 2-->crossing container which partially contains current comp.
        integer cLcompNo  !  current Light componet no
                          !   = Det.cmp(Cn).LightCompNo
                  	  !  Lcomp(cLcompNo)%comInfoNo is the index to
                          !  property info.  owned by the same property 
                          !  components
        integer cPtyNo    !  Current property No of the component
                          !  cPtyNo=Lcomp(cLcompNo)%comInfoNo. 
                  	  !  e.g, comInfo(cPtyNo)%refracIndex 
                          !  These are set whenever a ptlc move to a
                          !  new component
!       <<<<<<<<<<<<<<<<<<<<
        logical Firsti, FirstC ! to set  first  collision/interaction 
        logical ImperativeEmin ! if becomes T, residual range is not
                            ! considered when KE< Emin.
       type(epmedia):: FirstMedia  ! first collision media
        character*8  Proc1 ! to keep first interaction type
       type(epPos)::  FirstInt !  to keep first int. pos.
       type(epTrack):: FirstIntTrack ! to keep first int track ptcl
       type(ep3Vec)::  Bfield   !  Bfield.x, y, z  the B in  Tesla
                                !    in local coordinate 
        real*8   pathInB        !  total path length in B for a particle 
       type(ep3Vec)::  Efield   !  Efiled.x, y, z  the E in V/m
                                !   local coordinate.
        real*8 SumDe       !  sum  of all energy deposit

        real*8 MaxPath     !  MaxPathL of a comp. is converted into r.l
!!        real*8 Upsilon     !  Ee/me * B/Bcr for use with sync.
		               !  --> modepEMcontrol
!!        real*8 Xai         !  Eg/me * B/Bcr/2      !  --> modepEMcontrol
        real*8 Total1ryE   !  sum of K.E of 1ry
        integer Bndryerr   ! counter for unsuccessful boundary search in
	                   ! an  event.  Failure may occure if a particle
			   ! is just on a boundary of a component by chance.
  		           ! If this happens so many times in an event,
                           ! something wrong may be happing so that we
			     ! stop the executon ( in epbndry ).	       
        integer Making1ry  ! =1 during preparing 1ry for each event 
!
    !!       common /ZepTrackv/ Media, FirstMedia,
        common /ZepTrackv/ FirstMedia,
     *  Incident, cTrack, Move, 
     *  Bfield, Efield,
     *  FirstInt, FirstIntTrack, MaxPath, pathInB,
		     ! Upsilon, Xai,  
     *  Beta, Tcoefx, Tminx,  SumDe, Cn,  FirstCn,
!!     *  MediaNo,  NoOfMedia, Ethin, Total1ryE,
     *  Ethin, Total1ryE,
     *  Firsti, FirstC, Bndryerr, cLcompNo, cPtyNo,
     *  CrossMode, ImperativeEmin, Making1ry


       common /ZepTrackvC/ Proc1
