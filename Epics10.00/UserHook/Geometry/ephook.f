#include "ZepMain.f"
#include "Zepcondc.h"
#if defined (INTINFO)
#include "epUI.f"
#endif
!        **************************************************************
!        *
!        * ephook:  collection of subroutines which should be managed by
!        *      the user.  This is a sample program to count the  energy
!        *      loss at given layers.
!        *  
!        **************************************************************
!         
!
!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   +  Note that besides the subroutines listed below in this module,
!   +  the following inquirey subroutines are available for use:
!   +
!   +  call epqFirstI(pos);
!         #include "ZepPos.h"
!         record /epPos/ pos ; output.( pos.x, pos.y, pos.z )
!                              First interaction position in world
!
!      call epqfirstM(firstM)
!        #include "Zmedia.h"
!        record /epmedia/ firstM
!        firstM.colElem; integer. first interaction element number
!              For e/g incident, this will be 0 and you cannot use
!              other firstM.
!        firstM.colA;    integer. first interaction target mass number 
!        firstM.colZ;    integer. first //                atomic 
!        firstM.name;    character*8.  media name
!        firstM.elem(firstM.colElem).A is a real number and is not necessary
!              equal to firstM.colA (which may be, e.g,  36.4)
!        call epqFirstP(proc)
!          character*8 proc;  get first interaction type. (say, "col" "deca"
!          "pair" etc)
!       
!   +  call epqmat(i, mat): i; integer input.  component number
!   +                       mat: characeter*8 output.  mediea
!   +                            such as 'Si' of  the i-th component.
!   +  call epqstn(n):  n; integer output.  get current stack depth
!
!   +  call epqsTrack(n, sTrack): inquire stacked track.
!                    n; integer input. stack depth
!                   record /epTrack/ sTrack. Energy=0 means n is invalid
!   +  call epqevn(nev):  nev; integer output. Event number already created
!                              in this run.
!   +  call epqinc(aTrack)
!                 inquire incident particle.
!   +         record /epTrack/ aTrack
!   +             if multiple particles are incident, the first one is obtained

!   +  call epqncp(ncomp) inquire the total number of components
!                   integer ncomp  output.
!
!   +  call sqtevn(nev): nev; integer output.  total number of events
!   +                                           created so far.
!   +
!   +  call sqcont(icon):  icon; integer output. 
!   +                            icon == 0 ==>  first job
!   +                            icon != 0 ==>  continued job
!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        
!
!       *****************************
!       *        At this momemet, all the system level initialization 
!       *   has been ended.  Do your own initialization for all the
!       *   events.
!
        subroutine uiaev
        use epModGeometry
        implicit none
        call epGeom0
        end

!       **************************************************
!       *  At this moment, all the system level initialization for
!       *  an event has been ended.  Do your own initializaiton for
!       *  the event.   Note, at this moment, no incident particle
!       *  has been made yet.
!       *
        subroutine ui1ev
        implicit none
        end
!       **************************************************
!       *  If the default incident particle treatment does not
!       *  suffice you, make your own incident particle here.
!       * (you can  make even multiple particles as  the incident).
!       *  If default is ok, don't touch this.
        subroutine usetip(icon)
        implicit none
        integer icon
!            
!        If you make the incident yourself.
!             make icon = 0 on return, 
!             ^^^^^^^^^^^^^^^^^^^^^^^^
!        To make the incident, you may make the next call
!           as many times as you want.
! 
!        call epputTrack(aTrack)
! ******     if the number of incident partilces could be  ****
!            very large (> 3900), it's  safe to use
!        call epputTrk2(io, aTrack) where io is the
!            binary working file logical number; you must
!            open the file your self.
!           record /epTrack/ aTrack
!              
            icon=1
        end
        subroutine uafi1ev
        use epModGeometry
        implicit none
!          this is called when all the initialization for
!         an event has been ended.  (incident particle has
!         set).  Do your own final init. for the event
        call epGeom1
        end

!       ***********************************************
!       *   This is called 
!       *     when a charged particle loses
!       *     energy in a component for which you have specifed to do
!       *     energy loss counting.
!       * Or
!       *     when a photon energy becomes lower than a given minimum
!       *     energy.  It's your decision whether you discard such 
!       *     a photon or count as energy lost here.
!  
        subroutine userde(info, aTrack, Move, media)
        implicit none
#include "ZepTrack.h"
#include "Zmove.h"
#include "Zmedia.h"

!          next is needed if you want to judge particle code.
!  #include "Zcode.h"     


        integer info             ! input. 0--> particle is still active.
                                 !        1--> particle is dying.
                                 !        (i.e, after moving Move.dl, it 
                                 !         dies) .
       type(epTrack)::  aTrack  ! input. current track. before it is moved.
       type(epmove)::   Move      ! input. containes info of moved track.
!                                ! output. You have to set Move.Abort
!         Move.Abort = 0 if you want to continue the simulation of
!                           this event normally.
!                    = 3 if you want to discard this particle
!                        For example, you may use this if you use
!                        formula for further development due to this
!                        particle 
!                    = 1 if you want to abort the generation of 
!                           this event, but want to execute ue1ev.
!                    = 2 if you want to abort the genration of 
!                         this event, and skip ue1ev.
       type(epmedia):: media
      Move%Abort = 0
!          You can use the following info.
!   ================================================================
!     aTrack.p         : the same particle record as Cosmos, i.e.,
!     aTrack.p.fm.p(1), p(2), p(3): momentum of the particle in GeV
!     aTrack.p.fm.p(4)            : total energy of    //
!     aTrack.p.mass               : mass  //                 
!     aTrack.p.code               : particle code 
!     aTrack.p.subcode            : particle subcode 
!     aTrack.p.charge             : charge  //.
!     aTrack.pos.x, y, z          : local coordinate of the particle in cm
!                                   If you need world coordinate, use
!                          call epl2w(aTrack.cn, aTrack.pos, wpos) 
!                    where
!                   record /epPos/ wpos is the output  world coordinate
!     aTrack.w.x, y, z            : direction cosines in the local coordinate
!                                   If you need world coordinate values, use
!                          call epl2wd(aTrack.cn, aTrack.w, ww)     
!                    where
!                   record /epDirec/ ww is the output  values in world coord.
!
!     aTrack.t                    : time in sum of (cm/beta) 
!     aTrack.wgt                  : weight of the particle when thinning is
!                                   done (normall 1.0)
!     aTrack.cn                   : component number
!
!     Move.Track      :   Info of moved track.
!     Move.boundary   :   Component boundary position when the particle crosses
!                         the boundary. (Move.Cross =T)
!                         Move.boundary.x, y, z. in cm. 
!     Move.proc       :   If the particle is going to interact, the process
!                            name is set. proc= one of
!                        'brem': bremstralhng
!                        'knoc': knock on (bhabha or moller scattering)
!                                other heavy particle knock-on
!                        'anih': positron anihilation
!                        'pair' : pair creation
!                        'comp' : compton scattering
!                        'phot' : photo electric effect
!                        'photop': photo-hadron production
!                        'coll'  : hadron's nuclear interaction
!                        'decay' : decay 
!                        '    '  : no interaction yet.
!    Move.dl          :  Path length in cm
!    Move.dE          :  Energy loss during dl. in GeV
!    Move.dEeff       :  Effective energy loss (use this for counting) GeV
!    Move.dEion       :  Energy loss due to ionization loss (GeV)
!    Move.dx          :  path length in g/cm^2
!    Move.dt          :  path length in r.l
!    Move.Cross       :  Becomes T if the ptcl crosses the boundary.
!    Move.Trunc       :  Becomes T if the ptcl track is truncated.
!    Move.Abort       :  must be set by the user.
!    
!    If you want to use a particle code, like,
!        if(aTrack.p.code .eq.  kphoton) ...
!    You have to uncomment Zcode.h above; they are the same as
!    Cosmos. Some of popular ones are:
!         kphoton:  photon             kelec  :  electron
!         kmuon  :  muon               kpion  :  pion
!         kkaon  :  kaon               knuc   :  nucleon
!         kneue  :  electron neutrino  kneumu :  muon neutrino
!         kalfa  :  He                 klibe  :  LiBeB group
!         kcno   :  CNO   group        khvy   :  Na/Mg/Si group
!         kvhvy  :  S/Cl/Ar group      kiron  :  Fe group
!         kdmes  :  D  meson           keta   :  eta  meson
!    Subcode may sometimes be needed:
!         regptcl:  particle           antip  : anti-particle
!         k0s    :  k0short            k0l    : k0 long
! 
        end
!       ************************************************
!       *  This is called when a particle passes the boundary
!       *  of a component for which you have specifed to
!       *  count particle nubmer.
!
        subroutine userbd(info, aTrack, Move, media) 
        use epModGeometry
        implicit none
#include "ZepTrack.h"
#include "Zmove.h"
#include "Zmedia.h"

!          next is needed if you want to judge particle code.
!  #include "Zcode.h"     

        integer info              ! 0--> ptcl is exiting to void
                                  ! <0 --> ptcl is exiting to |info| comp.
                                  ! >0 --> ptcl is entering from info comp.
       type(epTrack):: aTrack    ! input. Current track before it is moved.
                                  !  If info <=0, the track position is 
                                  !  somewhere inside the component from which
                                  !  the track is existing.
                                  !  If info >0, aTrack.pos is the position 
                                  !  just before exiting the prvious component
       type(epmove)::  Move        ! input/output. 
                                  ! Move.Track is the track
                                  ! infomation of the current particle
                                  ! moved to a new position.
                                  ! Say, Move.Track.cn is the  current
                                  ! comp. number.  For other details,
                                  ! see userde.

!
!                      comp.1    comp.2
!          info<0    |  *-----x|          |        
!                                                   * is aTrack.pos
!                                                   x is Move.Track.pos
!          info>0    |        *|x         |       
!  
       type(epmedia):: media       ! input.
        call epGeomB(info, aTrack, Move, media) 
        Move%Abort = 0
        end
!       **************************************************
!       *  This is called when all the system level "end process"
!       *  for the event has been ended.  Do your own end process
!       *  for the event.
!       *
        subroutine ue1ev
        use epModGeometry
        implicit none
#include "Zmedia.h"
#include "ZepPos.h"
        call epGeomE1ev
        end
!       *************************************************
!       *  This is called when all the system level "end process"
!       * for all the events has been  ended. Do your own end
!       * process for the events

      subroutine ueaev
      implicit none
      end
