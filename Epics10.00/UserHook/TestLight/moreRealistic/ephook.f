#include "../../main.f"
#include "Zepcondc.h"
#if defined (INTINFO)
#include "epUI.f"
#endif
#include "modtemp.f"
c        **************************************************************
c        *
c        * ephook:  collection of subroutines which should be managed by
c        *      the user.  This is a sample program to count the  energy
c        *      loss at given layers.
c        *  
c        **************************************************************
c         
c
c   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   +  Note that besides the subroutines listed below in this module,
c   +  the following inquirey subroutines are available for use:
c   +
c   +  call epqFirstI(pos);
c         #include "ZepPos.h"
c         record /epPos/ pos ; output.( pos.x, pos.y, pos.z )
c                              First interaction position in world
c
c      call epqfirstM(firstM)
c        #include "Zmedia.h"
c        record /epmedia/ firstM
c        firstM.colElem; integer. first interaction element number
c        firstM.colA;    integer. first interaction target mass number 
c        firstM.colZ;    integer. first //                atomic 
c        firstM.name;    character*8.  media name
c        firstM.elem(firstM.colElem).A is real number and is not necessary
c                         equal to firstM.colA (which may be, e.g,  36.4)
c       
c   +  call epqmat(i, mat): i; integer input.  component number
c   +                       mat: characeter*8 output.  mediea
c   +                            such as 'Si' of  the i-th component.
c   +  call epqstn(n):  n; integer output.  get current stack depth
c
c   +  call epqsTrack(n, sTrack): inquire stacked track.
c                    n; integer input. stack depth
c                   record /epTrack/ sTrack. Energy=0 means n is invalid
c   +  call epqevn(nev):  nev; integer output. Event number already created
c                              in this run.
c   +  call epqinc(aTrack)
c                 inquire incident particle.
c   +         record /epTrack/ aTrack
c   +             if multiple particles are incident, the first one is obtained

c   +  call epqncp(ncomp) inquire the total number of components
c                   integer ncomp  output.
c
c   +  call sqtevn(nev): nev; integer output.  total number of events
c   +                                           created so far.
c   +
c   +  call sqcont(icon):  icon; integer output. 
c   +                            icon == 0 ==>  first job
c   +                            icon != 0 ==>  continued job
c   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c        
c
c       *****************************
c       *        At this momemet, all the system level initialization 
c       *   has been ended.  Do your own initialization for all the
c       *   events.
c
        subroutine uiaev
        use modtemp
        implicit none
#include "ZepTrackp.h"


        integer i,j

        integer::deact(15)=(/5963, 5967, 5975, 5979, 5983, 5987,
     *     6023, 6027, 6035, 6039, 6043,  6079,
     *     6095, 6099, 6103/ )


        call epqncp(Ncomp)   ! inquire the number of components
        



!         deactivate light transport for some component
!        j = sizeof(deact)/4   !  sizeof gives byte size
        j = size(deact)       ! this is array size
        write(0,*) ' array size of deact is ', j
        do  i = 1, j
           call epLightResetCountDE(deact(i), 2)      ! scintillator
           call epLightResetCountDE(deact(i)+1, 0)    ! sensor
           call epLightResetCountDE(deact(i)+2, 0)    ! sensor
        enddo
!         activate light  transport for some component; 
        call  epLightResetCountDE(6008, 3112)  ! mn = 11 not 10 for x-direction
!                        because surface # assignment is different from 10
!         but there is no
!         sensor in this case;  sensor ouput is not possible.
!         we make a photon count at all surface # for comp. 6008.
        call  epResetCountIO(6008,  2)  !  2 here means when photon is exiting
                            !  if sensor is spcified, 1 may be given to count
                            !  entering photons.

       ! ******after resetting countDE/countIO, we must call next*******
       !   Next is no more needed. If it is called withouth argument,
       !   (i.e,   call epLi...DE ), it's ok but redundant
!        call epLightCountDE( Light ) ! arg  is better to be the same as Light in
            !  sepicsfile.  (at least one of 11,12,21,22.)

        end

c       **************************************************
c       *  At this moment, all the system level initialization for
c       *  an event has been ended.  Do your own initializaiton for
c       *  the event.   Note, at this moment, no incident particle
c       *  has been made yet.
c       *
        subroutine ui1ev
        use modtemp
        implicit none



        end
c       **************************************************
c       *  If the default incident particle treatment does not
c       *  suffice you, make your own incident particle here.
c       * (you can  make even multiple particles as  the incident).
c       *  If default is ok, don't touch this.
        subroutine usetip(icon)
        implicit none
        integer icon
c            
c        If you make the incident yourself.
c             make icon = 0 on return, 
c             ^^^^^^^^^^^^^^^^^^^^^^^^
c        To make the incident, you may make the next call
c           as many times as you want.
c 
c        call epputTrack(aTrack)
c ******     if the number of incident partilces could be  ****
c            very large (> 3900), it's  safe to use
c        call epputTrk2(io, aTrack) where io is the
c            binary working file logical number; you must
c            open the file your self.
c           record /epTrack/ aTrack
c              
            icon=1
        end
        subroutine uafi1ev
c          this is called when all the initialization for
c         an event has been ended.  (incident particle has
c         set).  Do your own final init. for the event
        end

c       ***********************************************
c       *   This is called 
c       *     when a charged particle loses
c       *     energy in a component for which you have specifed to do
c       *     energy loss counting.
c       * Or
c       *     when a photon energy becomes lower than a given minimum
c       *     energy.  It's your decision whether you discard such 
c       *     a photon or count as energy lost here.
c  
        subroutine userde(info, aTrack, Move, media)
        use modtemp
        implicit none
#include "ZepTrack.h"
#include "Zmove.h"
#include "Zmedia.h"

c          next is needed if you want to judge particle code.
c  #include "Zcode.h"     


        integer info             ! input. 0--> particle is still active.
                                 !        1--> particle is dying.
                                 !        (i.e, after moving Move.dl, it 
                                 !         dies) .
        record /epTrack/ aTrack  ! input. current track. before it is moved.
        record /epmove/  Move      ! input. containes info of moved track.
c                                ! output. You have to set Move.Abort
c         Move.Abort = 0 if you want to continue the simulation of
c                           this event normally.
c                    = 3 if you want to discard this particle
c                        For example, you may use this if you use
c                        formula for further development due to this
c                        particle 
c                    = 1 if you want to abort the generation of 
c                           this event, but want to execute ue1ev.
c                    = 2 if you want to abort the genration of 
c                         this event, and skip ue1ev.
        record /epmedia/media
      real*8 r, dr
      integer idxr
      data dr/0.1d0/
      real*8 x,y 
      save

      if(aTrack.cn .eq.0) then
         write(0,*) ' error --------------------'
         stop
      endif




      Move.Abort = 0

c          You can use the following info.
c   ================================================================
c     aTrack.p         : the same particle record as Cosmos, i.e.,
c     aTrack.p.fm.p(1), p(2), p(3): momentum of the particle in GeV
c     aTrack.p.fm.p(4)            : total energy of    //
c     aTrack.p.mass               : mass  //                 
c     aTrack.p.code               : particle code 
c     aTrack.p.subcode            : particle subcode 
c     aTrack.p.charge             : charge  //.
c     aTrack.pos.x, y, z          : local coordinate of the particle in cm
c                                   If you need world coordinate, use
c                          call epl2w(aTrack.cn, aTrack.pos, wpos) 
c                    where
c                   record /epPos/ wpos is the output  world coordinate
c     aTrack.w.x, y, z            : direction cosines in the local coordinate
c                                   If you need world coordinate values, use
c                          call epl2wd(aTrack.cn, aTrack.w, ww)     
c                    where
c                   record /epDirec/ ww is the output  values in world coord.
c
c     aTrack.t                    : time in sum of (cm/beta) 
c     aTrack.wgt                  : weight of the particle when thinning is
c                                   done (normall 1.0)
c     aTrack.cn                   : component number
c
c     Move.Track      :   Info of moved track.
c     Move.boundary   :   Component boundary position when the particle crosses
c                         the boundary. (Move.Cross =T)
c                         Move.boundary.x, y, z. in cm. 
c     Move.proc       :   If the particle is going to interact, the process
c                            name is set. proc= one of
c                        'brem': bremstralhng
c                        'knoc': knock on (bhabha or moller scattering)
c                                other heavy particle knock-on
c                        'anih': positron anihilation
c                        'pair' : pair creation
c                        'comp' : compton scattering
c                        'phot' : photo electric effect
c                        'photop': photo-hadron production
c                        'coll'  : hadron's nuclear interaction
c                        'decay' : decay 
c                        '    '  : no interaction yet.
c    Move.dl          :  Path length in cm
c    Move.dE          :  Energy loss during dl. in GeV
c    Move.dEeff       :  Effective energy loss (use this for counting) GeV
c    Move.dEion       :  Energy loss due to ionization loss (GeV)
c    Move.dx          :  path length in g/cm^2
c    Move.dt          :  path length in r.l
c    Move.Cross       :  Becomes T if the ptcl crosses the boundary.
c    Move.Trunc       :  Becomes T if the ptcl track is truncated.
c    Move.Abort       :  must be set by the user.
c    
c    If you want to use a particle code, like,
c        if(aTrack.p.code .eq.  kphoton) ...
c    You have to uncomment Zcode.h above; they are the same as
c    Cosmos. Some of popular ones are:
c         kphoton:  photon             kelec  :  electron
c         kmuon  :  muon               kpion  :  pion
c         kkaon  :  kaon               knuc   :  nucleon
c         kneue  :  electron neutrino  kneumu :  muon neutrino
c         kalfa  :  He                 klibe  :  LiBeB group
c         kcno   :  CNO   group        khvy   :  Na/Mg/Si group
c         kvhvy  :  S/Cl/Ar group      kiron  :  Fe group
c         kdmes  :  D  meson           keta   :  eta  meson
c    Subcode may sometimes be needed:
c         regptcl:  particle           antip  : anti-particle
c         k0s    :  k0short            k0l    : k0 long
c 
        end
c       ************************************************
c       *  This is called when a particle passes the boundary
c       *  of a component for which you have specifed to
c       *  count particle nubmer.
c
        subroutine userbd(info, aTrack, Move, media) 
        implicit none

#include "ZepTrack.h"
#include "Zmove.h"
#include "Zmedia.h"
#include "Zcode.h"
c          next is needed if you want to judge particle code.
c  #include "Zcode.h"     

        integer info              ! 0--> ptcl is exiting to void
                                  ! <0 --> ptcl is exiting to |info| comp.
                                  ! >0 --> ptcl is entering from info comp.
        record /epTrack/aTrack    ! input. Current track before it is moved.
                                  !  If info <=0, the track position is 
                                  !  somewhere inside the component from which
                                  !  the track is existing.
                                  !  If info >0, aTrack.pos is the position 
                                  !  just before exiting the prvious component
        record /epmove/ Move        ! input/output. 
                                  ! Move.Track is the track
                                  ! infomation of the current particle
                                  ! moved to a new position.
                                  ! Say, Move.Track.cn is the  current
                                  ! comp. number.  For other details,
                                  ! see userde.

c
c                      comp.1    comp.2
c          info<0    |  *-----x|          |        
c                                                   * is aTrack.pos
c                                                   x is Move.Track.pos
c          info>0    |        *|x         |       
c  
        record /epmedia/media       ! input.

        Move.Abort = 0

        end
c       **************************************************
c       *  This is called when all the system level "end process"
c       *  for the event has been ended.  Do your own end process
c       *  for the event.
c       *
        subroutine ue1ev
        use modtemp
        implicit none
#include "Zmedia.h"
#include "ZepPos.h"
#include "Zptcl.h"
        record /epmedia/ firstM
        record /epPos/   firstP
        integer i, j, icon
        character*8 proc

        integer nptcl
        parameter (nptcl=1000)
        record /ptcl/ ptcls(nptcl)
        integer xyf   ! side face # for BGO (smaller #)
        integer nevent

        real(4)::avew, tracedp

        real(4)::edepo, ing(5), total
!              6 is the max surface numbers, i.e box in this case
!              1,2,.. corresponds to surface # 1, 2...
!              2 below is for scintillation and Cerenkov light
        real(4)::pc(2,6), pct(6)
        integer na  ! number of vol. attributes.  box-->3
        real(8):: vol(6)    ! 6 >= na
        real(4)::deT, deEff

        call epqevn(nevent)
        call epqFirstM(firstM)
        call epqFirstI(firstP)
        call epqFirstP(proc)
        write(*,*) ' event #=', nevent
        do i = 1, Ncomp
           call epqLightSensor(i, edepo,  ing,  total, icon )
           if(icon ==0 .and. total > 0) then
              write(0,'(i7, 1p,5g12.3)')
     *         i, edepo, ing(1:2), ing(5), total
           endif
        enddo

        do i = 1, Ncomp
           call epqLightPC(i, 0,  pc,  pct, icon)   ! 0 is for exiting 
                                                    ! 1 for entering
           if(icon == 0  .and. sum(pct(:)) > 0.) then
                 ! we don't show count at the side wall;  for y-BGO  
                 !  surface # 3 and 4  for x-BGO surface # 2 and 5
              !   
              call epqvolatr(i, na, vol)   !  get vol. attributes.  a,b,c of box
              if( vol(1) < 5.d0 ) then
                 ! boxa <  5cm so  this should be y-BGO
                 xyf = 3
              else
                 xyf = 2
              endif
              write(*, '(i6,1p, 3g12.3, 2x, 3g12.3)') 
     *             i, pc(:,xyf), pct(xyf), pc(:,7-xyf), pct(7-xyf)
           endif
        enddo
        do i = 1, Ncomp
           call epqEloss(i, deT, deEff)
           if( deT > 0. ) then 
              write(*, '(i7, 1p, 2g12.4)')  i, deT, deEff
           endif
        enddo
        call epqLightAveW(avew, tracedp)
        write(*,*) ' <w>', avew, ' traced # ',tracedp
        write(*,*) 
        end
c       *************************************************
c       *  This is called when all the system level "end process"
c       * for all the events has been  ended. Do your own end
c       * process for the events

      subroutine ueaev
      implicit none
         
      end

