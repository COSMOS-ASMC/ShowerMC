      subroutine epUI( info, loc1, loc2)
! Basic purpose of this routine is do something when a particle makes
! an interaction. 
      implicit none
#include "ZepTrackv.h"
#include "Zcode.h"
      integer,intent(inout):: info  !  Basically, the particle code
                ! of the particle which made an interaction. For rare
                ! particles, they may be given a some unique value
                ! diff. from acutal code (and branch is made to
                ! "case default" )
                ! Exact code may be obtained by seeing Move.track.p.code.  

                ! This routine is called whennever
                ! a particle interacts, i.e, just after the interaction.
                ! The interaction products  are stacked in the particle
                ! stacking area.  The user can access to the stack
                ! by calling a subroutine with a stack index.
                ! The product is stacked with index from loc1 to loc2.
                ! In the case of  brems, 
                ! a gamma and the scattered  electron are the products.
                ! 
                ! When this is called, info is !=0.
                ! If the user dose not give any value on return, 
                ! this routine will not be called ***for that particle***
                ! again, until the next event  simulation starts. 
                ! The default is default to the letter.  If you require
                ! this routine be called for that particle  again within 
                ! the current event simulation, make info = 0 on return.
                !
      integer,intent(in)::loc1,loc2   ! see above
!                    To get information of a stacked track, use
!              call  epgetTrack(i, aTrack, icon) 
!          where  i must be loc1<= i <= loc2.  
!                 record /epTrack/ aTrack is the ouput
!                 integer icon == 0 ==> i is valid
!                         icon != 0 ==> i is invalid.
!

!  The user can access following variables for the interaction point.

!   Move.proc: the interaction type. see each "case" below
!          Cn:  component number  
!   Move.track:  interacting particle info. e.g  Move.track.pos.x etc
!          is the  position info. in local coordinate of Cn.
!          to convert it to world  coordinate,
!          call epl2w(Cn, Move.track.pos, posw)  for position
!          call epl2wd(Cn, dir,  dirtemp)        for direction cos

!       Move.track.p.code : code  
!       Move.track.p.fm.p(1:4): 4 momentum        

!   Media(MediaNo): media info. E.g  Media(MediaNo).name is the media name
!              (see Epics/epics/Zmedia.h) 
!   Media(MediaNo).colA (also colZ, colXs): target of nuclear interaction.
!                (A,Z) and cross-section (in mb)
!   a(i).p.code etc can be used to know produced particle properties.
!          (i=1,n).
!
      record /epTrack/ aTrack
      integer icon, i

!      write(0,*) ' called with info=',info, loc1,loc2, Move.proc
!      write(0,*) 'Cn, code sub charge=',Cn, Move.track.p.code,
!     *  Move.track.p.subcode,Move.track.p.charge
!      write(0,*) 'Energy=',Move.track.p.fm.p(4)

      select case(info)
!             you may freely modify "case" ; 
!             say, case(kpion:kgnuc) instead of case(kpion:kkaon)...
!             etc below
         case(kphoton)  ! gamma interaction
            ! do something for gamma ineraction
            ! if you want this routine be called again for this
            ! event,  make info =0 else don't touch it. Then,
            ! this routine will not be called until next event
            ! simulation starts.  
            !  Possilbe Move.proc values are
            !     comp : Compton scattering
            !     pair : pair creation
            !     phot : photoelectric effect
            !     coh  : coherent scattering
            !     photop : photo-production of hadrons
            !     mpair : magnetic pair creation

         case(kelec)    ! electron interaction
             !   Move.proc
             ! brem : Bremstrahlung   
             ! knoc : knockon (Moller or Bhabka)
             ! anih : positron annihilation
             ! sync : synchroton emission
             ! Note-- Cerenkov light emission happens along the 
                     ! path of a charged particle and dose not
                     ! call this routine
         case(kmuon)     ! muon interaction
              ! knoc: knock-on
              ! decay:   decay   may include negative muon capture
              ! pair
              ! brem
              ! nuci : nuclear interaction
         case(kpion:kkaon)  ! pi, K
              ! knoc
              ! decay
              ! coll   collision
             
         case(knuc)       ! p,n
              ! knoc
              ! coll :  anti-prooton /anti-neutron annihilation  included

              ! To wait until nuclear collision takes place
!            if( Move.proc /= 'coll' ) then
!               info = 0
!            else
!               do i = loc1, loc2   ! print stacked products
!                  call epgetTrack(i, aTrack, icon) 
!                  write(0,*) i, aTrack.p.code
!              enddo
!            endif

         case(kgnuc)     ! heavy int
              ! knoc
              ! coll :  
!         case(klight)  !  light interaction.   This is probably
                        !  nonsense for the moment
          
         case default   ! others.  mainly eta decay
              ! knoc, coll, decay.    
      end select
      end

      subroutine epGUI(info, sinfo)
!           This is for future extension. 
      use modGUI
      implicit none
!  #include "ZepTrackv.h"
!  #include "Zcode.h"
      integer,intent(in):: info
      type(gui),intent(inout)::sinfo
      end
