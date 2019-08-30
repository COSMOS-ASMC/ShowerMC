#include "ZsubstRec.h"
!            treat interaction of MovedTrack
!
      subroutine cinteraction
!     use modXsecMedia
      use modSetIntInf
      use modEMcontrol
      use modColInfo
      implicit none

#include  "Zcode.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zevhnv.h"
#include  "Zincidentv.h"
#include  "Zmass.h"
#include  "Zpwork.h"

      integer i
      real*8 Ein1, Ein2, Eout, dEabs1, deabs2, dErel1, dErel2
      real*8  dErel, dEabs, Ein

!
!          used to judge if user hook should be called 
!          after MovedTrack interacted.
!
      integer neverNEP/0/,  neverE/0/,  neverG/0/
      save  neverNEP, neverE, neverG
      integer never
      integer stackpos
      save never
      integer icon
      integer,parameter::loopmax=100
      integer:: loopc



!c      record /ptcl/ fragA(maxHeavyMassN),  nonIntNucA(maxHeavyMassN)
!c      integer noOfFrag, noOfNonIntN
!
!      record /track/aTrack
!
!   **  ptcl stacking is done in each subroutine; should be changed
!           (except for hadronic interactions)
!       logic of employing events satisfying some conditions only 
!
!        generate 1 event ; data is in PWork(:)
!        call chookNEPInt(neverNEP)  ! the user may give 5 to neverNEP
!                                  ! to discard the current event
!        never = neverNEP
!        
!           
!

      loopc = 0
      do while (loopc < loopmax )
!          try until desired event is generated
!          (mainly for multiple production)
         loopc = loopc + 1
         Nproduced = 0

!
!c    assume MovedTrack is not changed in interaction routine
!     but may be added asflag by chookEing if it is elec.
!      (in the case of Job ='newskel' )
!     
!
!     cinteElec, cintePhoton, cinteNEP  will store generated particles in
!     Pwork(i), i=0,* 1, Nproduced  (in modColInfo).  Particle  momentum
!     is given in the sysemt where the projectile  (=MovedTrack) is defined.
!     When they are stored in the stack area,  track infomation (such as
!     direction coseindes) are computed.
!
         if(MovedTrack%p%code .eq. kphoton) then
!               call cintePhoton(MovedTrack%p)
            call epinteG( MovedTrack%p )
         elseif(MovedTrack%p%code .eq. kelec) then
!     call cinteElec(MovedTrack%p)
            call epinteElec( MovedTrack%p )
         else
            call cinteNEP( MovedTrack%p ) 
         endif   
         if(IntInfArray(ProcessNo)%process .eq. 'coll') then
               MovedTrack%pos%colheight = MovedTrack%pos%height
         endif

         if(MovedTrack%p%code .eq. kelec) then
            if(neverE .ne. 1) then
               call chookEInt(neverE)
               never = neverE
            endif
         elseif(MovedTrack%p%code .eq. kphoton) then
            if(neverG  .ne. 1)  then
               call chookGInt(neverG)
               never = neverG
            endif
         else
            if(neverNEP .ne. 1) then
               call chookNEPInt(neverNEP)
               never = neverNEP
            endif
         endif
         if( never /= 5) exit
         never = 0
      enddo   ! end of while;   
!     if(btest(Eabsorb(1), BitEconsv-1) ) then
      if(btest(Eabsorb, BitEconsv-1) ) then      
         if(IntInfArray(ProcessNo)%process .eq. 'coll' .or.
     *        IntInfArray(ProcessNo)%process .eq. 'photop' .or.
     *        IntInfArray(ProcessNo)%process .eq. 'munuci' ) then

!                                    last one not used yet
            call chookEabsorbC( MovedTrack, Nproduced,  Pwork, 0)
!               what is being done below is almost the same as
!               done above call.
            Ein1 = MovedTrack%p%fm%p(4)  
     *           + masn*(TargetNucleonNo-TargetProtonNo) +
     *           masp*TargetProtonNo
            Ein2 = MovedTrack%p%fm%p(4) + masp
            Eout = 0.
            do i = 1, Nproduced
               Eout = Eout + Pwork(i)%fm%p(4)
            enddo

            dEabs1 = Eout- Ein1
            dErel1 = Eout/Ein1 -1.0
            dEabs2 = Eout- Ein2
            dErel2 = Eout/Ein2 -1.0
            if( abs(dErel1) .lt. abs(dErel2)) then
               dErel=dErel1
               dEabs =dEabs1
               Ein = Ein1
            else
               dErel=dErel2
               dEabs =dEabs2
               Ein = Ein2
!                 no mass case in Eout
               Ein1 = MovedTrack%p%fm%p(4)
               dEabs1 = Eout- Ein1
               dErel1 = Eout/Ein1 -1.0
               if(abs(dEabs1) .lt. abs(dEabs) ) then
                  dErel = dErel1
                  Ein = Ein1
                  dEabs = dEabs1
               endif
            endif
            if( abs(dErel) .gt.  0.1 .or.
     *           abs(dEabs) .gt.  1.e5 ) then
               write(0,*) " code=",MovedTrack%p%code
               write(0,*) " chg=",MovedTrack%p%charge
               write(0,*) ' Moved E=',MovedTrack%p%fm%p(4)
               write(0,*) ' Ein =', Ein, ' Eout=',Eout
               write(0,*) ' Rerr =', Eout/Ein -1.0
               write(0,*) ' dEabscol= ',dEabs
               write(0,*) 'ActiveModel=', ActiveMdl
            endif
         else
!              possible process; compton, mscat, bscat,
!              anihi, decay, photoe, brems, pair cohs
            Ein = MovedTrack%p%fm%p(4)
            if(IntInfArray(ProcessNo)%process .ne. 'decay' .and.
     *           IntInfArray(ProcessNo)%process .ne. 'brems' .and.
     *           IntInfArray(ProcessNo)%process .ne. 'pair'  .and.
     *           IntInfArray(ProcessNo)%process .ne. 'cohs' )  then
               Ein = Ein + masele
            endif
            if(Ein .gt. MovedTrack%p%mass) then
               Eout = 0.
               do i = 1, Nproduced
                  Eout = Eout + Pwork(i)%fm%p(4)
               enddo
               dEabs = Eout- Ein
               dErel = Eout/Ein -1.0
               if( abs(dErel) .gt. 0.2 )  then
!            if( abs(dEabs) .gt. 1.e5) then
                  if( abs(dEabs) .gt. 1.e5 ) then
                     write(0,*) '****************************'
                  else
                     write(0,*) '----------------------------'
                  endif
                  write(0,*) 'proc=', IntInfArray(ProcessNo)%process
                  write(0,*) 'code=',MovedTrack%p%code, ' charge=',
     *                 MovedTrack%p%charge, ' E=',Ein
                  write(0,*) 'dEabs= ', dEabs, dErel, Nproduced
                  do i = 1, Nproduced
                     write(0,*) i, Pwork(i)%code, Pwork(i)%fm%p(4)
                  enddo
               endif
            endif
         endif
      endif
!///////////////////////      

      if(OneDim .eq. 0) then
!            3 dimensional
!                 stack the leading ptcl  first (to save stack area)
         call cmovePtcl3(MovedTrack, Pwork, Nproduced, Nstacked)
      else
         MovedTrack%vec = IncidentCopy%vec
         call cmovePtcl1(MovedTrack, Pwork, Nproduced, Nstacked)
      endif


      if(never .eq. 0 .or. never .eq. 1 ) then
!          user may set never=3 
      elseif(never .eq. 3) then
!          don't follow this and  child; reset stackpos
         call cgetCurrentStackpos(stackpos)
         stackpos=stackpos-Nstacked
         call cresetStackpos(stackpos)
      elseif(never .eq. 4) then
!             discard this event generated by the current primary
!             clear stack
         call cinitStack
      else
         call cerrorMsg('return value from chookE,G,NEPInt wrong', 1)
         write(0,*)  ' never=', never
         stop
      endif

      end
!          following is remnant of never=2 and ad-hoc model
!          when chookNEPInt is called before push is called.
!      The reason that we put chookNEPInt interface after push
!      is to have easy interface for skeleton making
!        
!      elseif(never .eq. 2) then
!              save only fragments and non interacting nucleons
!         if(MovedTrack.p.code .eq. kgnuc ) then
!c               get fragment and non interacting nuc.
!            call cqHvyIntF(fragA, noOfFrag)
!            call cqHvyIntNIN(nonIntNucA, noOfNonIntN)
!c
!            if(OneDim .eq. 0) then
!                call cmovePtcl3(MovedTrack, fragA, noOfFrag)
!                call cmovePtcl3(MovedTrack, nonIntNucA, noOfNonIntN)
!             else
!                MovedTrack.vec = IncidentCopy.vec
!                call cmovePtcl1(MovedTrack, fragA, noOfFrag)
!                call cmovePtcl1(MovedTrack, nonIntNucA, noOfNonIntN)
!             endif
!          endif
!
!     ************************************
!          move partcles in a given array to stack
!          3 dimensional case.

      subroutine  cmovePtcl3(iTrack, pw, n, npush)
      implicit none
!                put n ptcls in pw into stack.
!          if ThinSampling, 
#include  "Zcode.h"
#include  "Ztrackp.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"            
#include  "Ztrackv.h"
!
      integer,intent(in):: n
      type(track),intent(in)::iTrack  ! input. incident ptcl.
      type(ptcl),intent(in)::pw(n)

      integer,intent(out):: npush ! output. actual number of ptcls put in stack.
                   ! in case of ThinSampling, this may be <= Nproduced. 

      type(track)::aTrack  
      integer i
      integer loc1
      integer nact


      aTrack = iTrack
      npush = 0
      call cgetCurrentStackPos(loc1)  ! upto loc1 is already filled 

!     do i =  n, 1, -1          ! move leading ptcl first (< V9.0)
      do i = 1, n   ! now Pwork is from higher energy(except had int)
#ifdef SUBSTREC
         aTrack%p = pw(i)
#else
         aTrack%p%fm%p = pw(i)%fm%p
         aTrack%p%mass = pw(i)%mass
         aTrack%p%code = pw(i)%code
         aTrack%p%subcode = pw(i)%subcode
         aTrack%p%charge = pw(i)%charge
#endif

!               reset direction cos and related stuffs
         call cresetDirec(aTrack)


#if LABELING == 1
!
!              whennever an interaction occur,  update labelcounter
!             if info>0, clear the timer and infor counters.
!
         if(aTrack%info .gt. 0) then
!cc               aTrack.info = 0
!cc               aTrack.t = 0.
         endif
         Labelcounter = Labelcounter + 1
         aTrack%label = Labelcounter
#elif LABELING == 2
!            the above simple counter may be replaced by the
!            next sophisticated one. 
!            the same one is in the 1dim mode move routine below.
!
         call ctrickycount(iTrack, aTrack, pw, i)
!
#endif
         npush = npush + 1
         call cpush(aTrack)
      enddo 
      if(ThinSampling .and. npush .ge. 2 ) then
!         if(ThinSampling .and. npush .ge. 2  .and.
!     *   IntInfArray(ProcessNo).process .ne. 'photop') then
         call cthinStack(loc1+1, npush, iTrack, nact)
!               among npush from loc1+1, nact is accepted
!               adjust npush
         npush = nact
      endif
      end
      subroutine cthinStack(stackloc, n, iTrack, nout)
      implicit none
#include  "Zcode.h"
#include  "Ztrackp.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"            
#include  "Ztrackv.h"
#include  "Zstackv.h"
! 
      integer stackloc  ! first loc of stack where tracks of current int.
      integer n
      type(track)::iTrack  ! input. incident ptcl. of the interaction
      integer nout

      call cthinning(stack(stackloc), n, iTrack, nout)
      call cresetStackPos(stackloc-1+nout)
      end




#if LABELING == 2
!       ************************************
!             This routine updates label counters.
!          but for the survival particle form brems and knock-on
!          the label counter is not updated. 
!          For those ptcl with info > 0, timer and info counters
!          are cleared. 
      subroutine ctrickycount(iTrack, aTrack, pw, i)
      use modSetIntInf
      implicit none
#include  "Zcode.h"
#include  "Ztrackp.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"            
#include  "Ztrackv.h"
!
      integer n
      type(track)::iTrack  ! input. incident ptcl.
      type(track)::aTrack  ! input/output. i-th  secondary track 
      type(ptcl)::pw(*)    ! input.  secondary pool
      integer i              ! input.  i-th secondary (index)

      logical reset

      reset = .true.

      if( iTrack%pos%height .gt. 40.d3 ) then
         if(IntInfArray(ProcessNo)%process .eq. 'brems' ) then
!             for brem, electron has the same label as the incident
            if( aTrack%p%fm%p(4) .lt. iTrack%p%fm%p(4)*0.8) then
               aTrack%label = iTrack%label
               reset = .false.
            endif
         elseif(IntInfArray(ProcessNo)%process .eq. 'knock' .or.
     *        IntInfArray(ProcessNo)%process .eq. 'mscat' .or.
     *         IntInfArray(ProcessNo)%process .eq. 'bscat' ) then
!                for the knock-on, survival particle has the same
!                label.
            if(iTrack%p%code .ne. kelec  .or.
     *           iTrack%p%charge .ne. -1) then
!                     knockon by p,mu,pi..e+ (not by e-)
               if(aTrack%p%code .ne. kelec) then
!                       survival one has the same label
                  aTrack%label = iTrack%label
                  reset = .false.
               endif
            else
!                  electron; make the higher one has the same label
               if(pw(1)%fm%p(4) .gt. pw(2)%fm%p(4)) then
                  if(i .eq. 1) then
                     aTrack%label = iTrack%label
                     reset = .false.
                  endif
               else
                  if(i .eq. 2) then
                     aTrack%label = iTrack%label
                     reset = .false.
                  endif
               endif
            endif                     
         endif
      endif
      if(reset) then
         Labelcounter = Labelcounter +1
         aTrack%label = Labelcounter
!c         aTrack.t = 0.          ! timer reset
!c         aTrack.info = 0        ! cross counter reset
      endif
      end
#endif
!     ************************************
!          move partcles in a given array to stack
!          1 dimensional case.

      subroutine  cmovePtcl1(iTrack, pw, n, npush)
      implicit none

#include  "Zcode.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"                  
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zincidentv.h"


      integer n
      type(track)::iTrack   ! input. incident track
      type(ptcl)::pw(n)
      integer  npush   !  output.  actuall number of ptcls put in stack.
      real*8 temp, p

      integer i
      
      type(track)::aTrack 
      integer loc1
      integer nact

      aTrack = iTrack

      call cgetCurrentStackPos(loc1)

!      do i =  n, 1, -1          ! move leading ptcl last (< v9.0)
      do i =  1, n         !   now ok v9.0
#ifdef SUBSTREC
         aTrack%p = pw(i)
#else
         aTrack%p%fm%p = pw(i)%fm%p
         aTrack%p%mass = pw(i)%mass
         aTrack%p%code = pw(i)%code
         aTrack%p%subcode = pw(i)%subcode
         aTrack%p%charge = pw(i)%charge
#endif
!            see if angle of particle is larger than a lmit
         call cscalerProd(aTrack%p%fm%p, DcAtObsXyz, temp)
         call cpxyzp(aTrack%p%fm, p)
         if(p .gt. 0.) then
            temp = temp/p
         else
            temp = 1.
         endif
         if(temp .gt. BackAngLimit)  then
!              only take some limitted angle particles
!            call cresetMom(aTrack)  which is
            aTrack%p%fm%p(1) = p * aTrack%vec%w%r(1)
            aTrack%p%fm%p(2) = p * aTrack%vec%w%r(2)
            aTrack%p%fm%p(3) = p * aTrack%vec%w%r(3)

            call cgetZenith(aTrack, aTrack%vec%coszenith)


#if LABELING == 1
!                  whennever secondary particles are generated,
!               each of them get an updated label cocunter
!               if the particle has crossed the highest level
!               (info > 0),  timer and info counter is cleared
!
            if( aTrack%info .gt. 0) then
!cc                  aTrack.info = 0
!c                  aTrack.t = 0.
            endif
            Labelcounter = Labelcounter + 1
            aTrack%label = Labelcounter
#elif LABELING == 2
!                   this may be used if a tricky count is needed
!               in stead of above  counting
            call ctrickycount(iTrack, aTrack, pw, i)
#endif
            npush = npush + 1  
            call cpush(aTrack)
         endif
      enddo
      if(ThinSampling .and. npush .gt. 0 ) then
         call cthinStack(loc1+1, npush, iTrack, nact)
         npush = nact
      endif
      end

!     ****************************************************************
      subroutine cqIntePtcl(ptclA, num)
      implicit none
!          inquire the particle information that made interactions
!        to produce secondary particles.
!        If "MovedTrack" is a heavy,  ptclA will get interacting nucleons
!        otherwise, ptclA will have MovedTrack.p itself.
!
!           
#include  "Zcode.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
! #include  "Ztrackp.h"
#include  "Ztrackv.h"

      type(ptcl)::ptclA(*)   ! output. interacted particles. max size
                               ! should be maxHeavyMassN (= 56 =Fe)
      integer num              ! output. number of ptcls in ptclA
!
!
!c      if(MovedTrack.p.code .ge. kdeut .and. 
!c     *    MovedTrack.p.code .le. khvymax ) then
      if( MovedTrack%p%code .eq. kgnuc ) then
         call cqHvyIntIN(ptclA, num)
      else
         num = 1
#ifdef SUBSTREC
         ptclA(1) = MovedTrack%p
#else
         ptclA(1)%fm%p = MovedTrack%p%fm%p
         ptclA(1)%mass = MovedTrack%p%mass
         ptclA(1)%code = MovedTrack%p%code
         ptclA(1)%subcode = MovedTrack%p%subcode
         ptclA(1)%charge = MovedTrack%p%charge
#endif
      endif
      end




