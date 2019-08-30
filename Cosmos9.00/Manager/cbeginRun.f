!        Initialize simulation.  
#if  defined (KEKB) || defined (KEKA)
#define DOMPI
#endif
!
      subroutine cbeginRun
      use modEMcontrol
!      use modColInfo
      implicit none
#include  "ZmediaLoft.h"
#include  "Zmanagerp.h"
#include  "Zmanager.h"
#include  "Zelemagp.h"
#include  "Zevhnp.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zincidentp.h" 
! #include  "Zprimary.h"
! #include  "Zprimaryv.h"
#include  "Zcondc.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#if defined (DOMPI)
#include  "mpif.h"
#include  "Zmpi.h"
      integer intdata
#endif
      character*16 temp
      integer jold, icon
      real*8 s1
      integer:: i
!!!!!!!!!!!!
! after reading input namleist parameter, program flow
!     comes here.  the work dependent on LatObsSite and
!     LongitObsSite  must not be done before ciniAtmos,
!     if the NRL atmosphere model is to be used by specifing
!     AtmosFile="...." in which case, (lat, long) is given
!     in that file. So (lat,long) pair is forced to replace
!     LatitOfSite and LongitObsSite values.
      
      if(Cont) then
!           restore status at the end of previous run
         call crestoreStatus
         Cont =.true.
      elseif(Job .eq. 'flesh'  .or. Job .eq. 'newflesh') then
         temp = Job
         call cgetSkelFile 
         Job = temp    !  keet as now
      elseif(Job .eq. 'skeleton' .or. Job .eq. 'newskel' ) then
!            to be safe, parameters are written in a file first and read
!            again ; this will avoid the possible difference of internal
!            parameter values read from original param and SkeletonParam
!
#ifndef DOMPI
         temp = Job
         call copenNLfw(TempDev, SkeletonFile, icon)
         if(icon .ne. 0) then
            call cerrorMsg(SkeletonFile, 1)
            call cerrorMsg('File shown above cannot be opened',0)
         endif
!           at flesh time, you don't need to rewrite Job
         if(Job .eq. 'newskel') then
            Job = 'newflesh'
         else
            Job = 'flesh'
         endif
         call cwriteParam(TempDev, 1)
         close(TempDev)
!           macIFC/MACOSX cannot put correct namelist char data 
!           (' is missing) so cannot read SkeletonFile
!           avoid read it. the user must do the equiv.  by hand 
!            call copenNLf(TempDev, SkeletonFile, icon)
!            call creadParam(TempDev)
!            close(TempDev)
!            restore Job
         Job = temp
#endif
      endif
!        reset parameters
      if(BorderHeightH .eq. 0.0) then
         BorderHeightH = HeightOfInj + 1.0d0    ! from v6.10
      endif

      EventsInTheRun = 0      ! moved from ciniTracking0. 2001.09/30

      if(DestEventNo(2) .eq. 0) then
         DestEventNo(2) = DestEventNo(1)
      elseif(DestEventNo(1) .lt. 0) then
         DestEventNo(2) = -abs(DestEventNo(2))
      endif

!         to lower case letters.
      temp = Generate
      call c2lowerCase(temp, Generate)
      temp = Job
      call c2lowerCase(temp, Job)

      call rndsw(jold, 1)       ! specify random number generator 1.

!  Now (v7.651) next random number init is placed before cixsec
!     otherwise qgsjetII-03 generates the same event every run,
!     (ciniRN was placed inside cwhatJob bef. v7.651)
      call ciniRN
!         until  >>>   moved after ciniAtmos
!           set Target and do related xsec business; moved before
!           cintModels (for qgsjet dummy ptcl generation in cQGSini)
!!      call cixsec
!!      call cHowMCS  ! MCSmodel is analysed ...
!!      call ciniMCS  ! init for MCS and set doNewMCS (t/f)
!! >>>>>
!       cintModels probably reads disk file(s) internally, 
!       we may better to avoid simultaneous access to that disk 
!       by different ranks; reading  starts from rank0, rank1...
#if defined DOMPI
      if(mpirank .eq.  0) then
         call cintModels('cosmos') ! analysis of interaction models and init.
         if(mpisize .gt. 1) then
            call MPI_SEND(mpirank, 1, MPI_INTEGER, 1, 1,
     *           MPI_COMM_WORLD, mpierr)
         endif
      else
         call MPI_RECV(intdata, 1, MPI_INTEGER, mpirank-1, 1,
     *      MPI_COMM_WORLD, mpistat, mpierr)
         call cintModels('cosmos') ! analysis of interaction models and init.
         if( mpirank .lt. mpisize-1 ) then
            call MPI_SEND(mpirank, 1, MPI_INTEGER, mpirank+1, 1,
     *       MPI_COMM_WORLD, mpierr)
         endif
      endif
#else
      call cintModels('cosmos') ! analysis of interaction models and init.
#endif

!     init for Atmosphere:  This includes, as Atmosphere Media,
!     everything except knockon  will be performed
!
      call ciniAtmos
      call ccheckAtmosModel
!           set Target and do related xsec business; moved before
!           cintModels (for qgsjet dummy ptcl generation in cQGSini)
!!      call cixsec  ! no more needed
      call cHowMCS  ! MCSmodel is analysed ...
      call ciniMCS  ! init for MCS and set doNewMCS (t/f)
!!!!!!!!!
!      call cwriteParam(0, 1)
!!!!!!!!!      
!         init for geomag
      call crdGeomag(GeomagFile, YearOfGeomag)
!     init for LPM effect energy sampling.
!       Not needed now v9.0.  csetLPMcnst is called inside
!     epBrHSampE.  consts depend on media. The 3rd argument is
!     Emin/Ee and dependent on  Ee  if we fix Emin, but we
!     fix it (vmin=Emin/Ee) to be ~ 1e-5, since  small Eg brems is
!     strongly suppressed at very high energies, we can safely
!     neglect low energy gamma emission. So the const setting
!     may be done independently of the Ee and can be called
!     when we compute media dependent const. But we put it
!     in the Eg sampling routine (to prepare for varying vmin,
!     with fixed Emin).    Also csetLPMcnst is called from
!     epPrHSampE for pair cretion. with different vmin.(not
!     used for pari cre).
!       
!      s1 = (TargetAtomicN**(1./3.d0)/183.d0)**2
!      call csetLPMCnst(s1, log(s1), 1.d-4, X0)
!        initialize  observation  
      call cinitObs
!        init for   primary sampling
      call ciniSPrim(PrimaryFile)

      if(CutOffFile .ne. ' ') then
         call crigCut0(CutOffFile) ! read cutoff talbe and init.
      endif
!     init for muon interaction routines; specific for Air.
!            next two are removed; 
!         call cRdmuTab        ! set various consts for mu int.
!     call cSetMu(TargetAtomicN, TargetMassN)
!           equivalent management
!      call eprdMfile  LibLoft/PreInte/eprdMFile.f
!        call epReadTab   LibLoft/EM/FromEpi/epReadTab.f
!          call epRdmuTab  LibLoft/Mu/FromEpi/epRdmuTab.f
!            call epSetMu    LibLoft/Mu/FromEpi/epSetMu.f
!11                 FromEpics = .false.  ! muon interaction routines for Air are    not needed?
      
                           ! inside Cosmos 
      call ciniSPrimAng    !  this is in csPrimaAgn.f in Tracking dir.
!        check job
      call cwhatJob
!          this is moved here; before v6.10 it was before cwhatJob
      call ciniTracking0  !  init for  tracking for all events
!     init for knockon process
      do i = 1, NoOfMedia
         if(KnockOnRatio .lt. 1.d0) then
!              ********* Knockon is not used now ******
            if(Job .eq. 'newskel') then
               call epSetTcut( Media(i), KEminObs(1)*KnockOnRatio )
            elseif(KEminObs2(1)*KnockOnRatio .gt. 0.) then
!                 this must come after cwhatJob, since KEminObs2 must
!           be fixed.  For skeleton-flesh job, KnockOnRatio should be
!           small enough so that KEminObs2*KnockOnRatio < KEminObs at
!           flesh time  
               
               call epSetTcut(Media(i), KEminObs2(1)*KnockOnRatio)
            else
               call cerrorMsg('KnockOnRatio<1 and others mismatch', 0)
            endif
         else
!!             v7.643
            if( knockOnRatio == 1.0d0 .or. RecoilKineMinE == 0.) then
               RecoilKineMinE= KEminObs(1)
            endif
            call epSetTcut( Media(i),RecoilKineMinE ) 
         endif
      enddo
!         user hook
      call chookBgRun
      
      end

      subroutine ccheckAtmosModel
      implicit none
!      now we don't use system function.       Current model can be known via cq
      !  next is for "AtmosModel" variable
#include "Zmanagerp.h"
      integer:: modelnum
!!      character(128):: string="$COSMOSTOP/Scrpt/currentAtmos.sh"
!!      character(128):: command
!!      integer:: status
!!      integer,external:: system
!!      
!!!          env. var $COSMOSTOP is replaced by actual value
!!      call cgetfname(string, command)
!!!      write(0,*) " command: ", trim(command)
!!     see atmos model is same as given in AtmosModel
!!      status = system( command )
!!      write(0,*) ' status=', status
!!     !      if( status /= AtmosModel ) then


      call cqAtmosModel(modelnum)

      if( modelnum /= AtmosModel ) then
         write(0,*) 'Current atmosphere model # is ', modelnum
         write(0,*) 'But "AtmosModel" in namelist is ',AtmosModel
         write(0,*)
     *    "Make  them consistent by changing the value of AtmosModel"
         write(0,*) ' or by using "atmosModel.sh" command to change'
         write(0,*) " the atmosphere model"
         stop
      endif
      end  subroutine ccheckAtmosModel

      subroutine ciniTracking0
      use modEMcontrol
      implicit none
#include  "Zmanager.h"
#include  "Zmanagerp.h"
#include  "Zevhnp.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"            
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zprimary.h"
#include  "Zprimaryv.h"
#include  "Zobs.h"
#include  "Zobsv.h"
#include  "Zincidentp.h"
#include  "Zheavyv.h"

!

      integer i
       real*8 dstep 
!  
!
!
      ObserveAS = index(Generate, 'as')  .gt. 0 .or.
     *            index(Generate, 'lat') .gt. 0 
      if(index(Generate, 'qas') .gt. 0) then
         SkipPtclGen = 1  ! quick as generation for heavies.
      else
         SkipPtclGen = 0
      endif
      
!
!    d = min( dZ/2, geneal min) for e+/e- other charged.
!    general min = same as so far for e+/e-, p, heavy
!                  mu, pi, K, --> dE/Ek< 1%; dE=Ek/100.
!                  min=Ek/100/2e-3 g/cm2 = Ek/10/2e-3 kg/m2
!                     =Ek/2e-2= 50Ek 
!                  Ek=1  -->50 kg/m2  = 5g/cm2 ~5000cm ~50m
!                  Ek=0.1-->5 kg/m2=0.5g/cm2   ~ 5m
!                    0.01-->0.05kg/m2=0.05g/cm2 ~ 0.5m
!                    0.001-->                     0.5m  
!        
!             ////////                   
      do i = 1, NoOfSites
         if( i .eq. 1) then
            dstep = ObsSites(1)%pos%depth 
         else
            dstep = ObsSites(i)%pos%depth - ObsSites(i-1)%pos%depth
         endif
         if(dstep  .lt. 15. )  then
            StepControl=2
         else
            StepControl = dstep/25.0 
         endif
         maxstep(i) =   dstep/StepControl
      enddo
!         added v7.651
      maxstep(0) = maxstep(1)
      maxstep(NoOfSites+1) = maxstep(NoOfSites)

!          compute the offset point in 'xyz' system
!        the deepest detector origin + Offset is the point
!        to which the primary is directed.
!         offset in the detector system.
      Offset%r(1) = 0.
      Offset%r(2) = 0.
      Offset%r(3) = OffsetHeight 
!        convert it to xyz system.
      call cdet2xyz(ObsSites(NoOfSites)%pos%xyz, Offset, Offset)
!        make it offset 
      do i= 1, 3
         Offset%r(i) = Offset%r(i) -
     *    ObsSites(NoOfSites)%pos%xyz%r(i)
      enddo
!     next Eabsorb(1)--> Eabsorb
!          Eabsorb(2)--> EabsorbL   in modEMcontrol  from version9.00
!      if(Eabsorb(1) .ne. 0) then
!         if(Eabsorb(2) .le. 0) then
!            Eabsorb(2) = NoOfSites
!         elseif( Eabsorb(2) .gt. NoOfSites) then
!            call cerrorMsg("Eabsorb(2) > NoOfSites", 0)
!         endif
!     endif
      
      if(Eabsorb /= 0) then
         if(EabsorbL <= 0) then
            EabsorbL = NoOfSites
         elseif( EabsorbL > NoOfSites) then
            call cerrorMsg("EabsorbL > NoOfSites", 0)
         endif
      endif

      end
!     ************
      subroutine  ciniRN
!     initialize random # generator. This was a part of
!     cwhatJob before v7.651 but separated as ciniRN
!     and is called from before cixsec. 
!     ************
      implicit none
#include  "Zmanager.h"
#include  "Zmanagerp.h"
!
!
      character*190  msg
      real*8 u
      integer::i,   now(2)

!      move from ciniTracking0  
      RefreshIR = InitRN(1) .lt. 0 .and. 
     *        ( Job .ne. 'flesh' .and. Job .ne. 'newflesh')

      if(InitRN(1) .gt. 0 .and. InitRN(2) .gt. 0 ) then
         call rnd1r(InitRN)     ! init randeom number generator
!      *****************
      elseif(.not. RefreshIR .and. InitRN(2) .lt. 0) then
         call cmkSeed(0, now)   ! make seed using timer and hostname
         call rnd1r(now)
!           dummy use of 1000 times
         do i = 1, 1000
            call rndc(u)
         enddo
      endif
!          this is almost ok but later once more saved.
      call rnd1s(SeedSave)
      end
!     ******************
      subroutine  cwhatJob
!     ************
      implicit none
#include  "Zmanager.h"
#include  "Zmanagerp.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"            
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
!   
      integer klena, icon
      character*8  uid 
      character*16 temp

      integer i
      character*190  msg

      if(KEminObs(2) .ne. KEminObs(1)) then
         write(0,*) ' KEminObs(2) is forced to be the same as'
         write(0,*) ' KEminObs(1)=',KEminObs(1)
         KEminObs(2)= KEminObs(1)
      endif

      if(Job .eq. ' ' .or. Job .eq. 'skeleton' .or. 
     *   Job .eq. 'newskel'  ) then
         if(Job .ne. 'newskel') then
!             save present conditions
            do i = 1, 8
               KEminObs2(i) = KEminObs(i)
            enddo
            Generate2 = Generate
            EndLevel2 = EndLevel
         elseif(Job .eq. 'newskel') then
            if( KEminObs2(1) .ge. KEminObs(1) .and.
     *          EndLevel2 .le. EndLevel .and.
     *          index(Generate2,'as') .eq. 0 .and.
     *          index(Generate2,'lat') .eq. 0 ) then
               call cerrorMsg(
     *          'Doing newskel job seems nonsense', 1)
               call cerrorMsg(
     *          'Check Generate2, KEminObs2(1), EndLevel2',0)
            endif
         endif
         NoOfSites2 = NoOfSites    ! probably not needed
         if(Job .eq. ' ') then
            if(SeedFile .ne. ' ') then
!                  open seed file for output
               write(msg, *) 'opening SeedFile=',
     *         SeedFile(1:klena(SeedFile))
               call cerrorMsg(msg, 1)
               call copenfw(SeedFileDev, SeedFile, icon)
               if(icon .ne. 0) then
                  call cerrorMsg(SeedFile, 1)
                  call cerrorMsg('File shown above cannot be opened',0)
               endif
               if(Cont) then
                  call cskiptoEOF(SeedFileDev)
               endif
            endif               
         elseif(Job .eq. 'skeleton' .or. Job .eq. 'newskel' ) then
            call
     *      cerrorMsg('  ********** skeleton making **********', 1)
            write(msg, *)  '      Generate=', Generate
            call cerrorMsg(msg, 1)
!                 save skeleton inf. in skelotonFile file.
!             The file will be modified when the distjob command
!             processes Job = 'flesh' later. You need not modify
!             skeleton file if distjob is employed.
            if(.not. Cont) then
               temp = Job    ! save current Job
               if(Job .eq. 'newskel') then
                  Job = 'newflesh'
               else
                  Job = 'flesh'
               endif
               call copenNLfw(TempDev, SkeletonFile, icon)
               if(icon .ne. 0) then
                  call cerrorMsg(SkeletonFile, 1)
                  call cerrorMsg(
     *              'File shown above cannot be opened',0)
               endif
               call cwriteParam(TempDev, 1)
               close(TempDev)
               Job = temp
            endif

!                open SeedFile 
            if(SeedFile .eq. ' ') then
!                error. you need file; not needed actually
!               write(msg, *)
!     *         ' SeedFile must not be blank for skelton making'
!               call cerrorMsg(msg, 0)
            else
               write(msg, *) 'opening SeedFile=',
     *         SeedFile(1:klena(SeedFile))
               call cerrorMsg(msg, 1)
               call copenfw(SeedFileDev, SeedFile, icon)
               if(icon .ne. 0) then
                  call cerrorMsg(SeedFile, 1)
                  call cerrorMsg('File shown above cannot be opened',0)
               endif
               if(Cont) then
                  call cskiptoEOF(SeedFileDev)
               endif
            endif
         endif
      elseif(Job   .eq. 'flesh' .or. Job .eq. 'newflesh') then
!             don't worry about KEminObs2 etc.  They have been read
!             from &Hparam
          call cerrorMsg('  ********** fleshing job   *********', 1)
          if(Job .eq. 'flesh') then
             if(EndLevel .gt. EndLevel2) then
!                to deeper detph than skeleton
                write(msg, *)
     *          ' fleshing will be done to deeper depth than'//
     *               ' skeleton making time'
                call cerrorMsg(msg, 1)
                write(msg, *) ' No of old levels=', EndLevel2,
     *               ' No of new levels=', EndLevel
                call cerrorMsg(msg, 1)
             elseif(EndLevel .lt. EndLevel2) then
                call cerrorMsg('EndLevel must be >= skelton time', 0)
             endif
             write(msg, *)  '      Old Generate=', Generate2
             call cerrorMsg(msg, 1)
             write(msg, *)  '      New Generate=', Generate
             call cerrorMsg(msg, 1)
          else
             if(EndLevel .lt. EndLevel2) then
!                to deeper detph than skeleton
                write(msg, *)
     *          ' fleshing will be done to deeper depth than'//
     *           ' skeleton making time'
                call cerrorMsg(msg, 1)
             endif
!
!              copy  old Generate, KEminObs2 to current value
!              they are future values at newskel time
             do i = 1, 8
                KEminObs(i) = KEminObs2(i)
             enddo
             Generate = Generate2
             EndLevel = EndLevel2
          endif
!                open SeedFile
          if(SeedFile .eq. ' ') then
!c             write(*, *) ' SeedFile must not be blank for flesh job'
!c             call cerrorMsg(msg, 0)
!               seed will be read from Mdev.
          else
             write(msg, *) 'opening SeedFFile=', 
     *        SeedFile(1:klena(SeedFile))
             call cerrorMsg(msg, 1)
             call copenf(SeedFileDev, SeedFile, icon)
             if(icon .ne. 0) then
                call cerrorMsg(SeedFile, 0)
                call cerrorMsg('File shown above seems missing',0)
             endif
          endif
      else
           write(msg,*) ' Job=',Job, ' undefined'
           call cerrorMsg(msg, 0)
      endif
      if((Trace .gt. 0 .and. Trace .lt. 60) .or. Trace .gt. 100) then
!           defalut trace.  fix the dirctor
          if(TraceDir .eq. ' ') then
             call cgetLoginN(uid)
             TraceDir = '/tmp/'//uid(1:klena(uid))
          endif
       endif
      end
!        **************************************** read cont job info
      subroutine crestoreStatus
      implicit none
#include "Zmanagerp.h"
      integer icon

      call copenNLf(TempDev, ContFile,icon)
      if(icon .ne. 0) then
         call cerrorMsg(ContFile, 1)
         call cerrorMsg('File shown above seems missing',0)
      endif
      call creadParam(TempDev)
      close(TempDev)
      end
      subroutine cgetSkelFile
      implicit none
#include "Zmanagerp.h"
      integer icon
      call copenNLf(TempDev, SkeletonFile, icon)
      if(icon .ne. 0 ) then
         call cerrorMsg(SkeletonFile, 1)
         call cerrorMsg('File shown above seems missing',0)
      endif
!          read skelton parameters for flesing
      call creadParam(TempDev)
      close(TempDev)
      end


