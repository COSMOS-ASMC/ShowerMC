#include "ZsubstRec.h"
#include "Zunionmap.h"
      ! many are now obso.  should be arranged
!       ********************************************
!
!       interface routines between epics and cosmos
!       for hadronic interactons. see 
!                    ReadmeNI
!       ********************************************
!
!

       subroutine ep2cosPtcl(aPtclin)
!          this is called when a ptcl is popped up from
!          stack. (prog/Light/epLightStack.f)
!          at this moment, it is kept as aPtcl here
       implicit none
#include "Zepi2cos.h"
       type(ptcl):: aPtclin

       aPtcl = aPtclin

       end
!       *****************************
       subroutine ep2cosCond
        implicit none
#include "Zepi2cos.h"
#include "Zcode.h"
!c       real*8 tmass, zchg 
#ifdef UNIONMAP
       TrackBefMove%p = aPtcl
#else
       TrackBefMove%p%fm = aPtcl%fm
       TrackBefMove%p%mass = aPtcl%mass
       TrackBefMove%p%charge = aPtcl%charge
       TrackBefMove%p%code = aPtcl%code
       TrackBefMove%p%subcode = aPtcl%subcode
#endif

!             This is not used in m.f.p calculation from v7.
!c       TargetMassN = tmass
!c       TargetAtomicN = zchg
       FromEpics = .true.   !  no meaning now
       end
!       *********************
       subroutine ep2cosCondr
        implicit none
#include "Zepi2cos.h"
        FromEpics =.false.  ! no meaning now
        end

!      **********************************
       subroutine ep2cosCond2
       implicit none
#include "Zepi2cos.h"

#ifdef  UNIONMAP
       MovedTrack%p = aPtcl
#else
       MovedTrack%p%fm = aPtcl%fm
       MovedTrack%p%charge = aPtcl%charge
       MovedTrack%p%mass = aPtcl%mass
       MovedTrack%p%code = aPtcl%code
       MovedTrack%p%subcode = aPtcl%subcode
#endif
       end
!     ***************************
      subroutine eppushPtcl(cTrack)
!            push Cosmos made ptlcs in the stack area
      implicit none
#include "Zepi2cos.h"      
#define PTCL
#include "ZepTrack.h"
#undef PTCL
       type(epTrack)::  cTrack  ! input. projectile of the interacion
                           ! some track info is extracted
                           ! from here. 

       type(epTrack)::  nTrack 
      
      real*8  p, sump2, adj
      integer i

      nTrack%pos%x = cTrack%pos%x
      nTrack%pos%y = cTrack%pos%y
      nTrack%pos%z = cTrack%pos%z
      nTrack%t = cTrack%t
      nTrack%user = cTrack%user
      nTrack%cn = cTrack%cn
      do i = 1, Nproduced
         p = Pwork(i)%fm%p(4)**2 - Pwork(i)%mass**2
         sump2 = Pwork(i)%fm%p(1)**2
     *       +  Pwork(i)%fm%p(2)**2
     *       +  Pwork(i)%fm%p(3)**2
         if(p .le. 0. .or. sump2 == 0.) then
            p = 0.
            Pwork(i)%fm%p(4) = Pwork(i)%mass
#ifdef SUBSTREC
            nTrack%p = Pwork(i)
#else
            call epsubptcl(Pwork(i), nTrack%p)
#endif
            nTrack%w%x = 0.
            nTrack%w%y = 0.
            nTrack%w%z = 1.
         else
!            this may not be the same as p because, pi0<->pic
!            k0<->kch exchange to conserve total charge in Cosmos
!            ; because mass is changed. adjust it otherwise
!            direction cos may become inconsistent and
!            boundary error might happen
            adj = sqrt(p/sump2)
            Pwork(i)%fm%p(1) =Pwork(i)%fm%p(1)* adj
            Pwork(i)%fm%p(2) =Pwork(i)%fm%p(2)* adj
            Pwork(i)%fm%p(3) =Pwork(i)%fm%p(3)* adj
            p = sqrt(p)
            nTrack%w%x = Pwork(i)%fm%p(1)  / p
            nTrack%w%y = Pwork(i)%fm%p(2)  / p
            nTrack%w%z = Pwork(i)%fm%p(3)  / p
         endif
#ifdef SUBSTREC
         nTrack%p = Pwork(i)
#else
         call epsubptcl(Pwork(i), nTrack%p)
#endif
         call eppush(nTrack)
      enddo
      end

!         for ibm; probably obso.
      subroutine epsubptcl(inp, out)
      implicit none
#include   "Zptcl.h"
      type(ptcl)::  inp, out
      out = inp
      end

!           fix interaction type choosing the smallest path.
!           This is a modified version of cfixProc in Cosmos.
!           The difference is due to the atmosphere and other
!     constant density material
!!!      subroutine epfixProc(den, gramcm2, proc)
!!!      use modV1ry
!!!      implicit none
!!!#include  "Zglobalc.h"
!!!#include  "Ztrack.h"
!!!#include  "Ztrackv.h"
!!!      real*8  den     ! input.  denstiy g/cm^3
!!!      real*8  gramcm2 ! output. path length in g/cm2
!!!      character(*),intent(out):: proc ! output. process id characters
!!!                    ! <= 8 
!!!      real*8 len, minlen
!!!!     
!!!      integer i
!!!      
!!!      minlen = Infty
!!!
!!!      if( V1ry == 2 ) then
!!!         call epForceV1ryInt    !  make path for the target int.
!!!                   ! 0, so that next section select it. 
!!!      endif
!!!
!!!      
!!!      do i = 1, NumberOfInte
!!!         if(.not. IntInfArray(i)%decay) then
!!!!              convert kg/m2 into length in m
!!!!            IntInfArray(i).length = IntInfArray(i).thickness*0.1/den/100.
!!!            len = IntInfArray(i)%thickness*0.001d0/den
!!!         else
!!!            len = IntInfArray(i)%length
!!!         endif
!!!
!!!         if(i .eq. 1 .or.  len .lt. minlen) then
!!!               ProcessNo = i
!!!               IntInfArray(i)%length = len
!!!               minlen = len
!!!         endif
!!!      enddo
!!!
!!!      if(IntInfArray(ProcessNo)%decay) then
!!!         IntInfArray(ProcessNo)%thickness =
!!!!     *     IntInfArray(ProcessNo).length*100. * den*10. ! cm g/cm3*10-> kg/m2
!!!     *      IntInfArray(ProcessNo)%length*1000. * den   !  in kg/m2
!!!      endif
!!!      gramcm2 = IntInfArray(ProcessNo)%thickness * 0.1   ! in g/cm2
!!!      proc = IntInfArray(ProcessNo)%process
!!!      end
      subroutine epResetProcNoForV1ry
          reset ProcessNo to be for the "hadint"  of virtual 1ry.
      use modV1ry
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
      integer:: i

      do i = 1, NumberOfInte

         if( IntInfArray(i)%process == hadint ) then 
            ProcessNo = i
            return  ! *******
         endif
      enddo
      write(0,*) ' for hadoron interaction=', hadint
      write(0,*) ' not in the interaction candidates in IntInfArray'
      write(0,*)
     *  ' Virtual 1ry treatment error: msg from  epResetProcNoForV1ry'
      stop
      end
      
!     *********************************************
      subroutine epsampPtcl(aPtcl)
#include "Zptcl.h"
      integer fin
       type(ptcl)::  aPtcl  ! output. E, code, subcode, mass, charge

      call csampPrimary(aPtcl, fin)
      

      end
!          dummy routine which will never used but is needed
!          to bypass the problem of unresolved external ref.
      subroutine chookTrace
      end

      subroutine epsmpNEPIntL(media)
!             xmedia=>media is to avoid name                          
!       collision of media  in modXsecMedia and                       
!       media argument in the subroutine def. 
      use modXsecMedia, xmedia=>media
      implicit none
#include  "Zglobalc.h"
#include  "Zcode.h"
#include  "Ztrackp.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zheavyp.h"
#include  "Zelemagp.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zcmuint.h"
!     **************************************************
!
      type(xsmedia),intent(inout):: media

      real*8 mfp,  xs, mass

      real*8 collkgram, u, length

      call cdecayLeng(TrackBefMove, length)

      if(length .ne. Infty) then
         call csetIntInf(length, .true., 'decay')
      endif

      if( TrackBefMove%p%code /= kmuon)  then
!         call epgetxs(ActiveMdl2, TrackBefMove.p, media, xs, mfp)
         call cGetXsec(ActiveMdl2, TrackBefMove%p, media, xs, mfp)
         if(xs == smallxs ) then
            collkgram = Infty
         elseif(xs == largexs) then
            collkgram = 0.
         else
            call rndc(u)
            collkgram = -mfp*log(u)
         endif

         call csetIntInf(collkgram, .false., 'coll')            
      endif
      end
!      
!     next 2 sub are moved into cepSampNEPIntL.f
!      
!      subroutine epsaveFirstCol
!!        This is called only when the first nuclear interaction
!!        or decay to save the collision info.
!      implicit none
!#include "Zepi2cos.h"      
!#include "ZepMaxdef.h"
!
!
!       type(ptcl)::  ptclSave(EPMAX_STACK)
!      integer nptcls
!      common /stackSavec/ ptclSave, nptcls
!      
!      integer i
!      if( Nproduced >  EPMAX_STACK) then
!         write(0,*) "Warning, first interacton generated: "
!         write(0,*) "# of particls =",Nproduced, " > ",
!     *               EPMAX_STACK, " so "
!         write(0,*) " only first ", EPMAX_STACK, " ptcls are saved"
!         write(0,*) 
!     *     " for later retrieving (by calling epqFirstColPtcls)."
!         write(0,*) " Simulation itself are not affected. "
!         write(0,*)
!     *     " If you need all the generated particle information "
!         write(0,*) " please consider the use of epUI interface."
!      endif
!      nptcls = min(Nproduced, EPMAX_STACK)
!
!      do i = 1,  nptcls
!         ptclSave(i) = Pwork(i)
!      enddo
!      end
!      subroutine epqFirstColPtcls(ptcls, n, m)
!        implicit none
!#include  "Zptcl.h"
!#include "ZepMaxdef.h"
!       type(ptcl)::  ptclSave(EPMAX_STACK)
!      integer nptcls
!      common /stackSavec/ ptclSave, nptcls
!      integer n, m
!       type(ptcl)::  ptcls(n)
!      save
!      if(n .lt. nptcls) then
!         write(0,*)
!     *    " n must be > ", nptcls, " in epqFirstColPtcls"
!         stop  9999
!      endif
!      do m = 1, nptcls
!         ptcls(m) = ptclSave(m)
!      enddo
!      m = nptcls
!      end
