!   make next deinfed if first interaction iformation of nuclius break up
!   information is needed (make firstcollison=0 in epiev in ephook.f )
!
#define BREAKUPINFO
#undef BREAKUPINFO

#undef  DUMMYDATE
#if defined (NEXT486) || defined (KEKB)
#define DUMMYDATE
#elif defined PCLinux
#define DUMMYDATE
#elif defined CF_AlphaLinux
#define DUMMYDATE
#elif defined MACOSX
#define DUMMYDATE
#elif defined PCLinuxIFC
#define DUMMYDATE
#endif
#if defined DUMMYDATE
      subroutine IDATE()
      end
      subroutine ITIME()
      end
#endif
#undef DUMMYDATE

      subroutine cdpmjet(pj, ia, iz,  a, ntp)
      use dpmEmergency
      implicit none

#include "Zptcl.h"
#include "Zkfcode.h"
#include "Zmanagerp.h"
#include "Zevhnp.h"
#include "Zevhnv.h"
#include "Zcode.h"
#include "Zmass.h"
      
!    &&&&&&&&&&&&&&
!

      integer irsave(2)
      logical dodpm
      real*8    ELAB
      integer   IDP, NTMASS, NTCHAR, NPMASS, NPCHAR
      common /dpmjetcom/ ELAB, IDP, NTMASS, NTCHAR, NPMASS,NPCHAR
     *  ,irsave, dodpm
!    &&&&&&&&&&&
!
      type(ptcl):: pj  ! input. projectile particle
      integer ia     ! input. nucleon number of the target
      integer iz     ! input. charge no. of target
      type(ptcl):: a(*)  !  output. produced ptcls.  in cms.
      integer ntp   ! number of produced ptcls
!

      integer icon,  loopc, econloopc
      
      integer KKMAT
!///////////
      integer i
!/////////

      integer  NEVENT, ICASCA
      COMMON /DTEVNO/ NEVENT,ICASCA
      real*8 dErel, Ein
      integer econ
      real*8 sumTE, sumKE, sumRareTE, sumRareKE
      real*8 cdiff, diffmin
      real*8 cdiff0/0.03/

      integer  krarec
      logical first
      data first /.true./
      save 

#if defined (BREAKUPINFO)
!            use this when you want to get heavy ion breakup info.
!            at the first collision. (also put it in ephook --> epiev
!            to make firstcollision=0
      integer firstcollision
      common /checkglaub/  firstcollision
#endif
!///////////////
!////////////
!      logical deb
!      common /cdebug/ deb
!////////////

      if(first) then
         NEVENT = 0
!         ELAB = 3.0d0   
!         KKMAT = -1
!         call DT_KKINC(
!     *   1, 1, ia, iz,  1, ELAB, KKMAT, icon)
      endif

      KKMAT = -1
!      KKMAT = 1
      NEVENT = NEVENT + 1 

      cdiff = cdiff0  ! 3% relative error in Econsv. ok
      diffmin = 1.d20
      NTMASS = ia            
      NTCHAR = iz
      ELAB = pj%fm%p(4)
      if(pj%code .eq. kgnuc) then
         ELAB = ELAB/pj%subcode
      endif
      call ccos2dpjidp( pj, IDP, NPMASS, NPCHAR )
      dodpm = .false.

      if( IDP .eq. 0) then
!            dpm cannot treat this         
         if(pj%fm%p(4) .gt. 9.) then
            if(pj%code .eq. kkaon .and. LundPara(5) .ne. 0) then  ! not come now
               call chAcolAdhoc(pj, ia, iz, a, ntp)
            elseif(pj%code .eq. keta) then  ! not come now
               call chAcolAdhoc(pj, ia, iz, a, ntp)
            else
               call chAcolAdhoc(pj, ia, iz, a, ntp)    !  gzai,sigma 
!cc               call chANewLund(pj, ia, iz, a, ntp)
            endif
         else
            ActiveMdl ='byenergy'
!                 use fritiof 1.6 or nucrin
            call chALund(pj, ia, iz, a, ntp)
         endif
      else

         if(pj%code .eq. knuc .and.
     *      pj%subcode .ne. antip ) then
            dodpm = pj%fm%p(4) .gt. 1.6
         elseif(pj%code .eq. knuc .and. 
     *      pj%subcode .eq. antip ) then
!               for light element, 0.97 is o.k for PB
!               safe to use a higher value.
            dodpm = pj%fm%p(4) .ge. 0.99
         else
            if(pj%code .ne. kgnuc) then
!         
               if(pj%code .ne. kkaon .or. pj%charge .ne. 0) then
                  dodpm = pj%fm%p(4)  .gt. pj%mass + 0.6
!                      next 8 lies added in v.7.37 Oct25/2008
                  if(pj%code .eq. kpion .and. pj%charge .eq. -1
     *               .and. ia .eq. 1 ) then
                     dodpm = pj%fm%p(4)  .gt. 1.03
                  elseif( pj%code .eq. klambda .and. 
     *              pj%subcode .eq. antip 
     *              .and.  ELAB .lt. 4.0 .and. ia .eq. 1 ) then
                     dodpm = .false.
                  endif
               elseif(pj%subcode .ne. antip) then
                  dodpm  = pj%fm%p(4) .gt. 1.5
               endif
            else
!                 dpmjet cannot treat heavy int. below 5 GeV
              dodpm = ELAB .gt. 5.1
           endif
         endif
!           endif here       is moved to the last
         if(dodpm) then
            icon = 1
            loopc = 0 
            econ =1
            econloopc = 0 
            do while( econ .ne. 0)
!                 dpm; GeV/n
               do while (icon .ne. 0)
!&&&&&
!            call rnd1s(irsave)
! &&&&&
                  call DT_KKINC(
     *                 NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,ELAB,KKMAT,icon)
!               !!!!!!!! 2013.Apr.19
                  if( dpmEmergencyFlag == 1 ) then
                   ! should happen if proj = pi- and target is p
                     if( pj%code /= kpion .or. pj%charge /= -1) then
                        ! stragne
                        write(0, *)
     *                   ' ****warning:',
     *                   ' not pi-p but emergency in dpmjet3'
                        icon = 1
                     else
                         ! make charge +
!                        pj%charge = 1
!                        call ccos2dpjidp( pj, IDP, NPMASS, NPCHAR )
!                         NPCHAR is not pion charge 
!                         IDP=14 is pi- and IDP=13 is pi+-
                        IDP = 13
                        icon = 1
                     endif
                     dpmEmergencyFlag = 0
                  elseif(icon .ne. 0) then
!                   once the event is rejected, even if different
!                   random number cannot get rid of this loop. 
!                   so we change the energy by 1 %
                     ELAB=ELAB*1.01   
                     loopc = loopc + 1
!                  if(loopc .gt. 100) then    ! probably useless.
                     if(loopc .gt.  20) then
                        write(ErrorOut,*)
     *                       ' pj code, subcode, charge=',pj%code,
     *                       pj%subcode, pj%charge, ' E=',pj%fm%p(4),
     *                       ' Target A=', ia, ' Z=', iz,
     *                       ' loopc=',loopc
                        call cpdpmjetinp
                        call cerrorMsg(
     *                        'dpmjet3 too many loop in cdpmjet',1)     
!                            no ptcl produced
                        call cdpmNoPtcl
                        icon = 0
                        econ = 0
                     endif
                  endif
                  if(first) then
                      call cdpmRareDecay( dpmRareDecay ) ! control
                      if( dpmRareDecay == 1 .or.
     *                     dpmRareDecay == 2 ) then
                         icon = 1 ! skip first event which may have 
                                ! not decayed  rare particles
                      endif
                      first = .false. 
                   endif
                enddo

                if(econ .ne. 0) then  ! first time econ /= 0
                   call csetdpmresul(pj, ia, a, ntp,
     *               sumTE, sumKE, sumRareTE, sumRareKE, krarec,
     *               dpmRareDecay)

                  Ein = NTMASS*masp + pj%fm%p(4)
                  dErel = sumTE/Ein -1.
!////////////////
!                  if(firstcollision==0)then
!                     write(0,*) 
!     *                    ' dErel=',dErel, ' sumTE=',sumTE, ' Ein=',Ein
!                  endif
!////////////
                  if( abs(dErel) .lt. abs(diffmin)) then
                     diffmin = dErel
                  endif

                  if(abs(dErel) .gt. cdiff) then
!  
                     econloopc = econloopc +1
                     if(econloopc .gt. 10 
     *                  .and. cdiff .ne. cdiff0 ) then
                        write(0,*) 'EinT=',Ein, ' EoutT=',sumTE
                        write(0, *)'EinK=', pj%fm%p(4)-pj%mass
                        write(0,*) 'EoutK=', sumKE
                        write(0,*) 'rare p=',krarec
                        write(0,*) 'RareTE=', sumRareTE,' KE=',sumRareKE
                        write(ErrorOut,*) ' dpmjet3 Econsv difficult ',
     *                       ' proj code, subcode, charge=',pj%code,
     *                       pj%subcode, pj%charge, ' E=',pj%fm%p(4),
     *                       ' Target A=', ia, ' Z=', iz,
     *                       ' min Rel dE=',  dErel,
     *                       ' loopc=',econloopc
                        call cpdpmjetinp
                     endif
                     if(cdiff .ne. cdiff0) then
                        econ =0
                     else
                        cdiff = abs(diffmin)*1.2
                        econloopc = -5
                     endif
                  else
                     econ =0
                  endif
               endif
            enddo
         else
            if(pj%code .le. knuc) then
!               use fritiof 1.6 or nucrin
               ActiveMdl = 'byenergy'
               call chALund(pj, ia, iz, a, ntp)
            elseif(pj%code .eq. kgnuc) then
!  
               call cerrorMsg('Recursive call of heavy routine',0)
!                    dpmjet is bypaseed for heavy < 5GeV
!               if(ELAB .gt. 0.5) then
!                  ActiveMdl = 'byenergy'
!                  call cheavyInt(pj, ia, iz, a, ntp)
!               else
!                  a(1) = pj
!                  ntp = 1
!               endif
            else
!                  neglect this interaction
               ntp = 0
            endif
         endif
      endif
!/////////////////////
      do i = 1, ntp
         if(a(i)%code .eq. kphoton .and.  a(i)%charge .ne.  0) then
            write(0,*)
     *       'dpmerr  photon charge=' ,a(i)%charge, ' subc=',
     *       a(i)%subcode, ' E=',a(i)%fm%p(4)
            write(0,*)
     *      'dpmerr dodpm=',dodpm,' ntp=', ntp, ' pj=', pj%code,
     *       ' pj%charge=', pj%charge, ' ia=',ia, ' iz=', iz
         endif
      enddo
!/////////////////

      end
      subroutine cpdpmjetinp
      implicit none
!            print current dpmjet input
      character*100 msg
!    &&&&&&&&&&&&&&
      integer irsave(2)
      logical dodpm
      real*8    ELAB
      integer   IDP, NTMASS, NTCHAR, NPMASS, NPCHAR
      common /dpmjetcom/ ELAB, IDP, NTMASS, NTCHAR, NPMASS,NPCHAR
     *   ,irsave, dodpm
!    &&&&&&&&&&&
      character*16 model
      call cqActiveMdl(model)
      write(msg, *) ' dodpm=',dodpm, ' activemodel=', model
      call cerrorMsg(msg, 1)
      write(msg, *) ' dpmjet input: IDP=', IDP, ' Elab=',ELAB
      call cerrorMsg(msg, 1)
      write(msg, *) ' target: A=', NTMASS, ' Z=',NTCHAR,
     *   ' NPA=',NPMASS, ' NPZ=',NPCHAR
      call cerrorMsg(msg, 1)
!
      write(msg, *) ' ir=', irsave
      call cerrorMsg(msg, 1)
      end

      subroutine csetdpmresul(pj, tgA, a, ntp, sumTE, sumKE, sumRareTE,
     *      sumRareKE, krarec, dpmRareDecay)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "Zcode.h"
#include "Zptcl.h"
      type(ptcl):: pj    ! input. projectile 
      integer,intent(in)::tgA ! input. target A  added Nov.9,2013
      type(ptcl):: a(*)  ! output.  
      integer ntp  ! ouptut. number of ptcls put in a
      real*8  sumTE  ! output. produced ptcls' total energy 
      real*8  sumKE  ! //                      kinetic
      real*8  sumRareTE  ! //       rare ptcl's total E
      real*8  sumRareKE  ! //       rare ptcl's kinetic E
      integer krarec ! output rare particle counter
      integer,intent(in):: dpmRareDecay ! how to control rare 
                                      ! ptcls
! COMMON /DTEVT1/ :
!                   NHKK         number of entries in common block
!                   NEVHKK       number of the event
!                   ISTHKK(i)    status code for entry i
!                   IDHKK(i)     identifier for the entry
!                                (for particles: identifier according
!                                 to the PDG numbering scheme)
!                   JMOHKK(1,i)  pointer to the entry of the first mother
!                                of entry i
!                   JMOHKK(2,i)  pointer to the entry of the second mother
!                                of entry i
!                   JDAHKK(1,i)  pointer to the entry of the first daughter
!                                of entry i
!                   JDAHKK(2,i)  pointer to the entry of the second daughter
!                                of entry i
!                   PHKK(1..3,i) 3-momentum
!                   PHKK(4,i)    energy
!                   PHKK(5,i)    mass
!

!  The final state particles from the actual event (number NEVHKK)
!  can be found in DTEVT1 and identified by their status:
!
!     ISTHKK(i) = 1    final state particle produced in
!                      photon-/hadron-/nucleon-nucleon collisions or
!                      in intranuclear cascade processes
!                -1    nucleons, deuterons, H-3, He-3, He-4 evaporated
!                      from excited nucleus and
!                      photons produced in nuclear deexcitation processes
!                1001  residual nucleus (ground state)
!
!  The types of these particles/nuclei are given in IDHKK as follows
!
!     all final state part. except nuclei :
!       IDHKK(i)=particle identifier according to PDG numbering scheme
!     nuclei (evaporation products, and residual nucleus) :
!       IDHKK(i)=80000, IDRES(i)=mass number, IDXRES(i)=charge number
!
!  The 4-momenta and masses can be found in PHKK (target nucleus rest frame):
!                   PHKK(1..3,i) 3-momentum (p_x,p_y,p_z)
!                   PHKK(4,i)    energy
!                   PHKK(5,i)    mass
!
!
      PARAMETER (NMXHKK=90000)

      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

! extended event history
      COMMON /DTEVT2/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),
     &                IHIST(2,NMXHKK)

      integer i
      integer code, subcode, charge
#if defined (BREAKUPINFO)
      integer firstcollision
      common /checkglaub/  firstcollision
      character(2):: id
!            counters
      integer ngNc   !  generated nucleons
      integer neNPc  !  evaporated nucleons from projectile 
      integer neNTc  !  evaporated nucleons from target
      integer neHPc   ! evaporated heavies from projectile 
      integer neHTc   ! evaporated heavyies from target 
      integer nrPc   ! remnant  heavies from projectile
      integer nrTc   ! remnant  heavies from target 
      real(8):: pjKEn  ! projecitle kinetic energy /n
      real(8):: thresh   ! pjKEn*0.5. if KE/n of evaporation > this
                   ! it should be from projectile 
#endif
!      if(NHKK .gt. 1000 ) then
!         write(0,*) ' NHKK=',NHKK
!      endif
!///////////////
      ntp =0
      sumKE = 0.
      sumTE = 0.
      krarec = 0
      sumRareTE =0.
      sumRareKE = 0.
#if defined (BREAKUPINFO)
      ngNc = 0
      neNPc = 0
      neNTc = 0 !  evaporated nucleons from target
      neHPc = 0  ! evaporated heavies from projectile 
      neHTc  = 0  ! evaporated heavyies from target 
      nrPc  = 0  ! remnant  heavies from projectile
      nrTc  = 0 ! remnant  heavies from target 
      pjKEn = pj%fm%p(4) + pj%mass
      if(pj%code == kgnuc ) then
         pjKEn = pjKEn/pj%subcode
      endif
      thresh = pjKEn*0.5
#endif
      do i=1, NHKK
         IF (ABS(ISTHKK(i)).EQ.1) THEN
            if(IDHKK(i) .lt. 79000) then
               call ckf2cos(IDHKK(i), code, subcode, charge)
!                 neglect all rare patlces such as D0,D+...
               if( dpmRareDecay == 0 ) then
                  call cdpmRareDecay2(IDHKK(i), code)
               endif
#if defined (BREAKUPINFO)
               if( firstcollision == 0 ) then
                  if( code == knuc ) then
                     if(ISTHKK(i) > 0 ) then
                        id = 'G'
                        ngNc = ngNc + 1
                     else
                        if(( PHKK(4,i) - PHKK(5,i)) >   thresh) then
                           neNPc = neNPc + 1
                           id = 'PE'
                        else
                           id = 'TE'
                           neNTc = neNTc + 1
                        endif
                     endif
                     write(*,'(a,3i4, a, 1p, g14.4,i3)')
     *                  id, code, subcode, charge,
     *                 ' KE/n =', PHKK(4,i) - PHKK(5,i),ISTHKK(i)
                  endif
               endif
#endif
            else
               charge = IDXRES(i)
               code = kgnuc
               subcode = IDRES(i)
#if defined (BREAKUPINFO)
               if(firstcollision == 0 ) then
                  if( (PHKK(4,i) - PHKK(5,i))/subcode > thresh) then
                     neHPc =  neHPc + 1
                     id = "PE"
                  else
                     neHTc =  neHTc + 1
                     id = "TE"
                  endif
                  write(*,'(a,3i4, a, 1p, g14.4,i6)')
     *            id, code, subcode, charge,
     *            ' KE/n =', (PHKK(4,i) - PHKK(5,i))/subcode, ISTHKK(i)
               endif
#endif
            endif
            sumKE = sumKE + (PHKK(4,i) - PHKK(5,i))
            sumTE = sumTE + PHKK(4,i)
            if(code .ne. krare) then
               ntp = ntp +1
               call cmkptc(code, subcode, charge, a(ntp))
               a(ntp)%fm%p(1) = PHKK(1,i)
               a(ntp)%fm%p(2) = PHKK(2,i)
               a(ntp)%fm%p(3) = PHKK(3,i)
               a(ntp)%fm%p(4) = PHKK(4,i)
               if( code == kgnuc)  then
                  a(ntp)%mass =  PHKK(5,i)
               endif
            else
               krarec = krarec + 1 
               sumRareKE = sumRareKE +  (PHKK(4,i) - PHKK(5,i))
               sumRareTE = sumRareTE +  PHKK(4,i) 
            endif
!            write(*,*)
!     *      IDHKK(i), code, subcode, charge,
!     *      sngl(PHKK(4,i)),
!     *      sngl(PHKK(1,i)),sngl(PHKK(2,i)),sngl(PHKK(3,i))
         ELSEIF (ABS(ISTHKK(i)).EQ.1001) THEN
!             residual nucleus;
            charge = IDXRES(i)
            code = kgnuc
            subcode = IDRES(i)
            ntp = ntp +1
            call cmkptc(code, subcode, charge, a(ntp))
            a(ntp)%mass =  PHKK(5,i)
            a(ntp)%fm%p(1) = PHKK(1,i)
            a(ntp)%fm%p(2) = PHKK(2,i)
            a(ntp)%fm%p(3) = PHKK(3,i)
            a(ntp)%fm%p(4) = PHKK(4,i)
            sumKE = sumKE + (PHKK(4,i) - PHKK(5,i))
            sumTE = sumTE + PHKK(4,i)
#if defined (BREAKUPINFO)
            if( firstcollision == 0 ) then
               if((PHKK(4,i) - PHKK(5,i))/subcode > thresh) then
                  nrPc = nrPc + 1
                  id = "PR"
               else
                  nrTc = nrTc + 1
                  id = "TR"
               endif
               write(*,'(a,3i4, a, 1p, g14.4,i6)')
     *          id, code, subcode,charge,
     *          ' KE/n=', (PHKK(4,i) - PHKK(5,i))/subcode, ISTHKK(i)
            endif
!            write(0,*) ' p=',PHKK(1,i), PHKK(2,i), PHKK(3,i),PHKK(4,i)
#endif
         ENDIF
      enddo
!        aritificial symmetry recovering   
      call cASRforDPM(pj, tgA, a, ntp)

      call crot3mom(pj, a, ntp)
      
#if defined (BREAKUPINFO)
      firstcollision = 1
      id = "No"
      write(*,'(a)') '   ntp  ngNc  neNPc neNTc neHPc neHTc  nrPc  nrTc'
      write(*, '(a,8i6)') id,
     *  ntp, ngNc, neNPc, neNTc, neHPc, neHTc, nrPc, nrTc
#endif
      end
      subroutine cASRforDPM(pj, tgA, a, ntp)
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
#include "Ztrackp.h"
      type(ptcl):: pj  ! input.  projectile particle
      integer,intent(in):: tgA ! input Target A
      integer,intent(in):: ntp  ! # of produced particles
                                ! in a
      type(ptcl):: a(ntp)  ! in/out. 


      
      integer,parameter::minA=5 ! which may have asymm
      real(8),parameter::maxPz=2.0 ! GeV/c
      real(8):: cosf, sinf, pt, u
      integer:: i, j
      type(ptcl):: proton, b 
      
      if( ASRforDPM == "no") return !!!!!!
!       both target and progjectile is non heavy (< minA)
!     nothing to do
!     if( ( pj.code /= kgnuc  .or. pj.subcode < minA )
!           *     .and. tgA < minA )  return !!!!!!!!!!
!       changed to (2016. Jun. )
      if ( tgA < minA )  return !!!!!!!!!!    
      if( ASRforDPM == "m" ) then
!          reverse the sign of Px with 1/2 prob.
         call rndc(u)
         if(u < 0.5) then 
            a(:)%fm%p(1) = - a(:)%fm%p(1)
         endif
      elseif( ASRforDPM(1:1) == "r" )  then
           ! now r or r1. This business may better done in CMS
           ! but we use pj or tg rest frame for simplicity.
         if(  tgA >= minA ) then
            ! target backscater region both for r and r1
            do i = 1, ntp
               if( a(i)%code == kpion .or.
     *              a(i)%code == kkaon  ) then
                  if(  a(i)%fm%p(3) < maxPz ) then
                     pt = sqrt( a(i)%fm%p(1)**2 +a(i)%fm%p(2)**2 )
                     call kcossn(cosf, sinf)
                     a(i)%fm%p(1)= pt*cosf
                     a(i)%fm%p(2)= pt*sinf
                  endif
               endif
            enddo
         endif
         if( ASRforDPM /= "r1" .and.
     *       pj%code == kgnuc .and. pj%subcode >= minA) then
!                proj. very forward region. see it at pj%rest
!                system.
            j = 0
            do i = 1, ntp
               j = j + 1
               if( j == 1 ) then
                  call cmkptc(knuc, -1, 1, proton)
                  proton%fm%p(1:3) = pj%fm%p(1:3)/pj%subcode
                  proton%fm%p(4) =
     *             sqrt( sum(proton%fm%p(1:3)**2+proton%mass**2))
               endif         
                    ! boost to pj /n  rest frame
               call cbst1(j, proton, a(i), b)
               if(  b%fm%p(3) > - maxPz ) then
                  pt = sqrt( a(i)%fm%p(1)**2 +a(i)%fm%p(2)**2 )
                  call kcossn(cosf, sinf)
                  a(i)%fm%p(1)= pt*cosf
                  a(i)%fm%p(2)= pt*sinf
               endif
            enddo
         endif
      else
         write(0,*) ' ASRforDPM=', ASRforDPM, ' invalid'
         stop
      endif
      end  subroutine cASRforDPM

      subroutine cdpmNoPtcl
!         clear no. of produced partilce
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!
      PARAMETER (NMXHKK=90000)

      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

      NHKK= 0
      end


      SUBROUTINE DT_USRHIS(MODE)
      implicit none
      integer MODE
      end


!      IDP number
!     'PROTON  ' , 'APROTON ' , 'ELECTRON' , 'POSITRON' ,
!        1            2            3             4  
!     'NEUTRIE ' , 'ANEUTRIE' , 'PHOTON  ' , 'NEUTRON ' , 'ANEUTRON' ,
!        5             6           7             8             9
!     &'MUON+   ' , 'MUON-   ' , 'KAONLONG' , 'PION+   ' , 'PION-   ' ,
!       10             11          12           13           14
!     &'KAON+   ' , 'KAON-   ' , 'LAMBDA  ' , 'ALAMBDA ' , 'KAONSHRT' ,
!       15            16           17           18           19
!     &'SIGMA-  ' , 'SIGMA+  ' , 'SIGMAZER' , 'PIZERO  ' , 'KAONZERO' ,
!       20            21           22           23           24  
!     &'AKAONZER' , 'NEUTRIM ' , 'ANEUTRIM' , 'NEUTRIT ' , 'ANEUTRIT' ,
!       25            26           27           28           29 
!     
      subroutine ccos2dpjidp(pj, idp, npmass, npchar)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
!          cosmos code to idp in dpmjet
      type(ptcl)::  pj
      integer idp  ! output. dpm code. which can be a projectile
      integer npmass ! ouptut.  projectile code (mass) if pj is
                     ! used as a projectile
      integer npchar ! output. charge in the same meaning as above.

      if(pj%code .eq. kgnuc) then
         idp = 1
         npmass = pj%subcode
         npchar = pj%charge
      else
         npmass = 1
         npchar = 1
         if(pj%code .eq. knuc) then
            if(pj%charge .eq. 1) then
               idp=1
            elseif(pj%charge .eq. 0) then
               if(pj%subcode .eq. antip) then
                  idp = 9
               else
                  idp = 8
               endif
            else
               idp = 2
            endif
         elseif( pj%code .eq. kpion) then
            if(pj%charge .eq. 1) then
               idp = 13
            elseif(pj%charge .eq. -1) then
               idp = 14
            else
               idp = 23
            endif
         elseif( pj%code .eq. kkaon) then
            if(pj%charge .eq. 1) then
               idp = 15
            elseif(pj%charge .eq. -1) then
               idp = 16
            elseif(pj%subcode .lt. 0) then
!                 k0; disregard  short or long but see ptlc or anti
               idp = 24
            else
!                 k0_bar
!ccc               idp = 25
               idp = 24   ! tentative setting. k0bar will die easily
            endif
         elseif(pj%code .eq. klambda) then
            if(pj%subcode .eq. antip) then
               idp = 18
            else
               idp = 17
            endif
         elseif(pj%code .eq. kphoton) then
            idp = 7
         elseif(pj%code .eq. ksigma) then
            if(pj%charge .eq. 1) then
               idp = 21
            elseif(pj%charge .eq. -1) then
               idp = 20
            else
               idp = 22
            endif
         else
!              cannot be a projectile
            idp = 0
         endif
      endif
      end
      subroutine cinidpmjet(file)
      implicit none
#include "Zmanagerp.h"
#include "Zevhnp.h"

      character*(*)  file  ! input. first input file to initialize
                           ! dpmjet3. 
      
      integer icon
      integer NEVTS,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IEMU
      real*8  EPN
!!!!!!!!!!
      write(0,*) ' cinidpmjet: top file=', trim(file)
!!!!!!!!!!
      call copenf(TempDev, file, icon)
      if(icon .ne. 0) then
         call cerrorMsg( file, 1)
         call cerrorMsg('above file cannot be opened', 0)
      endif

!          init dpm; Thru this call, input data is read via TempDev
!     Glauber initialization data location is given in that
!     data and the call to cdpmOpen is made from within dpmjet
      IEMU = 0
      CALL DT_DTUINI(NEVTS,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IEMU)

      close(TempDev)
      end
      subroutine cdpmOpen(io, file)
      implicit none
#include "Zmanagerp.h"      
!      since dpmjet cannot use long file name, we may encounter
!      difficulty in dealing with full path which is needed to
!      specify the files when distributing jobs over a 
!      number of  workstations.  So I had to make a
!      correction in the dpmjet program;Say, the open statements
!      for reading GLAUB-PAR data

!         OPEN(LLOOK,FILE=CGLB//'.GLB',STATUS='UNKNOWN')
!         OPEN(LLOOK,FILE=CGLB(1:I-1)//'.GLB',STATUS='UNKNOWN')
!
!      are replaced  by 
!      
!           call cdpmOpen(LLOOK, CFILE)
      
      integer io  ! input.  file descripter to be opened
      character*12 file  ! input. file name, say, containig glauber
                   !  calculation result, such as atmos.GLB. 
                   ! the full path is assmued to be TopDir/file


      character*100 path


      integer icon

!!!!!!!!!!
      write(0,*) ' cdpmOpen: top io, file=',io, trim(file)
!!!!!!!!!!
      if( DpmFile == " " ) then
         call cformFullPath(file, path)
      else
         path = file
      endif
!!!!!!!!!!
      write(0,*) ' after fullapth: path=', trim(path)
!!!!!!!!!!
      call copenf(io, path, icon)
      if(icon .ne. 0) then
         call cerrorMsg(path, 1)
         call cerrorMsg('above file cannot be opened', 0)
      endif
      end


      subroutine cdpmOpen2(io, file)
      implicit none
#include "Zmanager.h"
!      This opens a file with Cosmos/Data/DPM as a prefix
!    
      integer io  ! input.  file descripter to be opened
      character*(*) file  ! input. file name.

      character*100 path
      integer icon, klena
!!!!!!!!!!
      write(0,*) ' cdpmOpen2: top io, file=', io, trim(file)
!!!!!!!!!!

      path = ' '
      path = TopDir(1:TopDirLeng)//'/Data/DPM/'
     *       //file(1:klena(file))

      
      call copenf(io, path, icon)
!!!!!!!!!!
      write(0,*) ' cdpmOpen2: path=',trim(path)
!!!!!!!!!!

      if(icon .ne. 0) then
         call cerrorMsg(path, 1)
         call cerrorMsg('above file cannot be opened', 0)
      endif
      end
      module moddpmRareDecay
      integer,parameter::nrare=19   ! # of rare ptcls
      integer,parameter::nrare2=11 ! first 11 ptcls
                     ! are treatable in Cosmos decay routine
                ! and collision is made by as proton or Kaon
      integer,save::todecay(nrare)=      ! there kf code
     *  (/ 411, -411,  421, 3222, -3222   3212,  3322, 
    !      D+    D-     D0  Sigma+ Sigma- Sigma0, Xi0,
     *  3312,  -3312,  3334,   4122, 
    !    Xi+    Xi-    Omega- LamdaC+

     *   431, -431,  441,  443,  4112,   
    !     Ds+  Ds-   etaC J/psi SigmaC0,
     *  4132,   4212,
    !    Xi0C  SigmaC+
     *    4222,    4232/)
    !    SigmaC++   XiC+ 
      end module moddpmRareDecay

      subroutine cdpmRareDecay( dpmRareDecay )
!       To make very rare short life particles decay,
!       LUND-MDCY input  in dpmjet.inp may be given, but
!       it seems reset inside  dpmjet (DT_INITJS).
!       So we reset it again. should be called after
!       dpmjet init. once for all.  
      use moddpmRareDecay
      implicit none
      integer,intent(in):: dpmRareDecay

      integer MDCY, MDME, KFDP
      real(8):: BRAT
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)

      integer:: i
      integer:: kc, PYCOMP

      if( dpmRareDecay == 1 ) then
         do i = nrare2+1, nrare
            kc =PYCOMP(todecay(i))
            MDCY(kc, 1) = 1     ! force to decay (0--> non decay)
         enddo  
      elseif( dpmRareDecay == 2 ) then
         do i =1, nrare
            kc =PYCOMP(todecay(i))
            MDCY(kc, 1) = 1     ! force to decay (0--> non decay)
         enddo  
      elseif( dpmRareDecay == 3 ) then
         ! nothing to do; completely as old versions 
      elseif( dpmRareDecay == 0 ) then    
         ! here nothing to do; neglect all rare ptcls; later 
         ! do so.
      else
         write(0,*) 'call cdpmRareDecay( dpmRareDecay )'
         write(0,*) ' dpmRareDecay =',dpmRareDecay, ' invalid'
         stop
      endif
      end
      subroutine cdpmRareDecay2(kf, code)
      use moddpmRareDecay
      implicit none
#include "Zcode.h"
      ! to be called only  when dpmRareDecay = 0 and
      ! when converting dpmjet code into cosmos code
      ! This routine is to neglect all rare particles
      ! including D+,D0  .. as defined i=1, nrare2
      ! in todecay(i).  To do so, code is given krare.
      integer,intent(in):: kf  ! dpmjet ptcl code (kf code)
      integer,intent(inout):: code ! if if code is /= krare and 
                          !  one of todecay(i), krare is given
                          ! else unchanged.
      integer::i
      
      if( code /= krare ) then
         do i = 1, nrare2
            if( kf == todecay(i) ) then
               code = krare
               return
            endif
         enddo
      endif
      end
      subroutine cqDpmIntType(inttype)
      implicit none
!
! general process information; only for dpmjet3
      INTEGER IPROCE,IDNODF,IDIFR1,IDIFR2,IDDPOM,IPRON
      COMMON /POPRCS/ IPROCE,IDNODF,IDIFR1,IDIFR2,IDDPOM,IPRON(15,4)
!      IPROCE
!             1 non-diffractive inelastic
!             2 elestic 
!             3 quasi elestic vector meson prod. (photon)
!             4 central diffraction
!             5 single diff. ptcl 1
!             6 //           ptcl 2
!             7 double diff. 
!             8 direct photo-hadron
! For moore detail, see manual in Documents/CosmicRays/phojetShort.pdf
!               say, IDIFR1 classifies IPROCE=5
      
      integer,intent(out):: inttype

      inttype = IPROCE
      end
      
