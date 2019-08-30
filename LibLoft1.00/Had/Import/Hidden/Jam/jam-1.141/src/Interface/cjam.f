!        if JamFragment=1 and if need to debug fragment infomation make
!       next defined
#undef ShowFragment

      subroutine cjamEvent(pj, ia, iz, sig, a, ntp)
      implicit none
#include "Zcode.h"      
#include "Zptcl.h"      
#include "Ztrackp.h"
      type(ptcl):: pj  !  input. projectile

      integer,intent(in)::ia, iz ! target A and Z
      real(8),intent(in):: sig  ! cross seccion in mb on this target
      type(ptcl):: a(*)  !output. generated ptcls
      integer,intent(out):: ntp  ! # of generated ptcls

      integer::ntp2
      logical ok
      integer inela, loop1, loop2

      integer,parameter::maxloop=10
 
      type(ptcl):: tg  !
      type(ptcl):: tgspec(209), pjspec(209)
      integer:: ntgspec, npjspec, nout
      real(8):: impactp 

!/////////
      integer:: i
!//////////
      if(ia == 1 ) then
         call cmkptc(knuc, -1, iz, tg)
      elseif(ia > 1 ) then
         call cmkptc(kgnuc, ia, iz, tg)
      else
         write(0,*) 'in cjamEvent: strange target A,Z',ia,iz
         stop
      endif
!        if anti-p has  0 K.E, 
!        only elastic events are obtained; avoid loop
      if(  pj%code== knuc .and. pj%charge==-1
     *     .and.  pj%fm%p(4) - pj%mass <=0.0 ) then
         pj%fm%p(4) = 0.001d0 + pj%mass
!         write(0,*) ' stopping pbar '
!           force to p target; otherwise infinite loop.
         call cmkptc(knuc, -1, 1, tg)
      endif
      tg%fm%p(4) = tg%mass
      tg%fm%p(1:3)= 0.
      loop1 = 0
      loop2 = 0
      ok = .false.
      do while (.not. ok)
         call cjamini(pj, tg, sig)
         call cjam(pj, ia, iz,  a, ntp2, 
     *     pjspec, npjspec, tgspec, ntgspec, inela, impactp)

!          jam cannot disable elastic collision.
!        for non heavy, we can use total x-section
!        which includes elastic one, but for heavy
!        we cannot get elastic cross-section. So we have to 
!          wait until inelastic collsion for heavy.
!        
!//////////////
!         if( inela == 0 ) then
!            sumTE= 0.
!            write(0, *)
!            write(0,'(a,3i4, 1p,g12.3)')  'inela ',inela, pj.code, ia,
!     *                     pj.fm.p(4)-pj.mass
!            write(0,'(a,3i4,1p, g12.3, a,1p,g12.3)')
!     *        'pj=', pj.code, pj.subcode, pj.charge, pj.mass,
!     *        'KE=', pj.fm.p(4)- pj.mass
!            write(0,*) '# of  pj spec=',npjspec
!            write(0,*) 'target A,Z =',ia, iz, tg.mass
!            write(0,*) '# of  tg spec=',ntgspec
!
!            do i = 1, ntp2
!               write(0,'(a,3i5,1p,g12.4)') 'in',
!     *              a(i).code, a(i).subcode, a(i).charge, 
!     *              a(i).fm.p(4) - a(i).mass
!               sumTE = sumTE + a(i).fm.p(4)
!            enddo
!            do i = 1, npjspec
!               write(0,'(a,3i5,1p,g12.4)') 'pjsp',
!     *           pjspec(i).code, pjspec(i).subcode, pjspec(i).charge, 
!     *           pjspec(i).fm.p(4) - pjspec(i).mass
!               sumTE = sumTE + pjspec(i).fm.p(4)
!            enddo
!            do i = 1, ntgspec
!               write(0,'(a,3i5,1p,g12.4)') 'tgsp',
!     *           tgspec(i).code, tgspec(i).subcode, tgspec(i).charge, 
!     *           tgspec(i).fm.p(4) - tgspec(i).mass
!               sumTE = sumTE + tgspec(i).fm.p(4)
!            enddo
!            write(0,*) ' sumTE=',sumTE
!            write(0,*) ' pj E + tgmass=', pj.fm.p(4) + tg.mass
!         endif
!/////////////// 
         if( JamXs == 1 ) then
            ok = pj%code /= kgnuc .or.  inela /= 0
         elseif(JamXs == 0 ) then
            ok = inela /= 0
         else
            write(0,*) ' JamXs=',JamXs, ' not usable in cjam.f'
            stop
         endif
         if(.not. ok ) then
            loop1 = loop1 + 1
            if(loop1 > maxloop) then
               write(0,*)
     *           'In jam event gen. JamXs=',JamXs, ' inela=',inela
               write(0,*) 'pj=',pj%code, pj%subcode,pj%charge
               write(0,*) 'pj KE=',pj%fm%p(4)- pj%mass
               write(0,*) 'tg A,Z=', ia,iz
               write(0,*)' generated # of ptcls =',ntp2
               write(0,*) 'event is neglected'
               ntp = 0
               exit
            endif
         endif
!         # of npjesec may become 1 for pion/koan
!           for elastic event.  (JamXs=1-->ok= f
!                                JamXs=0-->ok= t. but 
!                           cjamPjSpec dose not generate any frgment
!                           so don't worry)
#if defined ShowFragment 
         write(0,'(a,3i4,a,1p,g12.3,a, l)')
     *        'pj=', pj%code, pj%subcode, pj%charge,
     *        'KE=', pj%fm%p(4)- pj%mass, ' ok=',ok
         write(0,*) '# of  pj spec=',npjspec
         write(0,*) 'target A,Z =',ia, iz
         write(0,*) '# of  tg spec=',ntgspec

         do i = 1, ntp2
            write(0,'(a,3i5,1p,g12.4)') '---',
     *          a(i)%code, a(i)%subcode, a(i)%charge, 
     *         a(i)%fm%p(4) - a(i)%mass
         enddo
#endif

         if(ok) then
            if(npjspec > 0 ) then
               call cjamPjSpec(pj, tg,
     *            pjspec, npjspec, impactp, inela, a(ntp2+1), nout)
            else
               nout = 0
            endif
!/////////
#if defined ShowFragment
            write(0,*) 'pj spec-->frag=',nout
            do i = 1, nout
               write(0,'(a,3i4,1p,g12.3)')
     *          '***p fr', a(ntp2+i)%code, 
     *         a(ntp2+i)%subcode,  a(ntp2+i)%charge, 
     *         a(ntp2+i)%fm%p(4)-a(ntp2+i)%mass
            enddo
#endif
!/////////// 
            ntp2 = ntp2 + nout
            if(ntgspec > 0 ) then
               call cjamTgSpec(pj, tg,
     *           tgspec, ntgspec, impactp, inela, a(ntp2+1), nout)
            else
               nout = 0
            endif
!/////////
#if defined ShowFragment
            write(0,*) 'tg spec-->frag=',nout
            do i = 1, nout
               write(0,'(a,3i4,1p,g12.3)')
     *          '^^^t fr', a(ntp2+i)%code, 
     *         a(ntp2+i)%subcode,  a(ntp2+i)%charge, 
     *         a(ntp2+i)%fm%p(4)-a(ntp2+i)%mass
            enddo
            write(0,*)
#endif
!//////////////
            ntp2 = ntp2 + nout
            call cjamEconsv(pj, tg, a, ntp2, ntp, ok)
            if(.not. ok ) then
               loop2 = loop2 + 1
               if(loop2 > maxloop) then
                  exit 
                 ! use current event with 10 % E error
               endif
            endif
         endif
      enddo
      call cjamElaInfo(0, inela)  ! inform inela
      call crot3mom( pj, a, ntp )   ! rotate to the current  cooord
      end
      subroutine cjamElaInfo(inout, inela)
      implicit none
      integer,intent(in)::   inout ! 0--> inela is input
                                   ! 1--> //    output
      integer,intent(inout):: inela  ! 0--> jam interaction is elastic
                                  !  1--> //                  inelastic
      integer,save::saveinela  ! keep input inela
      
      if(inout == 0) then
         saveinela = inela
      else
         inela = saveinela
      endif
      end
      subroutine cjamEconsv(pj, tg, a, ntp2, ntp, ok)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zkfcode.h"
      type(ptcl):: pj  ! input.  projectile
      type(ptcl):: tg  ! input.  target
      integer,intent(in):: ntp2   ! # of ptcls in a
      integer,intent(out):: ntp   ! # of ptcls in a (after dropping
                                 ! spectator n
      type(ptcl)::a(ntp2)   ! in/out. produced ptcles
      logical,intent(out):: ok  ! if E consv, ok=.true. else f
      
      integer i


!           total energy conseration;  if target is proton, very accurate
!           if target is A,  and projectile kinetic energy/n is small
!           binding  energy of A affects the apparent conseration
!           (if the mass of A is taken to be (massN*(A-Z) + massP*Z)
!           if KE/n goes higher, apparent conservation becomes very good.
!           at 2.5 GeV/n, relative error becomes about ~3 %. so we set
!           it eps=0.035; some rare events with error> eps  may be retried

      real(8):: sumTE, projA
      real(8),parameter:: eps=0.035d0

      sumTE = 0.
      ntp = 0
      do i = 1, ntp2
         sumTE = a(i)%fm%p(4)  + sumTE 
         ntp = ntp + 1
         a(ntp) =  a(i)
      enddo
      if( abs( sumTE/(pj%fm%p(4) + tg%mass) - 1.0d0) > eps ) then
         if( pj%code == kgnuc ) then
            projA = pj%subcode
         else
            projA= 1.
         endif
!         write(0,'(a, 4i4, 1p,2g12.4)')
!     *    'E consv err;retry', tg.code, tg.subcode, pj.code, pj.subcode,
!     *            (pj.fm.p(4)-pj.mass)/projA,
!     *    sumTE/(pj.fm.p(4) + tg.mass) - 1.0d0
         ok = .false.
      else
         ok = .true.
      endif
      end

      subroutine cjamini(pj, tg, sig)
!      This must be called before generating one exclusive 
!      event (call jamev) 
!      implicit none  ! this cannot be used since jam1.inc uses
!      implicit double (...)
#include "jam1.inc"
#include "jam2.inc"

#include "Zptcl.h"
#include "Zkfcode.h"
#include "Zmanagerp.h"
#include "Zevhnp.h"
#include "Zevhnv.h"
#include "Zcode.h"
#include "Zmass.h"
      type(ptcl)::pj, tg  ! input 
      real(8),intent(in)::sig  ! input  cross-sction on target tg in mb
                         ! used to calculate the impact parameter

      integer ia, iz    ! input target nucleus

      character frame*8, proj*8, targ*8, cwin*15
      real*8 bmin, bmax, dt
      integer mevent, nstep
      integer first /1/
      save
      
      if(first .eq. 1) then
         first = 0
         mevent = 1
         bmin = 0.
         bmax = -7.50D0           ! maximum impact parameter (1530mb)Au  
!            for pp collison. bmax is not used.?
!            for AA bmax=0 make loop
         dt = 100.d0
         nstep =1
         mstc(8) = 0       ! job mode
         mstc(156) = 0     ! no call to jamanacl
         mstc(161) = 0  ! no output on JAMINFO.DAT, JAMRUN.DAT
         mstc(38) = 0   ! stderr. If default 3 is used, open
        ! in jam.f (line 767) is called ; this "open" may
        ! result in a fatal error when the program runs under SGE.
        ! reason is unknown.
         mstc(39) = 0     ! no ouptut on JAMMULTI.DAT  
         mstc(81)=0  ! neglect hard scattering (for meson)
         frame = 'lab'
         fname(1) = '0'
         fname(2) = '/dev/null'  ! IO mstc(37)'s target.
!     fname(3) = '/dev/null'   ! IO mstc(38)'s target
                           ! Gfortran cannot open simult.
                 ! with mstc(37), so .
         fname(3)='/dev/stderr' 
      endif
      bmax = - sqrt( sig/3.141592/10. )  ! in fm  must be < 0
      call ccos2jamsymbl(pj, proj)  ! make projectile symbol
!               such as proj = 'k+'
      call ccos2jamsymbl(tg, targ)

!	targ = '14:7'
!	targ = 'p'
!	targ = 'n'
!	targ='197Au   '
!	targ = '4:2'
!	targ = '2:1'
      if(pj%code == kgnuc ) then
         write(cwin, '(f12.3,"gev")')  (pj%fm%p(4)-pj%mass)/pj%subcode
      else
         write(cwin, '(f12.3,"gev")')  (pj%fm%p(4)-pj%mass)
      endif
      call jaminit(mevent, bmin, bmax, dt, nstep, frame, proj, targ,
     *  cwin)
      end
      subroutine  ccos2jamsymbl(pj, symb)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      type(ptcl):: pj  ! input. code, subcode, charge 
                        ! iouput. incase of Lambda, pj is
                        ! modified to be a neutron
      character*8 symb  ! output. symbol used by jam

      integer nc

      if(pj%code .eq. kpion) then
         if(pj%charge .eq. 1) then
            symb = 'pi+'
         elseif(pj%charge .eq. 0) then
            write(0,*)
     *       'pi0 cannot be col. projectile in jam'
            stop
         else
            symb = 'pi-'
         endif
      elseif(pj%code .eq. kkaon) then
         if(pj%charge .eq. 1) then
            symb = 'k+'
         elseif(pj%charge .eq. -1) then
            symb = 'k-'
         else
           symb = 'k0'
!            write(0,*)  ' k0 cannot be proj. in jam'
!            stop
         endif
      elseif(pj%code .eq. knuc ) then
         if(pj%charge .eq. 1) then
            symb = 'p'
         elseif(pj%charge .eq.  -1 ) then
            symb = 'pbar'
         elseif(pj%subcode .eq. regptcl) then
            symb = 'n'
         else
            symb = 'nbar'
         endif
      elseif(pj%code .eq. kgnuc) then
         write(symb, '(i3,":", i3)' ) pj%subcode, 
     *       pj%subcode-pj%charge
!            remove blanks
         call kseblk(symb, "*", nc)  !  A:N say 207:125
      elseif( pj%code >= klambda .and. pj%code <= klambdac ) then
!          regard it as neutron                                                
         call cmkptc(knuc, -1, 0, pj)
         call cadjm(pj, pj)
         symb = 'n'
      elseif( pj%code == kbomega ) then
!         very rare: regard it as proton
         call cmkptc(knuc, -1, 1, pj)
         call cadjm(pj, pj)
         symb = 'p'
      else
         write(0,*) 'code=',pj%code, ' charge=',pj%charge
         write(0,*) ' not supported in jam'
         stop
      endif
      end
      subroutine cjam(pj, ia, iz,  a, ntp, 
     *     pjspec, npjspec, tgspec, ntgspec, inela, impactp)
!      implicit none  ! this cannot be used since jam1.inc uses
!      implicit double (...)
#include "jam1.inc"
#include "jam2.inc"
#include "Zptcl.h"
#include "Zcode.h"
#include "Ztrackp.h"

      type(ptcl)::pj  ! input incident ptcl 
      integer,intent(in)::ia  !  target A
      integer,intent(in)::iz  !  target Z
      type(ptcl)::a(*)
      integer,intent(out):: ntp ! # of non spectator ptcls stored in a
      type(ptcl)::pjspec(*) !output. projectile spector
      integer,intent(out):: npjspec ! # of ptcls in pjspec
      type(ptcl):: tgspec(*) !output  target spectator
      integer,intent(out):: ntgspec ! # of ptcls in tgspec
      integer,intent(out):: inela ! =0  elastic event
                                  ! =1  inelastic event
      real(8),intent(out):: impactp ! impact parameter used
      integer i
      integer code, subcode, charge
      integer iev/1/
      real(8)::pprobt, u, pprobp, probp
      type(ptcl):: aptcl
      integer::nn, hn, rn  ! # of free N, Heavy ions and Remnant
      save
      pprobt = iz
      pprobt = pprobt/ia   !  rough proton probabilty of target spectator
      if( pj%code == kgnuc ) then
         pprobp = pj%charge
         pprobp = pprobp/pj%subcode
         probp = (iz+pj%charge)
         probp = probp/(ia+pj%subcode)
      else
         probp = pprobt
      endif
      call jamevt(iev)
      ntp = 0
      inela = 0
      npjspec = 0
      ntgspec = 0
!////////////
!      tgTE =0.
!      pjTE = 0.
!///////////
!      if( ia == 1  ) then
!!         if target is p/n,
!         if( pj.code<=knuc .and. pj.code >= kpion  ) then
!!          and projectile is pi/K/p/n, inelastic coll. is
!!          defeined to be some meson production (or Lambda 0 etc 
!!            productin which is later judged)
!            if(nv > 2 ) then
!               inela = 1
!            endif
!         endif
!      endif

      do i = 1, nv
         if( k(1,i) <= 11 .and. k(1,i) >= 0  ) then  ! avoid dead ptcl 
!            if( ia /= 1 .or. pj.code >= kgnuc ) then  ! target is A or proj iS A 
!               if( abs(k(7,i)) > 1 ) then
!                  inela = 1
!               endif
!            endif
            call ckf2cos(k(2,i), code, subcode, charge)
            if( code /=  krare ) then 
               if( abs( k(7,i) ) > 1 ) then
                  inela = 1
               endif
!///////////////////
!               write(0,'(a, 4i4,2x,4i4)')
!     *          'k(7,i)=',k(7,i), code, subcode, charge,
!     *           pj.code, pj.subcode, pj.charge,  ia
!///////////
!                 target spectator treatment.  (original Jam 
!                 breaks  spectator heavy target  into nucoleons;
!                 and seems to be almost neutrons. We reasign the chage
!               Nara san's spec is -1 for target spectator but seems to be
!               incorrect and 1 must be used
               call cmkptc(code, subcode, charge, aptcl)
               aptcl%fm%p(1:4) = p(1:4,i)  ! jam is T.E  (manual==>K.E=> false).
!//////////////
!               write(0,'(a, 4i4, 1p,4g12.3)')
!     *            'k(7,=',k(7,i), code, subcode, charge,
!     *            p(1:4,i)
!/////////////
               if( k(7,i) == 1 ) then  ! target spec.
!                      re-assign the charge****NO NEED NOW****
!                  call rndc(u)
!                  if(u < pprobt ) then
!                     aptcl.charge = 1
!                  else
!                     aptcl.charge = 0
!                  endif
                  if( JamFragment == 0 ) then
                     ntp = ntp + 1
                     a(ntp) = aptcl
                  elseif( JamFragment == 1) then
                     ntgspec = ntgspec + 1
                     tgspec(ntgspec) = aptcl
                  else
                     write(0,*) ' JamFragment=',JamFragment,' invalid '
                     write(0,*) ' detected in cjam'
                     stop 00000
                  endif
               elseif( k(7,i) == -1) then  !proj. spec
!                     projectile spectator
!                    re-assign the charge; *****NO NEED NOW******
!                  call rndc(u)
!                  if(u < pprobp ) then
!                     aptcl.charge = 1
!                  else
!                     aptcl.charge = 0
!                  endif
                  if( JamFragment == 0 ) then
                     ntp = ntp + 1
                     a(ntp) = aptcl
                  else
                     npjspec = npjspec + 1
                     pjspec(npjspec) = aptcl
                  endif
               else
                  ntp = ntp + 1
                  a(ntp) = aptcl
!      **** NO NEED NOW****
!                     jam charge assignment for N is stragne, we randomly
!                     re-assign it
!                  if( code == knuc .and. subcode /= antip ) then
!                     call rndc(u)
!                     if(u < probp) then
!                        a(ntp).charge = 1
!                     else
!                        a(ntp).charge = 0
!                     endif
!                  endif
               endif
            endif
!            write(*,*) i, k(2,i), p(1,i), p(2,i), p(3,i),
!     *             p(4,i), p(5,i)
         endif
      enddo

!///////////////
!      if( pj.fm.p(4) > inci.fm.p(4)*0.98 ) then
!         write(*,*) 'pjTE=',pjTE, ' pj. =',pj.fm.p(4)
!         write(*,*) 'tgKE=',tgTE, ' Mt=', pjTE+tgTE-pj.fm.p(4)
!         write(*,*)
!         write(*,*) "====================="
!         write(*,*)
!      endif
!/////////////////
      impactp = pard(2) ! is in jam2.inc
      end
      subroutine  cjamPjSpec(pj, tg,
     *            pjspec, npjspec, impactp, inela, a, nout)
      use modTargetFrag
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      type(ptcl):: pj ! input. proj.
      type(ptcl):: tg ! //     taget
      type(ptcl)::pjspec(*) !input. projectile spector
      integer,intent(in):: npjspec ! # of ptcls in pjspec
      real(8),intent(in):: impactp ! impact parameter in fm of the interaction
      integer,intent(in):: inela   ! inela or ela event
      type(ptcl)::a(*)  ! output.  fragment from the spectator
      integer,intent(out):: nout ! # of ptcls in a

      integer::iat, izt  !  target A, Z

      integer:: pjA, pjZ 
      integer:: nw  ! woonded N
      integer:: oA(209),  oZ(209)
      integer:: i, j, n
      integer:: nn, hn, rn
      integer,save:: nwrong=0

      logical:: break

      break = .true.
      nout = 0
      n = 0
      if(tg%code == kgnuc) then
         iat = tg%subcode
      else
         iat = 1
      endif
      izt = tg%charge
      if( inela == 0 .and. pj%code /= kgnuc ) then  !  elastic event. proj. is only 1  spectator
         nout = nout + 1
         a(nout) = pjspec(1)
      elseif( iat == 1) then  ! fragm cannot treat this case
         break=.false.
      elseif( npjspec > 0 ) then  
         pjA = pj%subcode
         if( pjA  > 56 .and. iat > 56) then
            call csampNucFrag(pj, iat, izt, npjspec,
     *             a, nout)
         elseif(pjA > 56 ) then  ! iat<=56
            !        use NHR by exchaning pj and target
            call cuseFragNHR(iat, izt,  pjA, pjZ, pjspec,  npjspec,
     *            oA, oZ, a, nout)
         else
            nw = pj%subcode - npjspec
            !    use fragm (target comes 1st)
            call cusefragm(iat, pjA, nw, pjspec,
     *           impactp, oA, oZ,   a, nout)
         endif  
      endif
      if( .not. break ) then
!              spectator nucleons do not breakup. so 
!              make them an ion
         if( npjspec > 0 ) then 
            nout = nout + 1
            call cjamNoFrag(pj, npjspec, pjspec, a(nout))
         endif
      endif
      end
      
      subroutine cuseFragNHR(pjA, pjZ, iat, izt, tgspec,  ntgspec, 
     *            oA, oZ, a, nout)
      use modTargetFrag
!           this is to sample fragmentation particles from
!           target.  target  A must be  >= 56.  proj A <=56
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      integer,intent(in):: pjA, pjZ  ! projectil A,Z  A<=56
      integer,intent(in)::iat, izt  !  target A,Z  A>=56
      type(ptcl)::tgspec(*)  ! input.  spectator nucleon info.
      integer,intent(in):: ntgspec ! # of target spectator in tgspec
      integer,intent(out):: oA(209),  oZ(209) ! working array
      type(ptcl)::a(*)  ! output.  fragment from the target
      integer,intent(out):: nout ! # of ptcls in a


      integer:: i, j, n
      integer:: nn, hn, rn

      call cSampFragNHR(pjA, pjZ, iat, izt,  ntgspec,
     *     nn, hn, rn, oA, oZ)
      n = 0
      nout =  nn+hn+rn
      do i = 1, nout
         if( oA(i) == 1 ) then
            n = n + 1           ! n/p
            a(i) =  tgspec(n)
            a(i)%charge = oZ(i)
         else
            a(i)%fm%p(1:4) = 0.   
            do j = 1, oA(i)
               n = n + 1
               a(i)%fm%p(1:4) = a(i)%fm%p(1:4)
     *              + tgspec(n)%fm%p(1:4)
            enddo
            call cmkptc(kgnuc, oA(i), oZ(i), a(i))
            a(i)%mass = sqrt( a(i)%fm%p(4)**2 
     *           -   dot_product(a(i)%fm%p(1:3), a(i)%fm%p(1:3)) )
         endif
      enddo
      end subroutine cusefragNHR

      subroutine cusefragm(iat, pjA, nw, pjspec, 
     *   impactp, oA, oZ, a, nout)
      use modfragment, only: fragm
!            this is to call fragm routine
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      integer,intent(in):: iat  !  target A
      integer,intent(in):: pjA ! projectile A
      integer,intent(in):: nw  ! woonded N
      type(ptcl)::pjspec(*)  ! input projectile spectator info
      real(8),intent(in):: impactp ! impact parameter in fm of the interaction
      integer,intent(out):: oA(209),  oZ(209)
             !  this is used inside of this prog. 
      integer,intent(out):: nout ! # of fragment in oA, oZ and a 
      type(ptcl):: a(*)  ! output.  fragment from the spectator


      integer:: i, j,  n
      n = 0
!              target comes first 
      call fragm(iat, pjA, nw, impactp, nout, oA, oZ)
      do i = 1, nout
         if( oA(i) == 1 ) then
            n = n + 1           ! n/p
            a(i) =  pjspec(n)
            if( oZ(i) > 1 ) then
               oZ(i) = 1  ! never seen, but for safety.  
            endif
            a(i)%charge =  oZ(i)
         else
            a(i)%fm%p(1:4) = 0.   
            do j = 1, oA(i)
               n = n + 1
               a(i)%fm%p(1:4) = a(i)%fm%p(1:4)
     *              + pjspec(n)%fm%p(1:4)
            enddo
            call cmkptc(kgnuc, oA(i), oZ(i), a(i) )
!!            call cmkptc(kgnuc, oA(i), 1, a(i))
            a(i)%mass = sqrt( a(i)%fm%p(4)**2 
     *           -   dot_product(a(i)%fm%p(1:3), a(i)%fm%p(1:3)) )
            if(oA(i) == 2) then ! d
               a(i)%charge = 1
            elseif( oA(i) <= 4 ) then ! He3 or He4
               a(i)%charge = 2
            elseif( oA(i) < 20 ) then
               a(i)%charge = max(3.0, oA(i)*0.5)
            else
               a(i)%charge = oA(i)*0.4
            endif
         endif
      enddo
      end subroutine cusefragm

      subroutine cjamNoFrag(p0, nspec, spec,  a)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      type(ptcl):: p0 ! input. nucleus before interaction
      integer,intent(in):: nspec ! among p0, nspec is spectator nuc.
      type(ptcl):: spec(nspec) ! inut. spectator nucleons
      type(ptcl)::a   ! output. npec nucleon is made to be a nucleus

      integer Zsum 
      integer i
      a%charge = 0
      a%fm%p(:) = 0
      if( nspec > 1 ) then
         do i =1, nspec
            a%charge = a%charge + spec(i)%charge
            a%fm%p(:) = a%fm%p(:) + spec(i)%fm%p(:)
         enddo
         a%subcode = nspec
         a%code = kgnuc
         if( a%charge <= 0 ) then 
!               very very rare; happend when target is H
!              and pj is He.  
!            write(0,*) 'In cjamNoFrag.  p0=',p0.fm.p(1:4) 
!            write(0,*) 'nspec =', nspec
!            write(0,*) ' charge sum=',0
!        if Z=0 is kept, dpmjet will die by NEGATIVE ENERGY.
            a%charge = 1
!              next is for safety
         elseif( a%charge >= nspec  .and. nspec > 1 ) then
            a%charge = max(nspec/2, 1)
         endif
         a%mass =  a%fm%p(4)**2 -
     *    dot_product(a%fm%p(1:3), a%fm%p(1:3))
!         if(a.mass < 0. ) then
!            write(*,*) 'errp mass2',   a.mass
!         endif
         a%mass = sqrt(a%mass)
      else
         a = spec(1)
      endif
      end 

      subroutine cjamTgSpec(pj, tg,
     *            tgspec, ntgspec, impactp, inela, a, nout)
      use modTargetFrag
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      type(ptcl):: pj ! input. proj.
      type(ptcl):: tg ! input. target.
      type(ptcl)::tgspec(*) !input. target spectator
      integer,intent(in):: ntgspec ! # of ptcls in tgspec
      real(8):: impactp  ! impact parameter in fm
      integer,intent(in):: inela  ! ela/inela flag
      type(ptcl)::a(*)  ! output.  fragment from the spectator
      integer,intent(out):: nout ! # of ptcls in a

      integer::iat, izt  !  target A,Z

      integer:: pjA, pjZ 
      integer:: oA(209),  oZ(209)
      integer:: nw  ! # of woonded target nucleons
      integer:: i, j, n
      integer::nn, hn, rn  ! # of free N, Heavy ions and Remnant
!       if target 1 < A <=56 use fragm by interchanging proj and target
!       if target A>56 use cSampFragNHR by using proj and targ as they are
      integer,save:: nwrong=0
      logical:: break
      if( tg%code == kgnuc) then
         iat = tg%subcode
      else
         iat = 1
      endif
      izt = tg%charge
      if(pj%code == kgnuc ) then
         pjA = pj%subcode
      else
         pjA=1
      endif
      pjZ = pj%charge

      break = .true.
      nout = 0
      n = 0
      if( inela == 0 .and. iat == 1 .and. ntgspec > 0  ) then  ! p/n target and elastic
         nout = nout + 1                      !   inela == 0 is redunant
         a(nout) = tgspec(1)
      elseif( ntgspec > 0 ) then  
         if( iat > 56 .and. pj%subcode > 56) then
            call csampNucFrag(tg, int(pj%subcode),
     *           int(pj%charge), ntgspec, a, nout)
         elseif( iat > 56 ) then !  pjA <=56
            call cuseFragNHR(pjA, pjZ, iat, izt, tgspec,  ntgspec,
     *            oA, oZ, a, nout)            
         elseif(  iat > 1  ) then
            nw = iat - ntgspec
            if(pjA == 1 ) then  !  fragm cannot treat this case
               break = .false.
            else
!                     target (=pjA) comes first
               call cusefragm(pjA, iat, nw, tgspec,
     *           impactp,  oA, oZ,   a, nout)
            endif
         else
            break = .false.
         endif
         if(.not. break ) then
            if( ntgspec > 1 ) then 
               nout = nout + 1
               call cjamNoFrag(tg, ntgspec, tgspec, a(nout))
            endif
         endif
      endif
      end
