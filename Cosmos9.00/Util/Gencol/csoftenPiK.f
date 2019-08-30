
      module SoftenPiK
!         This part is common to Cosmos(modifyX) and Cosmos/Util/Gencol
!     Used to soften the x-distribution.     
      implicit none

!         The parameters are listed above the --------------- line
! mode=0)      don't do any cut/modification of X
! mode=1)      discard all pi/K with x>xth. 
!              However, if the incident is pi/K the highest
!              energy pi/K, resp, is untouched.
! mode=2)      discard only pi0 with x> xth
! mode=3)      soften pi/K/eta  spectrum; 
!              If their x is >  x>xth, this is applied.
!              (except for the leading pi/K when
!              incident is pi/K).
! mode=-3      probability of softening depends on the user's algorithm
!              The user must modify cJudegeSplit.
!    trial version method          
!              If we find a pi/K/eta with x>xth,
!              apply softening with the prob. of (x-xth)**a (a=0.1 default).
!              Softening is done by splitting  the particle into two by
!              calling csplitMeson.  The sum of two particle is made to
!              be the same as the original particle. However, momentum 
!              cannot be conserved: the two particle will have the same
!              Pt as the original one and 3 memeta are adjuested to 
!              satisfy P^2 + m^2 =E^2. Let the two particles' energy be
!              E1 and E2.  If E2>E1, we split E2 into two again.
!             (At present, if E1> E2, we don't  split E1). 
!              (so, x>xth pi/K/eta may split into two or three).
!            This method results in some unnatural small bump at x<xth.
!

      integer,save:: mode = 3      
      real(8),save:: xth = 0.05    ! over this x, we apply cut  or modification
      real(8),save:: pw = 0.5 ! softening is determined by
!                  if(  (x-xth)**pw  < u  ) cycle  
!                  call csplitMeson(pstore(i), E1, E2, icon)
!                  so if  pw in (x-xth)**pw is smaller, stronger softening
      real(8),save:: repeat = 2.5  ! apply the simple softenning repeat times
                !   for one event if the event has pi/K/eta with X> Xth. 
                !   the odd number is probabilistic. If negative,
                !   positive number is taken  and understood as poisson
                !   average.

      real(8),save:: E0th = 500.   !  500 GeV  over this we apply cut/mod.
                           !  if make this 10 TeV, effect after shower max
                           !  decreases even for  E1ry=10^19 eV.  
                           !  This is in lab. frame
      integer,save:: fwbw = 3    ! Used in Gencol; modification is applied to
                                 ! 1--> forward only 2--> backword only
                                 ! 3--> both ; However, 1 is used when
                     ! Cosmos output (in cms) is the target, irespectively of
                     ! fwbw.  
      logical,save:: special=.false.  ! This is used olny in modifyX of
                     !   Cosmos. but not used in Gencol.
                     !   In Cosomos, it is used to see the x-dist. of
                     !   pi, K, eta  at the first interaction; 
                     ! For that, make this  .ture. and 
                     ! ***set next in parm*** (don't forget to restore
                     ! them to the default if special=.false. is reset.)
                     !  =====================================
                     !  Generate="em" (don't use "em/as")
                     !  EminObs = 8*xxxx,
                     !  =====================================
                     !  where xxxx is a value (GeV) little bit smaller then
                     !  the 1ry  energy. E.g 0.99999e8 if primary is 10^17 eV.
                     !
                     ! output will look like
                     ! xd  4  1  12.3e-4 
                     ! xd  4  0  2.89e-1
                     !  where xd is id, next two are code and charge
                     !  last one is X. in lab. frame.
             ! For special=.false., i.e, to generate air showers,
             ! you should have
             !  =====================================
             !  Generate="em/as" 
             !  EminObs = 500e-6,
             !  =====================================
        
      logical,save:: leadingPiK=.true.  ! for Pi/K incident case,
                                ! treat highest energy Pi/K of the same
                                ! charge as leading Pi/K and don't 
                                ! apply the method here
                  ! If used in Cosmos and 
                  !  if make this f, effect increase a bit esp. at
                  !  after max.  
                  !  leadingPiK=T --> Xmax ==> 35 g/cm^2
                  !            =F     Xmax ==> 25 g/cm^2
      logical,save:: useXinCMS=.false.  !This is used only in Cosmos.
                !  In Cosmos, generated particle's E is in Lob. so 
                !  x= E/E0 is also in Lob.  If this is .false., softening by 
                !  csoftenPiK is applied to the Lab.X. 
                !  if this is .true.,  softening is performed after
                !  converting the energy into cms system, and then
                !   re-boosted to Lab.  
                !  
      real(8),save:: k1u=0.25    !   
      real(8),save:: xthl=1.d0
      real(8),save:: xthu=10.d0
      real(8),save:: rejpw=0.
      logical,save:: modified   ! not input. used in Cosmos
      real(8),save:: E0lab=1.e12   ! not input. current E0 in lab. 
                                ! frame. Used in Cosmos

!             next one is about 0.6 of maxn in Epics/Util/Gencol/Zprivate.h
!             but this  cannot inherit maxn
      integer,parameter::half=12000  
!--------------------------------------------------------------------

      contains

      subroutine cJudgeSplit(pj, x,  split) 
      implicit none
!         judge if this x should be split
#include "Zptcl.h"
      type (ptcl):: pj    ! input. projectile. you may use this
      real(8),intent(in):: x  ! KE ratio
      logical,intent(out):: split  ! .true.  or .false.. if true, split meson

      real(8):: u

      if( mode == 3) then
         call rndc(u)   ! uniform rn in (0,1)
!          if pw in (x-xth)**pw is smaller, stronger softening
         split = (x-xth)**pw  > u  
      elseif(mode == -3 ) then
!            comment out next lines and supply your judgement
!         and fix split (.true. or .false.)
!         you may use xth, pw etc which are defined before "contains"
         write(0,*) ' you have to supply your owon code here'
         write(0,*) ' to judge if this particle should be '
         write(0,*) ' split or not' 
         write(0,*)
     *    ' Place is cJudgeSplit; first subroutine in  csoftenPiK'
         stop
      endif
         
      end     subroutine cJudgeSplit
      
      subroutine csoftenPiK(inciptcl, pstore, nin, nout)
!         inciptcl made a collsion and generated nin particles in pstore
!         particle information may be at CMS or Lab.
!     If we want to  make the x-dist. softer,
!      A) For Cosmos (air target),  we may use the Lab system
!         (for large X, the distribution would be the same 
!           as CMS)
      implicit none
#include  "Zcode.h"
#include  "Zptcl.h"
#include  "Zmanagerp.h"
      type (ptcl)::inciptcl  ! input. incident ptcl at collision
!                  
      integer,intent(in)::nin ! # of ptcls in pstore 
      type (ptcl)::pstore(*)  ! input/output.   ptcls to be softened.
                           !     Softened ones are put here
      integer,intent(out)::nout ! # of ptcls re-stored in pstore 
      logical,save::first=.true.
      integer:: jcon, i
      integer:: nrepeat
      real(8):: u
      integer:: exist
      integer:: nnow
      
      if( first ) then
         call copenf(TempDev, 
     *       "$COSMOSTOP/UserHook/modifyX/softenparam%dat", jcon)
         if( jcon /= 0 ) then
            write(0,*) 'Data file for csoftenPiK could not be found'
            write(0,*) 'in $COSMOSTOP/Util/Data/softenpiK%dat'
            write(0,*) 'Forgot to set COSMOSTOP ?'
            stop
         else
            call creadSoftenPara(TempDev)
         endif
         call cwriteSoftenPara(0)
         first = .false.
      endif
      if( repeat < 0. ) then
         call kpoisn(-repeat, nrepeat)
      else
         call rndc(u)
         nrepeat = repeat
         if( u < repeat-nrepeat)  then
            nrepeat = nrepeat + 1
         endif
      endif

      nnow = nin
      xth = csoftenFixXth(E0lab)
      k1u = ck1u(E0lab)
      xthu =cxthu(E0lab)
      do i = 1, nrepeat
         call csoftenPiK0(inciptcl, pstore, nnow, nout, exist)
         if( exist == 0 ) exit 
         nnow = nout
      enddo
      do i = 1, 1
         call csoftenPiK1(inciptcl, pstore, nnow, nout, exist)
         if( exist == 0 ) exit 
         nnow = nout
      enddo
      end   subroutine csoftenPiK


      subroutine csoftenPiK0(inciptcl, pstore, nin, nout, exist)
!         inciptcl made a collision and generated nin particles in pstore
!         particle.  information may be at CMS or Lab.
!     If we want to  make the x-dist. softer,
!      A) For Cosmos (air target),  we may use the Lab system
!         (for large X, the distribution would be the same 
!           as CMS)
      implicit none
#include  "Zcode.h"
#include  "Zptcl.h"

      type (ptcl)::inciptcl  ! input. incident ptcl at collision
!                  
      integer,intent(in)::nin ! # of ptcls in pstore 
      type (ptcl)::pstore(*)  ! input/output.   ptcls to be softened.
                           !     Softened ones are put here
      integer,intent(out)::nout ! # of ptcls re-stored in pstore 
      integer,intent(out)::exist  ! # of particcls with X> Xth
! ----------------------------------------------
      type (ptcl)::work(2)
      logical  ok
      integer::maxi
      real(8)::maxE
      integer::i, j 
      real(8):: x, E0, Ex, u, u1, u2     
      real(8):: E1, E2
      real(8)::temp
      integer:: nc  ! counting # of ptcls in pstore
      integer:: nw, icon
      logical:: split 

      E0 = inciptcl%fm%p(4) - inciptcl%mass
      nout = nin
      exist = 0
      if( E0lab < E0th) return      !!!!!
      nc = nin
!            if pi/K is incident, regards the highest one as leading
!            and avoid to modify 
      if(leadingPiK .and. (inciptcl%code == kpion .or.
     *     inciptcl%code == kkaon ) ) then
!              find max energy ptcl index and E; eta cannot be incident
         call cgetmaxptcl(pstore, nc, inciptcl%code, inciptcl%charge,
     *        maxi, maxE)
      else
         maxi = 0
      endif

      do i = 1, nin
         x = (pstore(i)%fm%p(4)-pstore(i)%mass)/E0
         if( abs(mode) == 3 ) then
            if( i /= maxi ) then
!                softening
               if(x > xth ) then
!                  modify high X of pi/K/eta
                  if(pstore(i)%code == kpion .or.
     *                 pstore(i)%code == kkaon 
     *                 .or.   pstore(i)%code ==  keta ) then
                     nw = 0
                     exist = exist + 1
!                       judge to split this
                     call cJudgeSplit(inciptcl, x,  split)
                     if( .not. split ) cycle   

                     call csplitMeson(pstore(i), E1, E2, icon)
                     if(icon == 0) then
                        modified = .true.
                     ! split;  if impossilbe, do nothing
                        nw = nw + 1
                        work(nw) = pstore(i)
                        work(nw)%fm%p(4) = E1 + work(nw)%mass 
                        call cae2p(work(nw)) ! adjust momentum
                        nw = nw + 1
                        work(nw) = pstore(i)
                        work(nw)%fm%p(4) = E2 + work(nw)%mass 
                        call cae2p(work(nw))
!                       original one is replaced by E1
                        pstore(i) = work(1)
!                       others are appended to  pstore. 
                        do j = 2, nw ! altough always nw=2
                           nc = nc + 1
                           pstore(nc) = work(j)
                        enddo
                     endif
                  endif
               endif
            endif
         elseif( mode == 1 ) then
            if( i /= maxi ) then
!               discard all pi/K/eta with x>xth (except leanding)
               if(pstore(i)%code == kpion .or. pstore(i)%code == kkaon 
     *              .or.   pstore(i)%code ==  keta ) then
                  if( x > xth ) then
!                   to neglect, put mass as E
                     pstore(i)%fm%p(4) = pstore(i)%mass
                     pstore(i)%fm%p(1:3) = 0.
                     modified = .true.
                  endif
               endif
            endif
         elseif( mode == 2 ) then 
            if( i /= maxi) then
!              discard only pi0/eta
               if( ( pstore(i)%code == kpion .and.
     *              pstore(i)%charge  == 0 ) .or.
     *              pstore(i)%code == keta ) then
                  if(x> xth) then
                     pstore(i)%fm%p(4) = pstore(i)%mass
                     pstore(i)%fm%p(1:3) = 0.    ! put zero energy
                     modified = .true.
                  endif
               endif
            endif
         elseif( mode == 0 ) then
!               do nothing
         else
            write(0,*) ' mode err=',mode, ' in csoftenPiK0'
            stop
         endif
      enddo
      if(abs(mode) == 3) then
         nout = nc
      endif
      end subroutine csoftenPiK0


      subroutine csoftenPiK1(inciptcl, pstore, nin, nout, exist)
!         inciptcl made a collision and generated nin particles in pstore
!         particle.  information may be at CMS or Lab.
!     If we want to  make the x-dist. softer,
!      A) For Cosmos (air target),  we may use the Lab system
!         (for large X, the distribution would be the same 
!           as CMS)
      implicit none
#include  "Zcode.h"
#include  "Zptcl.h"

      type (ptcl)::inciptcl  ! input. incident ptcl at collision
!                  
      integer,intent(in)::nin ! # of ptcls in pstore 
      type (ptcl)::pstore(*)  ! input/output.   ptcls to be softened.
                           !     Softened ones are put here
      integer,intent(out)::nout ! # of ptcls re-stored in pstore 
      integer,intent(out)::exist  ! # of particcls with X> Xth
! ----------------------------------------------
      type (ptcl)::work(2)
      logical  ok
      integer::maxi
      real(8)::maxE
      integer::i, j 
      real(8):: x, E0, Ex, u, u1, u2     
      real(8):: E1, E2
      real(8)::temp
      integer:: nc  ! counting # of ptcls in pstore
      integer:: nw, icon
      logical:: split 

      real(8):: xcent, centval
      xcent = sqrt(xth/xthu*xth*xthl)
      centval =  crejK1(xcent)


      E0 = inciptcl%fm%p(4) - inciptcl%mass
      nout = nin
      exist = 0
      if( E0lab < E0th) return      !!!!!
      nc = nin
!            if pi/K is incident, regards the highest one as leading
!            and avoid to modify 
      if(leadingPiK .and. (inciptcl%code == kpion .or.
     *     inciptcl%code == kkaon ) ) then
!              find max energy ptcl index and E; eta cannot be incident
         call cgetmaxptcl(pstore, nc, inciptcl%code, inciptcl%charge,
     *        maxi, maxE)
      else
         maxi = 0
      endif

      do i = 1, nin
         x = (pstore(i)%fm%p(4)-pstore(i)%mass)/E0
         if( x > xth*xthl ) cycle
         if( abs(mode) == 3 ) then
            if( i /= maxi ) then
!                softening
               if(x > xth/xthu) then
!                  modify high X of pi/K/eta
                  if(pstore(i)%code == kpion .or.
     *                 pstore(i)%code == kkaon 
     *                 .or.   pstore(i)%code ==  keta ) then
                     nw = 0
                     exist = exist + 1
!                       judge to split this
!!                     call cJudgeSplit(inciptcl, x,  split)
                     call rndc(u)
                     temp = crejK1(x)
!                     split = u < 0.25*(temp/centval)
!                     split = u < 0.125*(temp/centval)
                     split = u < k1u*(temp/centval)
                     if( .not. split ) cycle   

                     call csplitMeson(pstore(i), E1, E2, icon)
                     if(icon == 0) then
                        modified = .true.
                     ! split;  if impossilbe, do nothing
                        nw = nw + 1
                        work(nw) = pstore(i)
                        work(nw)%fm%p(4) = E1 + work(nw)%mass 
                        call cae2p(work(nw)) ! adjust momentum
                        nw = nw + 1
                        work(nw) = pstore(i)
                        work(nw)%fm%p(4) = E2 + work(nw)%mass 
                        call cae2p(work(nw))
!                       original one is replaced by E1
                        pstore(i) = work(1)
!                       others are appended to  pstore. 
                        do j = 2, nw ! altough always nw=2
                           nc = nc + 1
                           pstore(nc) = work(j)
                        enddo
                     endif
                  endif
               endif
            endif
         elseif( mode == 1 ) then
            if( i /= maxi ) then
!               discard all pi/K/eta with x>xth (except leanding)
               if(pstore(i)%code == kpion .or. pstore(i)%code == kkaon 
     *              .or.   pstore(i)%code ==  keta ) then
                  if( x > xth ) then
!                   to neglect, put mass as E
                     pstore(i)%fm%p(4) = pstore(i)%mass
                     pstore(i)%fm%p(1:3) = 0.
                     modified = .true.
                  endif
               endif
            endif
         elseif( mode == 2 ) then 
            if( i /= maxi) then
!              discard only pi0/eta
               if( ( pstore(i)%code == kpion .and.
     *              pstore(i)%charge  == 0 ) .or.
     *              pstore(i)%code == keta ) then
                  if(x> xth) then
                     pstore(i)%fm%p(4) = pstore(i)%mass
                     pstore(i)%fm%p(1:3) = 0.    ! put zero energy
                     modified = .true.
                  endif
               endif
            endif
         elseif( mode == 0 ) then
!               do nothing
         else
            write(0,*) ' mode err=',mode, ' in csoftenPiK0'
            stop
         endif
      enddo
      if(abs(mode) == 3) then
         nout = nc
      endif
      end subroutine csoftenPiK1

      subroutine csplitMeson(p, E1, E2, icon) 
      implicit none
#include  "Zcode.h"
#include  "Zptcl.h"
      type (ptcl)::p   ! input, a high energy ptcl
      real(8),intent(out):: E1  ! split ptcl energy KE
      real(8),intent(out):: E2  ! split ptcl
      integer,intent(out):: icon ! 0--> split ok, 1--> no split
      real(8):: u,  Emin, Em
      logical ok
      integer::count

      ok = .false.

      count = 0
      do while(.not. ok)
         call rndc(u)
         u = u*(1.-xth) + xth
         E1 = u*( p%fm%p(4) - p%mass)
         E2 = (p%fm%p(4) -p%mass) - E1
!         if(E1 > p.mass .and.
!     *      E2 > p.mass ) then
            ok = .true.
!         else
!            count = count + 1
!            if(count > 20) then
!               icon = 1
!               return
!            endif
!         endif
      enddo
      icon = 0
      end subroutine csplitMeson
 
      subroutine  cae2p( pc )
!             adjust mpmentum by refering to changed E
!         keep the Pt same, if possible
      implicit none
#include  "Zcode.h"
#include  "Zptcl.h"


      type (ptcl):: pc
      
      real(8):: E,  P2, cf

      E = pc%fm%p(4)
!      
!         Pt^2 + Pz^2 +m^2= E^2
!         so new Pz =sqrt( E^2-Pt^2-m^2) 
!          the sign is the same as original one
      P2 =E**2- ( pc%fm%p(1)**2  + pc%fm%p(2)**2 + pc%mass**2)
      if( P2 > pc%mass**2  ) then
         pc%fm%p(3) =  sign(sqrt(P2), pc%fm%p(3))
      else
!         keep the direction and shirnk the magnitude of p
         P2 = pc%fm%p(1)**2  + pc%fm%p(2)**2  + pc%fm%p(3)**2
         if( E <= pc%mass .or. P2 == 0. ) then
            pc%fm%p(1:3) = 0.         
         else
            cf = sqrt( (E**2 - pc%mass**2) / P2 )
            pc%fm%p(1:3) = pc%fm%p(1:3)*cf
         endif
      endif
      end  subroutine cae2p 
      subroutine cgetMaxptcl(pstore, nin,  pcode, pcharge, maxi, maxE)
!        get max energy ptcl with the same code / charge as incident
!        if meson is incident, most probably, it is leading.
      implicit none
#include  "Zcode.h"
#include  "Zptcl.h"

      integer,intent(in):: nin ! # of ptcls in pstore
      type (ptcl):: pstore(nin) 
      
      integer(2),intent(in):: pcode    ! incident code
      integer(2),intent(in):: pcharge  ! //       charge
      integer,intent(out):: maxi  !   index of maxE in pstore  // 
      real(8),intent(out):: maxE   ! max Energy with the same code/charg as
                              ! incident.  if there is no such, 0


      integer i

      maxi = 0
      maxE = 0.
      do  i = 1,  nin
         if( pstore(i)%code /=  pcode ) cycle
         if( pstore(i)%charge /= pcharge ) cycle
         if( maxE < pstore(i)%fm%p(4) ) then
            maxE =  pstore(i)%fm%p(4) 
            maxi = i
         endif
      enddo
      end subroutine cgetMaxptcl

      subroutine csoftenFE(inci, fwbwin, a, nin, nout)
!       front end of softening when it is to be  done at CMS
!         x is defined simply by E/E0, 
      implicit none 
#include "Zptcl.h"      
      type (ptcl):: inci  ! incident ptcl (in cms); symmetric case 
      integer,intent(in):: fwbwin  ! modification is applied to
                                 ! 1--> forward only 2--> backword only
                                 ! 3--> both ; However, 1 is used when
                     ! Cosmos output (in cms) is the target, independently of
                     ! fwbw.  
      type (ptcl):: a(*)  ! array containing ptcl info.
      integer,intent(in)::nin ! # of ptcls in w
      integer,intent(out)::nout ! # of ptcls put in w after modification

      type (ptcl):: w1(half)  ! working array
      type (ptcl):: w2(half)  ! working array

      integer::i, nc1, nc2, ncout1, ncout2

!       do modification  extract Pz>0  and Pz<0
      nc1 = 0
      nc2 = 0
      do i = 1, nin
         if(a(i)%fm%p(3) > 0.) then
            nc1 = nc1 + 1
            w1(nc1) = a(i)
         else
            nc2 = nc2 + 1
            w2(nc2) = a(i)
         endif
      enddo
      if( nc1 > 0 .and. IBITS(fwbwin,0,1)>0 ) then  ! LSB pos=0
         call csoftenPiK(inci, w1, nc1, ncout1)
      else
         ncout1 = nc1
      endif
      if( nc2 > 0 .and. IBITS(fwbwin,1,1)>0 ) then  ! 2nd bit pos=1
         call csoftenPiK(inci, w2, nc2, ncout2)
      else
         ncout2 = nc2
      endif

      nout = 0
      do i = 1, ncout1
         nout = nout + 1
         a(nout) = w1(i)
      enddo
      do i = 1, ncout2
         nout = nout + 1
         a(nout) = w2(i)
      enddo
      end subroutine csoftenFE

      subroutine creadSoftenPara(io)
      implicit none
      integer,intent(in):: io   ! logical dev. #
      character*24 vname
      character*100 vvalue


       call cskipsep(io)
       do while( cgetParmN(io, vname, vvalue ) )
          select case(vname)
          case('mode')
             call creadParaI(vvalue, mode)
          case('xth') 
             call creadParaR(vvalue, xth)
          case('E0th') 
             call creadParaR(vvalue, E0th)
          case('fwbw')
             call creadParaI(vvalue, fwbw)
          case('pw')
             call creadParaR(vvalue, pw)
          case('repeat')
             call creadParaR(vvalue, repeat)
          case('special')
             call creadParaL(vvalue, special)
          case('leadingPiK')
             call creadParaL(vvalue, leadingPiK)
          case('useXinCMS')
             call creadParaL(vvalue, useXinCMS)
          case('k1u')
             call creadParaR(vvalue, k1u)
          case('xthl')
             call creadParaR(vvalue, xthl)
          case('xthu')
             call creadParaR(vvalue, xthu)
          case('rejpw')
             call creadParaR(vvalue, rejpw)
          end select
       enddo
       end       subroutine creadSoftenPara
!      *************
       subroutine  cwriteSoftenPara(io)
       implicit none
       integer,intent(in):: io

       write(io,*)'----------------------'
       call cwriteParaI(io,'mode', mode)
       call cwriteParaR(io,'xth', xth)
       call cwriteParaR(io,'E0th', E0th)
       call cwriteParaI(io,'fwbw', fwbw)
       call cwriteParaR(io,'pw', pw)
       call cwriteParaR(io,'repeat',repeat)
       call cwriteParaL(io,'special',special)
       call cwriteParaL(io,'leadingPiK',leadingPiK)
       call cwriteParaL(io,'useXinCMS', useXinCMS)
       call cwriteParaR(io,'k1u', k1u)
       call cwriteParaR(io,'xthl', xthl)
       call cwriteParaR(io,'xthu', xthu)
       call cwriteParaR(io,'rejpw', rejpw)
       end       subroutine  cwriteSoftenPara


       subroutine cskipsep(io)
       implicit none
       integer io
       character(10)  sep
       do while (.true.)
          read(io, '(a)') sep
          if(sep(2:10) == '---------') exit
       enddo
       end  subroutine cskipsep
!        ************************* real*8 data
       subroutine creadParaR(vvalue, x)
        implicit none
        integer io
        character*(*) vvalue
        real*8 x
!        read(vvalue, *)   x, x
        read(vvalue, *)   x
        end       subroutine creadParaR
       subroutine creadParaR2(vvalue, x, n)
        implicit none
        integer io
        character*(*) vvalue
        integer n
        real*8 x(n)
        read(vvalue, *)   x
        end       subroutine creadParaR2

!     ************************* complex data
      subroutine creadParaCx(vvalue, c)
      implicit none
      character*(*) vvalue
      complex*8 c
      read( vvalue, *)   c
      end      subroutine creadParaCx
!     ************************ integer data
      subroutine creadParaI(vvalue, i)
      implicit none
      character*(*) vvalue
      integer i
      read(vvalue, *)   i
      end      subroutine creadParaI
!        ************************* character data
      subroutine creadParaC(vvalue, cha)
      implicit none
      character*(*) vvalue
      character*(*) cha
      read(vvalue, *)  cha
      end      subroutine creadParaC
!        ***************************** logical data
      subroutine creadParaL(vvalue, logi)
      implicit none
      character*(*) vvalue
      logical logi
      read(vvalue, *)  logi
      end           subroutine creadParaL
!        ---------------------------------------------
      subroutine cwriteParaR(io, vname, x)
      implicit none
      integer io
      character*(*) vname
      real*8  x
      
      write(io, *) ' ', vname,' ', x,' /'
      end      subroutine cwriteParaR
      subroutine cwriteParaR2(io, vname, x, n)
      implicit none
      integer io
      integer n  ! arra size of x
      character*(*) vname
      real*8  x(n)
      
      write(io,*) ' ', vname,' ', x,' /'
      end      subroutine cwriteParaR2

      subroutine cwriteParaCx(io, vname, c)
      implicit none
      integer io
      character*(*) vname
      complex*8  c
      write(io,  *) ' ', vname,' ', c,' /'
      end subroutine cwriteParaCx

      subroutine cwriteParaI(io, vname, i)
      implicit none
      integer io
      character*(*) vname
      integer i
      
      write(io,  *) ' ', vname,' ', i,' /'
      end      subroutine cwriteParaI

      subroutine cwriteParaC(io, vname, cha)
      implicit none
      integer io
      character*(*) vname
      character*(*) cha
      integer klena
      character*2 qmk/" '"/             ! ' 
      if(klena(cha) .gt. 0) then
         write(io,  *) ' ', vname, qmk, cha(1:klena(cha)),
     *        qmk,' /'
      else
         write(io, *) ' ', vname, qmk, ' ', qmk, ' /'
      endif
      end      subroutine cwriteParaC
      subroutine cwriteParaL(io, vname, logi)
      implicit none
      integer io
      character*(*) vname
      logical  logi

      write(io,  *) ' ', vname,' ', logi,' /'
      end      subroutine cwriteParaL

      function crejK1( x ) result(ans)
      implicit none
      real(8),intent(in):: x
      real(8):: ans
      ans = ( log(x/(xth/xthu)) * log((xth*xthl)/x) )**rejpw
      end  function crejK1

      function cxthu( E0 ) result(ans)
      implicit none
      real(8),intent(in):: E0 ! lab E0 in Gev
      real(8):: ans   ! xthu
      ans = 12.0*(E0/1.e8)**0.1
      end      function cxthu

      function ck1u( E0 ) result(ans)
      implicit none
      real(8),intent(in):: E0 ! lab E0 in Gev
      real(8):: ans   !  k1u
      ans = 0.3*(E0/1.e8)**0.09
      end      function ck1u

      function cgetParmN( io,  vname, vv ) result(ans)
!          get parameter variable name and given value(s)  from io
       implicit none
       integer io
       character*(*)  vname, vv  ! output
       logical ans

       integer linel
       parameter( linel = 100)
       character*(linel)  line
       integer loc, loc2
       vname = " "
       do while(.true.)
          read(io, '(a)', end=100 ) line
          if( line(1:1) .eq. " " .and. line(2:2) .ne. " ") then
             loc = index( line(3:linel), " ")  + 2
             vname = line(2:loc-1)
             loc2 = index( line, "/")
             if(loc2 .eq. 0 ) then
                write(0,* ) ' "/"  is missing in the input data file '
                write(0,*)  ' The line is: ', line
                stop 1234
             endif
             vv = line(loc+1:linel)  !  some data containes '/' such as '../../Media' so put all
                                     ! data.
             goto 50
          endif
       enddo
 50    continue
       ans = .true.
       return
 100   continue
       ans =.false.
       end function cgetParmN

      function csoftenFixXth(E0)  result(xth)
      implicit none
      real(8),intent(in)::E0 ! proton/pi/K incident E. in GeV
      real(8):: xth !
!        fix the xth above which softening is performed
!       at 10^12 eV:  0.1 
!          10^13      0.063
!          10^14      0.04
!          10^16     0.01585
!          10^17      0.01
!          10^19      0.004
      xth = 0.1d0/(E0/1000.d0)**0.2
      end  function csoftenFixXth

      end module SoftenPiK
