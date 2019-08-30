#include "Zintmodel.h"

#ifdef MODIFYX
#include "csoftenPiK.f"      
#endif

#include "ZcosmosBD.h"
#ifdef MODIFYX
      use SoftenPiK
#endif
      implicit none
#include  "Zptcl.h"
#include  "Ztrackp.h"
#include  "Zevhnp.h"

      include "Zprivate.h"

      integer i, nev, j, ntp, eof, ntpo
      type (ptcl):: w(maxn)

      call init
      do j = 1, nevent
         if(inpfileno .gt. 0) then
            call readinpfile(eof)
            if(eof .eq. 1) then
               write(0,*) 
     *        ' number of events generated is ',j-1
               goto 100
            endif
            call formpjtg(0)
         endif
         call gencol(w, ntp)
#ifdef  MODIFYX
!                 in cms we work
         call csoftenFE(pj, fwbw,  w, ntp, ntpo)
         ntp = ntpo
#endif
         call cutbyangle(w, ntp, ntpo)
         ntp = ntpo
         call sortbyke(w, ntp)  ! sort by kinetic energy 
         if(Trace .gt. 0) then
            call outtrace(j, w, ntp)
         endif
         if(outresul1 == 0 ) then
            call outresul(w, ntp)
         elseif(outresul1 == 1 ) then
            call outresulB(w, ntp)
         endif
      enddo
      write(0,*) 
     *  ' number of events generated is ',nevent
 100  continue
      write(0,*)
     * "Equivalent Lab. Energy was ", LabEquivE," GeV/n"
      end

      subroutine init
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zmass.h"
#include  "Zmanager.h"
#include  "Zmanagerp.h"
#include  "Ztrackp.h"

      include  "Zprivate.h"
      character*200 input, file
      character*20 uid
      integer klena, icon, eof

      external cblkManager
      external cblkEvhnp

      call creadParam(5)
      

      if(TraceDir .eq. ' ') then
         call cgetLoginN(uid)
         TraceDir = '/tmp/'//uid(1:klena(uid))
      endif

      if(DestEventNo(2) .eq. 0) then
         nevent =abs( DestEventNo(1) )
      else
         nevent = abs( DestEventNo(2) )
      endif
      call cmkSeed(InitRN(1), InitRN)
      call rnd1i(InitRN)        ! random number init.
      call cqUhookr(1, wzmin)
      call cqUhookr(2, wzmax)
      write(0,*) ' cos* min  max=', wzmin, wzmax
      call cqUhookr(3, trackl)
      call cqUhooki(1, outzero)
      call cqUhooki(2, outresul1)
!       make projectile going +z
      call cqUhookc(1, input) 
      if(input(1:4) .eq. "file") then
         read(input(5:10), *) inpfileno
         xyz=index(input, "xyz")
         call cqUhookc(2, input) 
         file = ' '
         file=input(1:klena(input))
         call copenfw2(inpfileno, file, 1, icon)
         if(icon .ne. 1) then
            write(0,*)
     *      ' file=', file, ' cannot be opened in Gencol'
            stop 9999
         endif
         call readinpfile(eof)
!          once rewind to read successively for event generation
         rewind inpfileno
      else
         inpfileno=0
         read(input, *) 
     *        pjcode, pjsub, pjchg, pjpx, pjpy, pjpz  !  proj. mom/n
         call cqUhookc(2, input) 
         input = trim(input)//" /"
         read(input, *) 
     *        tgcode, tgsub, tgchg, tgpx, tgpy, tgpz, targetName  ! target. mom/n
         call cqUhookc(3, input)
         if(input .ne. ' ') then
!            read(input, *) xpos, ypos, zpos
            read(input, *) xpos(1:3)
            xyz = 1
         else
            xyz = 0
         endif
      endif
      call formpjtg(1)    ! form proj. and target
      call cfixPrefix('configDummy')
      call csetCosOrEpi('gencol') ! this is probably not used


      if( index( IntModel,'qgsjet1') .ne. 0 ) then
#ifdef QGSJET1
         call qgs01init
         ActiveMdl = 'qgsjet1'
#else
         write(0,*) 'to use qgsjet1,  define it  in Zintmodel.h'
#endif
      elseif(index (IntModel, 'sibyll') .ne. 0 )  then
!cc  #ifdef  SIBYLL
         call csibyllinit
         ActiveMdl = 'sibyll'
!c #else
!cc        write(0,*) 'to use sibyll, define it in Zintmodel.h'
!cc #endif
      else
         call cintModels('gencol')
         call cfixModel( plab )
      endif

      write(0, *) 'Active int. model=',ActiveMdl

      if( plab%code == kgnuc ) then
         LabEquivE= plab%fm%p(4)/plab%subcode
      else
         LabEquivE=plab%fm%p(4)
      endif
      if( ActiveMdl == "dpmjet3" ) then
         call checkEnergy
         call checkPtcls
      endif
!             for muon N.I, we must read  mutab for which we must read 
!         all  elemag table too.
      if(pj%code == kmuon ) then
!c         call epReadMediaForMuon  !////////////
      endif
      end
      subroutine  epReadMediaForMuon
#include  "Zcode.h"
#include  "Zptcl.h"
      
      include  "Zprivate.h"
      integer icon

      if( targetName /= " " ) then
         MediaNo = 1
!ccc         call eprdMFile(targetName, icon)  !!////////////
         if( icon /= 0 ) then
            write(0,*) " targetName =", targetName
            write(0,*) " not acceptable "
            stop
         endif
      else
         write(0,*) " projectile code =", plab%code, "= muon"
         write(0,*) "   In this case, media name of the target must"
         write(0,*) " follow the target partcile momentum, e.g,"
         write(0,*) " UserHookc = '1  -1 0   0. 0     500',"
         write(0,*) "             '9  16 8   0  0   0. LiqO2',"
         write(0,*) "             "
         stop
      endif                       
      end


      subroutine checkEnergy
      implicit none
#include "Zptcl.h"
#include "Zprivate.h"
#include "Zcode.h"
      integer icon

      character(80):: input
      character(8):: bb(3)
      integer:: nr
      real(8):: Einput

      if( pj%code == kphoton ) then
         if( tg%code == kgnuc ) then
            if( pj%fm%p(4) < 7. ) then
               write(0,*)
     *         'photon E must be >= 7 for A target'
               stop
            endif
         else
            if( pj%fm%p(4) < 6. ) then
               write(0,*)
     *         'photon E must be >= 6 for non-A target'
               stop
            endif
         endif
      endif
         
      call copenf(11, "dpmjet.inp", icon)
      if(icon /= 0 ) then
         write(0,*) ' dpmjet.inp  cannot be opened'
         stop
      endif
      do while(.true.)
         read(11,'(a)') input
         call ksplit(input,8, 3, bb, nr)
         if( bb(1) == "ENERGY") then
            read(bb(2),*) Einput
            write(0,*) "ENERGY in dpmjet.inp =",  Einput
            if( Einput < LabEquivE) then
               write(0,*) " is  too small; give a value" 
               write(0,*) "little bit larger than ",LabEquivE
               write(0,*) "   if taget is A, ~1.5 % "
               write(0,*) "   if taget is p/n, very little"
                  write(0,*) 
     *           "E.g ",LabEquivE*1.001, " with ~ 4 digit accuracy" 
                  stop
!            elseif( Einput > LabEquivE*1.002 ) then
            elseif(  tg%code == kgnuc .and. 
     *              Einput > LabEquivE*1.040 ) then
               write(0,*) " is  too large; give a vlaue" 
               write(0,*) "little bit larger than ",LabEquivE
               write(0,*) 
     *           "E%g ",LabEquivE*1.03, " with ~ 4 digit accuracy" 
               stop
            elseif( tg%code /= kgnuc .and.
     *             Einput > LabEquivE*1.002 )  then
               write(0,*) " is  too large; give a vlaue" 
               write(0,*) "little bit larger than ",LabEquivE
               write(0,*) 
     *           "E%g ",LabEquivE*1.001, " with ~ 4 digit accuracy" 
               stop
            else
               write(0,*) " is close to ", LabEquivE, " and OK"
            endif
            close(11)
            exit
         endif
      enddo
      end

      subroutine checkPtcls
!       check projectile and target  ;
!       Are they same in param and dpmjet.inp ? 
      implicit none
#include "Zptcl.h"
#include "Zprivate.h"
#include "Zcode.h"
      integer icon

      character(80):: input
      character(8):: bb(5)
      integer:: nr
      character(8)::btype
      integer code, subc, charge

      call copenf(11, "dpmjet.inp", icon)
      if(icon /= 0 ) then
         write(0,*) ' dpmjet.inp  cannot be opened'
         stop
      endif
      do while(.true.)
         read(11,'(a)', end=100) input
         call ksplit(input,8, 5, bb, nr)

         if( bb(1) == "PROJPAR") then
            if( nr == 2 ) then
               btype=bb(2)
               call cbtype2cos(btype, code, subc, charge)
               if( code /= pj%code  .or. charge /= pj%charge ) then
                     ! subcode not checked 
                  write(0,*) 
     *                 ' projectile mismatch bw dpmjet.inp and param'
                  write(0,*) ' dpmjet.inp :', btype
                  write(0,*) ' param: ',  pj%code, pj%charge
                  stop
               endif
            elseif( nr == 3 ) then
               read( bb(2), *) subc
               read( bb(3), *) charge
               if( subc > 1 ) then
                  if( subc /= pj%subcode .or. charge /= pj%charge ) then
                     write(0,*) 'proj  mismatch bw dpmjet.inp and param'
                     stop
                  endif
               else
                  if( pj%code /= knuc .and. charge /= pj%charge ) then
                     write(0,*) 'proj  mismatch bw dpmjet.inp and param'
                     stop
                  endif
               endif
            endif
         elseif(bb(1) == "TARGETPAR") then
            if( nr == 2 ) then
               btype=bb(2)
               call cbtype2cos(btype, code, subc, charge)
               if( code /= tg%code  .or. charge /= tg%charge ) then
                     ! subcode not checked 
                  write(0,*) 
     *         ' target mismatch bw dpmjet.inp and param'
                  write(0,*) ' dpmjet.inp :', btype
                  write(0,*) ' param: ',  tg%code, tg%charge
                  stop
               endif
            elseif( nr == 3 ) then
               read( bb(2), *) subc
               read( bb(3), *) charge
               if( subc > 1 ) then
                  if( subc /= tg%subcode .or. charge /= tg%charge ) then
                     write(0,*)
     *                  'target mismatch bw dpmjet.inp and param'
                     stop
                  endif
               else
                  if( tg%code /= knuc .and. charge /= tg%charge ) then
                     write(0,*) 
     *                  'target mismatch bw dpmjet.inp and param'
                     stop
                  endif
               endif
            endif 
         endif
      enddo
 100  continue
      close(11)
      write(0,*)
     *  ' proj targ consistency btw dpmjet.inp and param: OK'
      end
      subroutine  cbtype2cos(btype, code, subc, charge)
      implicit none
#include "Zcode.h"
      character(8),intent(in):: btype
      integer,intent(out):: code, subc, charge
      
      select case(btype)
         case("PROTON") 
           code = knuc
           charge = 1
           subc = -1
        case("APROTON") 
           code = knuc
           charge = -1
           subc = 1
        case("PHOTON")
           code = kphoton
           charge = 0
           subc = 0
        case("NEUTRON")
           code = knuc
           charge = 0
           subc = -1
        case("ANEUTRON") 
           code = knuc
           charge = 0
           subc = 1
        case("PION+")
           code = kpion
           charge = 1
           subc = -1
        case("PION-")
           code = kpion
           charge = -1
           subc = 1
        case("KAON+")
           code = kkaon
           charge = 1
           subc = -1
        case("KAON-")
           code = kkaon
           charge = -1
           subc = 1
        case("LAMBDA")
           code = klambda
           charge = 0
           subc = -1
        case("ALAMBDA")
           code = klambda
           charge = 0
           subc = 1
        case("PIZERO")
           code = kpion
           charge = 0
           subc = 0
        case default
           write(0,*) ' not acceptable ptcl=',btype
           stop
        end select
!     'KAONSHRT' ,
!     &'SIGMA-  ' , 'SIGMA+  ' , 'SIGMAZER' , 'PIZERO  ' , 'KAONZERO' ,
!     &'AKAONZER' 
      end subroutine cbtype2cos

      subroutine readinpfile(eof)
      implicit none
#include "Zptcl.h"      
      include "Zprivate.h"

      integer eof ! output . data read--> 0
                  !   eof reached --> 1
      read(inpfileno,*, end=100)
     *     pjcode, pjsub, pjchg, pjpx, pjpy, pjpz
      read(inpfileno,*, end=100)
     *     tgcode, tgsub, tgchg, tgpx, tgpy, tgpz
      if(xyz .gt. 0 ) then
!         read(inpfileno,*, end=100) xpos, ypos, zpos
         read(inpfileno,*, end=100) xpos(1:3)
      endif
      eof = 0
      return
 100  continue
      eof = 1
      end
!     *******************
      subroutine formpjtg(confirm)
!     ******************
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zmass.h"
#include  "Zmanager.h"
#include  "Zmanagerp.h"
#include  "Ztrackp.h"

      include "Zprivate.h"
      type (ptcl)::resttg, Zcmsp, Zplab
      integer confirm  !  input. if 0, root s is  not printed.
                      !         else  printed
      real*8   roots, s
!         form projectile  and target
      integer::icon

      call cmkptc(pjcode, pjsub, pjchg,   pj)
!       pjmnum:    proj. mass number in integer
      if(pjcode .ne. kgnuc) then
         pjmnum = 1
      else
         pjmnum = pjsub
      endif
      pj%fm%p(1) = pjpx*pjmnum     ! total mom.
      pj%fm%p(2) = pjpy*pjmnum
      pj%fm%p(3) = pjpz*pjmnum
      pj%fm%p(4) =
     *     sqrt(pj%fm%p(1)**2 + pj%fm%p(2)**2 + pj%fm%p(3)**2
     *     + pj%mass**2)
         
!       make taget (rest or moving -z or ...)
      call cmkptc(tgcode, tgsub, tgchg, tg)
      if(tgcode .ne. kgnuc) then
         tgmnum = 1
      else
         tgmnum = tgsub
      endif
      tg%fm%p(1) = tgpx*tgmnum   ! total mom.
      tg%fm%p(2) = tgpy*tgmnum
      tg%fm%p(3) = tgpz*tgmnum
      tg%fm%p(4) =
     *     sqrt(tg%fm%p(1)**2 + tg%fm%p(2)**2 + tg%fm%p(3)**2 
     *   +   tg%mass**2)
!       
      if(tgpx .eq. 0. .and. tgpy .eq. 0. .and.
     *        tgpz .eq. 0.)  then
!     target is at rest;  s = (Ep+Et)^2 - (Pp+Pt)^2
!                           = (Ep+Mt)^2 - Pp^2
!                           =  Mp^2 +2EpMt +Mt^2
!   
         s= 2*pj%fm%p(4)*tg%mass +tg%mass**2 + pj%mass**2
      else
!         by  general formula;
!               Mp^2 + Mt^2 +2(Ep*Et - Pp.Pt)
         s = pj%mass**2 + tg%mass**2 +
     *       2*(pj%fm%p(4)*tg%fm%p(4) -
     *         dot_product(pj%fm%p(1:3), tg%fm%p(1:3) ) )
      endif
      roots = sqrt(s)
!      if(confirm .ne. 0) then
         write(0, *) ' roots/2=', roots/2
         write(0, *) 's,roots above are total not by /n'
!      endif
!c           boost to target rest system
      call cbst1(1, tg, pj, plab)
!         make rest target
      call cbst1(2, tg, tg, resttg)
       !   get /n cms
      call cgeqm2(plab, resttg, Cmsp, icon)   
!////////////
      write(0,*) ' Cmsp= (/n)',Cmsp%fm%p(:), Cmsp%mass
      write(0,*) ' plab=',plab%fm%p(4)
      write(0,*) ' tg=', tg%fm%p(4)
!///////////
      if(icon /= 0 ) then
         write(0,*) ' input 2 ptcls cannot form CMS'
         stop
      endif

!       Next part  will be used only by EPOS 
!       get equiv. CMS particle (Zcmsp) when projectile Pt is 0.
!       If crossing angle is /= 0, Cmsp and Zcmsp  
!       differ and we must use Zcmsp here!!!!.   2016/01/24
      zplab = plab
      zplab%fm%p(1:2) = 0.
      zplab%fm%p(3) = sqrt( zplab%fm%p(4)**2- zplab%mass**2)
       !   get /n cms
      call cgeqm2(zplab, resttg, Zcmsp, icon)   
!         inform Zcmsp to cintModels.f 
!      call cputGencolCMS(Cmsp) 
      call cputGencolCMS(Zcmsp) 
      end
!     ************************
      subroutine outresul(a, ntp)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zmass.h"
#include  "Zmanagerp.h"
#include  "Zmanager.h"
#include  "Ztrackp.h"
      include  "Zprivate.h"

      integer ntp
      type (ptcl):: a(ntp)
      integer  i, j
      real*8  p
      real*8  wx(3)
      integer nw, difcode(20)

      call getDiffCode(nw, difcode)

      do j = 1, ntp
         i = indx(j)
!         p= sqrt( a(i).fm.p(1)**2 + a(i).fm.p(2)**2
!     *        +      a(i).fm.p(3)**2 )
         p = sqrt(
     *  dot_product( a(i)%fm%p(1:3), a(i)%fm%p(1:3)) )
         wx(1:3)= a(i)%fm%p(1:3)/p
!         wx = a(i).fm.p(1)/p               
!         wy = a(i).fm.p(2)/p               
!         wz = a(i).fm.p(3)/p               
         if(xyz .eq. 0) then
            write(*,
     *     '(i3,i5,i4,1p, g14.5,0p, 3f17.13, i8)')
     *        a(i)%code, a(i)%subcode, a(i)%charge,
     *        a(i)%fm%p(4)-a(i)%mass, wx(1:3), j

!            write(*,
!     *     '(i3,i5,i4,1p, g14.5,0p, 3f17.13, 1p,3g12.4, i8)')
!     *        a(i).code, a(i).subcode, a(i).charge,
!     *        a(i).fm.p(4)-a(i).mass, wx(1:3), 
!     *        a(i).fm.p(1:3),   j

         else
            write(*,'(i3,i5,i4,g14.5,1p3E11.3,0p3f17.13, i8)')
     *        a(i)%code, a(i)%subcode, a(i)%charge,
!     *        a(i).fm.p(4)-a(i).mass, xpos, ypos, zpos, 
     *        a(i)%fm%p(4)-a(i)%mass, xpos(1:3), wx(1:3), j
!     *        wx, wy, wz, j
         endif
      enddo
      if(ntp .gt. 0 .or. outzero .eq. 1) then
         write(*, *) 
      endif
      end

      subroutine  gencol(a, ntp)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnv.h"
#include  "Zevhnp.h"
#include  "Zmanagerp.h"
      include "Zprivate.h"

      type (ptcl)::  a(*)
!             projectile and target information (both befor
!             and after collision ) in different system.
!
      integer  ntp
      integer j
      integer tZ, tA
      real*8  xs
!     
      tA = tgmnum
      tZ =  tg%charge
      if(pj%code == kmuon ) then
!         sepcial treatment
         call epGencolByMuon(plab, tA, tZ,  a, ntp)
!!      elseif( pj.code == kphoton ) then
!         sepcial treatment
!!         call epGencolByPhoton(plab, tA, tZ,  a, ntp)
      elseif(pj%code == kelec .or.
     *      ( pj%code /= kgnuc .and. pj%code >knuc ) ) then
         write(0,*) ' ptcl code =', pj%code, ' not '
         write(0,*) ' supported in Gencol'
         stop
      else
!           init for 1 event 
         if(ActiveMdl .eq. 'qgsjet2' ) then
            call cxsecQGS(plab, tA, xs)
         elseif(ActiveMdl .eq. 'epos' ) then
            call ceposIniOneEvent(plab, tg, xs)
         endif

!           generate event
         if(ActiveMdl .eq. 'qgsjet1') then
#ifdef QGSJET1
            call qgs01event(plab, tA, tZ, a, ntp)
#endif
         elseif(ActiveMdl .eq. 'sibyll') then
!cc #ifdef SIBYLL
            call csibyllevent(plab, tA, tZ, a, ntp)
!cc #endif
         else
            call chAcol(plab, tA, tZ,  a, ntp)
         endif
      endif
      call getImpactParam(b)
      do j = 1, ntp
!               boost to original target system
         call cibst1(j, tg,  a(j), a(j))
      enddo
      end
      subroutine  epGencolByMuon
      implicit none
!      call epInfoPhotoP(incGp);  incGpをHowPhotopに焼き直し
!      call ep2cosCond
!      call cfixModel( cTrack.p )
!      call ciniSmpIntL   ! not related.
!      call epsmpNEPIntL(Media(MediaNo))  !st decay length only
!      call ep2cosCondr  !set FromEpics = .false.
!  if(   cTrack.p.fm.p(4) .gt.  Media(MediaNo).cnst.muNEmina )
!      call epmuNsmpP( Media(MediaNo),
!     *    cTrack.p.fm.p(4), prob, path)
!      epfixProc    !  ProcessNo 
!
!      epmuinte:
!       call epmuNsmpE(Media(MediaNo), Move.Track.p.fm.p(4),
!     *                   Et)
!      cTrack.p.fm.p(4) = cTrack.p.fm.p(4) - Et -->muon
!        if(Et .gt. 152.d-3) then
!            if(Media(MediaNo).mu.MuNI .eq. 3 ) then
!c             generate gamma-N interaction; employ gamma interaction  
!c             routine                                                 
!               cTrack.p.fm.p(4) = Et
!               call cmkptc( kphoton, 0, 0, cTrack.p)
!               call epe2p(cTrack) ! adjust momentum                   
!c                                                                     
!               call ep2cosPtcl( cTrack.p )
!c                 for small basic cross section case.                 
!               call epfixTarget2(ActiveMdl, Media(MediaNo))
!               call ep2cosCond2(Media(MediaNo).colA,
!     *                          Media(MediaNo).colZ)
!               call cphotop     ! Cosmos function     
!         call eppushPtcl(cTrack) ! use pos. info from this ptcl

      end

      subroutine  epGencolByPhoton(plab, tA, tZ,  a, ntp)
      implicit none
!  #include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zelemagp.h"

!       include "Zprivate.h"

      type (ptcl):: plab  ! input proj. in Lab
      integer,intent(in):: tA   ! target A
      integer,intent(in):: tZ   ! target Z
      type (ptcl):: a(*)        ! ouptut. produced ptcl
      integer,intent(out):: ntp ! # of produced ptcls

      real(8):: TargetA, xs
      if( HowPhotoP > 0 .and.  plab%fm%p(4) > MinPhotoProdE ) then
!         call ep2cosCond
!         call cfixModel( pj )
!         call cgpXsec(Media(MediaNo).A,  E, xs)
!         prob = xs*Media(MediaNo).mbtoPX0
!         call rndc(u)
!         tgp = -log(u)/prob
!　　　−ーーーーーーーーーーー above is probably not needed

!////         call ep2cosPtcl(plab)
!         call epfixTarget2(ActiveMdl, Media(MediaNo))
!         call ep2cosCond2(Media(MediaNo).colA, 
!       *         Media(MediaNo).colZ, Media(MediaNo).colXs)
         TargetA= tA
         call cgpXsec(TargetA,  plab%fm%p(4), xs)
!////         call ep2cosCond2(tA, tZ, xs)   ! xs will not be used
         call cphotop        ! Cosmos function                     
!!         call eppushPtcl(cTrack)  ! puch Cosmos made  ptcl into Epics
         a(1:Nproduced) = Pwork(1:Nproduced)
         ntp = Nproduced
      else
         write(0,*) ' HowPhotoP=', HowPhotoP, ' Eg=', plab%fm%p(4)
         write(0,*) ' Either of above NG for photon projectile'
         stop
      endif

      end

      subroutine cutbyangle(a, ntp0,  ntp)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnv.h"
#include  "Zevhnp.h"
#include  "Zmanagerp.h"
      include "Zprivate.h"
      type (ptcl)::  a(*)
      integer ntp0 ! input. number of ptcls. in a
      integer ntp  ! output. could be the same as ntp0
      integer j 
      integer  i
      real*8 p, wz
      j = 0
      do i = 1, ntp0
         p = a(i)%fm%p(1)**2 + a(i)%fm%p(2)**2 +
     *       a(i)%fm%p(3)**2
         p = sqrt(p)
         wz = a(i)%fm%p(3)/p 
         if( wz .ge. wzmin .and. wz .le. wzmax ) then
            j = j + 1
            a(j)=a(i)
         endif
      enddo
      ntp = j
      end
      subroutine sortbyke(a, ntp)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"

      include "Zprivate.h"
      integer  ntp
      type (ptcl)::  a(*)
!             projectile and target information (both befor
!             and after collision ) in different system.
      integer  i
      do i = 1, ntp
         ke(i) = a(i)%fm%p(4) - a(i)%mass
      enddo
      call kqsortd(ke, indx, ntp)
      call ksortinv(indx, ntp)  
!       ke( indx(1) ) is the highest energy
      end
      subroutine outtrace(nev, a, ntp)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zmass.h"
#include  "Zmanagerp.h"
#include  "Zmanager.h"
#include  "Ztrackp.h"
      include  "Zprivate.h"

      integer ntp, nev
      type (ptcl):: a(ntp)
      integer  i, j, leng, icon, klena
      real*8  p
!     , wx, wy, wz 
!      real  x1, y1, z1, x2, y2, z2
      character*100  tracefile
      real(8):: x1(3), x2(3)

      write(tracefile, *) TraceDir(1:klena(TraceDir))//'/trace', nev
      call kseblk(tracefile, ' ', leng)
      call copenfw(TraceDev, tracefile(1:leng), icon)
      if(icon .ne. 0) then
         call cerrorMsg('tracefile could not be opened',0)
      endif
      if(xyz .eq. 0) then
         x1(1:3)=0.
      else
         x1(1:3)=xpos(1:3)
      endif

!          colliding partilces
      if(Trace == 1 ) then
         call puttrace(x1, pj, -2*trackl)
         call puttrace(x1, tg, -2*trackl) 
      endif

      do j = 1, ntp
         i = indx(j)
         call puttrace(x1, a(i), trackl)
      enddo
      close(TraceDev)
      end
      subroutine puttrace(x1, a, leng)
      implicit none
#include  "Zmanagerp.h"
#include  "Zptcl.h"
#include  "Ztrackp.h"
      include  "Zprivate.h"

      real(8),intent(in)::x1(3)  !   origin
      type (ptcl):: a   ! a partilce
      real(8),intent(in):: leng    ! length of track to be drawn

      real(8)::wx(3), x2(3), p

      p = sqrt( dot_product(a%fm%p(1:3), a%fm%p(1:3)))
      if( p >0. ) then
         wx(1:3)= a%fm%p(1:3)/p
      else
         wx(1:3) = 0.
      endif
      x2(1:3) = x1(1:3) + wx(1:3)*leng
      write(TraceDev,'(3g14.5, i3, g14.4, i3, i2)') 
     *     x1(1:3),
     *     a%code,  a%fm%p(4) - a%mass, a%charge,
     *     0
      write(TraceDev, '(3g14.5, i3, g14.4, i3, g14.4)' )
     *     x2(1:3),
     *     a%code,  a%fm%p(4) - a%mass, a%charge,
     *     trackl 
      write(TraceDev, *) 
      write(TraceDev, *) 
      end

      subroutine outresulB(a, ntp)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zmass.h"
#include  "Zmanagerp.h"
#include  "Zmanager.h"
#include  "Ztrackp.h"
      include  "Zprivate.h"
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

      integer ntp
      type (ptcl):: a(ntp)
      integer  i, j

      real*8  p, wx, wy, wz, pt, y, eta
      integer npzm, npzp, Ncht, Nchpzm, Nchpzp
      integer nw, diffcode(20) 
!//////////////
      type (ptcl):: inlab
      real(8):: ylab, etalab
      integer icon
      type (ptcl):: pcms, pl
      type (fmom):: gc
!/////////////
      npzm = 0
      npzp = 0
      Nchpzm = 0
      Nchpzp = 0

      do j = 1, ntp
         i = indx(j)
!         pt = sqrt( a(i).fm.p(1)**2 + a(i).fm.p(2)**2)
!         p= sqrt( pt**2  +    a(i).fm.p(3)**2 )
         if( a(i)%fm%p(3) .gt. 0.) then
            npzp = npzp + 1
            if(a(i)%charge .ne. 0) then
               Nchpzp = Nchpzp +1
            endif
         else
            npzm = npzm + 1
            if(a(i)%charge .ne. 0) then
               Nchpzm = Nchpzm + 1
            endif
         endif
      enddo
      Ncht = Nchpzm + Nchpzp
      call getDiffCode(nw, diffcode)
      write(*,'("h ",i3, 6i6)' ) 
     *   diffcode(1), ntp, npzm, npzp, Ncht, Nchpzm, Nchpzp

      do j = 1, ntp
         i = indx(j)
         pt = sqrt( a(i)%fm%p(1)**2 + a(i)%fm%p(2)**2)
!         p= sqrt( pt**2  +    a(i).fm.p(3)**2 )
!
!         wx = a(i).fm.p(1)/p               
!         wy = a(i).fm.p(2)/p               
!         wz = a(i).fm.p(3)/p               
!         if(xyz .eq. 0) then
!            write(*,'(3i3,g14.5,3f17.13, i5)')
!     *        a(i).code, a(i).subcode, a(i).charge,
!     *        a(i).fm.p(4)-a(i).mass, wx, wy, wz, j
!         else
!            write(*,'(3i3,g14.5,1p3E11.3,0p3f17.13, i5)')
!     *        a(i).code, a(i).subcode, a(i).charge,
!     *        a(i).fm.p(4)-a(i).mass, xpos, ypos, zpos, 
!     *        wx, wy, wz, j
!         
!         endif
         call cyeta(a(i),  y, eta)
!////////////
!c           boost to target rest system
         call cbst1(1, tg, a(i), inlab)
         call cyeta(inlab,  ylab, etalab)
!///////////////
!/////////////
!         call clorez(gc, a(i), pl)
!         write(*,'(2i4, 1p, g13.4, 0p, f10.4, 1p,6g13.4)')
!     *   a(i).code,  a(i).charge,  pt,  eta, pl.fm.p(4)-pl.mass,
!     *     a(i).fm.p(1:4)
!//////////////////
!         write(*,'(3i4, 1p, g13.4, 0p, f10.4,1p, 6g13.4)')
!     *   a(i)%code,  a(i)%subcode,  a(i)%charge,
!     *     pt, y,  eta, a(i)%fm%p(4)-a(i)%mass,  a(i)%fm%p(3),
!     *     ylab, etalab
         write(*,'(3i4, 1p, 2g13.4)')
     *   a(i)%code,  a(i)%subcode,  a(i)%charge,
     *      a(i)%fm%p(4)-a(i)%mass,  pt

         
      enddo
      if(ntp .gt. 0 .or. outzero .eq. 1) then
         write(*, *) 
      endif
      end
