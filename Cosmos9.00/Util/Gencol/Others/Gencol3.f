#include "Zintmodel.h"
#include "ZcosmosBD.h"
      implicit none
#include  "Zcode.h"
#include  "Zptcl.h"
#include  "Ztrackp.h"
      include "Zprivate.h"

      integer i, nev, j, ntp, eof, ntpout
      type (ptcl):: w(maxn)
!      write(*,*) ' this is dummy'
      integer,save:: Jdecay(klast)
!///////////
!      real(8),save:: tau=0.3e-10
      real(8),save:: tau=0.95e-10  ! k0s decay
      Jdecay(:) = 0
!        make Jdecay non zero to see if it is needed to decay,if
!        the  life is < tau.
!          data jcode/kbomega, kdmes, keta, kgzai,  klambda, klambdac,
!     *              kpion, ksigma, kkaon/
!        kpion--> pi0  (pi+/- is never decayed)
!        kkaon--> k0s  (only k0s is examined)
      Jdecay(kbomega) = 1
      Jdecay(kdmes) = 1
      Jdecay(kgzai) = 1
      Jdecay(klambda) = 1
      Jdecay(klambdac) = 1
      Jdecay(ksigma) = 1
      Jdecay(kkaon) = 1
      Jdecay(keta) = 1
      Jdecay(krho) = 1
      Jdecay(kphi) = 1
      Jdecay(komega) = 1
      Jdecay(kpion) = 1
!////////////
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
!//////////
         call cmydecay(Jdecay, tau,  w, ntp, ntpout)
         ntp = ntpout
!////////////
         call cutbyangle(w, ntp, ntp)
         call sortbyke(w, ntp)  ! sort by kinetic energy 
         if(Trace .gt. 0) then
            call outtrace(j, w, ntp)
         endif
         call outresul(w, ntp)
      enddo
      write(0,*) 
     *  ' number of events generated is ',nevent
 100  continue
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
      call cqUhookr(3, trackl)
      call cqUhooki(1, outzero)
!
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
     *        pjcode, pjsub, pjchg, pjpx, pjpy, pjpz
         call cqUhookc(2, input) 
         read(input, *) 
     *        tgcode, tgsub, tgchg, tgpx, tgpy, tgpz
         call cqUhookc(3, input)
         if(input .ne. ' ') then
            read(input, *) xpos, ypos, zpos
            xyz = 1
         else
            xyz = 0
         endif
      endif

      call formpjtg(1)

      call cfixPrefix('configDummy')
      call csetCosOrEpi('epics')
      if( index( IntModel,'qgsjet1') .ne. 0 ) then
#ifdef QGSJET1
         call qgs01init
         ActiveMdl = 'qgsjet1'
#else
         write(0,*) 'to use qgsjet1,  define it  in Zintmodel%h'
#endif
      elseif(index (IntModel, 'sibyll') .ne. 0 )  then
#ifdef  SIBYLL
         call sibyllinit
         ActiveMdl = 'sibyll'
#else
         write(0,*) 'to use sibyll, define it in Zintmodel%h'
#endif
      else
         call cintModels('epics')
         call cfixModel( plab )
      endif

      write(0, *) 'Active int. model=',ActiveMdl
      write(0, *) ' equiv. lab E=', plab%fm%p(4)
      if(xyz .eq. 0) then
!         write(*, '(a)') '#  mulsubKEdir  user '
      else
!         write(*, '(a)') '#  mulsubKExyzdir  user '
      endif
!      write(*, '(a)') '#--------------------------------'
      end
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
         read(inpfileno,*, end=100) xpos, ypos, zpos
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

      integer confirm  !  input. if 0, root s is  not printed.
                      !         else  printed
      real*8   roots, s
!         form projectile  and target

      call cmkptc(pjcode, pjsub, pjchg,   pj)
         
      if(pjcode .ne. kgnuc) then
         pjmnum = 1
      else
         pjmnum = pjsub
      endif
      pj%fm%p(1) = pjpx*pjmnum
      pj%fm%p(2) = pjpy*pjmnum
      pj%fm%p(3) = pjpz*pjmnum
      pj%fm%p(4) =
     *     sqrt(pj%fm%p(1)**2 + pj%fm%p(2)**2 + pj%fm%p(3)**2
     *     + pj%mass**2)
         
!       make taget (rest of moving -z)
      call cmkptc(tgcode, tgsub, tgchg, tg)
      if(tgcode .ne. kgnuc) then
         tgmnum = 1
      else
         tgmnum = tgsub
      endif
      tg%fm%p(1) = tgpx*tgmnum
      tg%fm%p(2) = tgpy*tgmnum
      tg%fm%p(3) = tgpz*tgmnum
      tg%fm%p(4) =
     *     sqrt(tg%fm%p(1)**2 + tg%fm%p(2)**2 + tg%fm%p(3)**2 
     *   +   tg%mass**2)
!       
      if(tgpx .eq. 0. .and. tgpy .eq. 0. .and.
     *        tgpz .eq. 0.)  then
!     target is at rest; 
         s= 2*pj%fm%p(4)*tg%mass +tg%mass**2 + pj%mass**2
      else
!         by  general formula
         s = (pj%fm%p(4)+tg%fm%p(4))**2 -
     *    (pjpx+tgpx)**2 + (pjpy+tgpy)**2+(pjpz+tgpz)**2
      endif
      roots = sqrt(s)
      if(confirm .ne. 0) then
         write(0, *) ' roots/2=', roots/2
      endif
!c           boost to target rest system
      call cbst1(1, tg, pj, plab)
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
!////////
!      call cgetcm(tg, plab, gc, icon)
!      write(0,*) ' gcbeta=', gc
!//////////////////

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
!/////////////
!         call clorez(gc, a(i), pl)
!         write(*,'(2i4, 1p, g13.4, 0p, f10.4, 1p,6g13.4)')
!     *   a(i).code,  a(i).charge,  pt,  eta, pl.fm.p(4)-pl.mass,
!     *     a(i).fm.p(1:4)
!//////////////////
         write(*,'(2i4, 1p, g13.4, 0p, f10.4,1p, g13.4)')
     *   a(i)%code,  a(i)%charge,  pt,  eta, a(i)%fm%p(4)
      enddo
      if(ntp .gt. 0 .or. outzero .eq. 0) then
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
      if( tg%code .eq. knuc ) then
         tA = 1
      elseif( tg%code .eq. kgnuc ) then
         tA = tg%subcode
      else
         write(0,*) ' target code=', tg%code, 'invalid'
         stop 9999
      endif
      tZ =  tg%charge

#ifdef QGSJET1
      if(ActiveMdl .eq. 'qgsjet1') then
         call qgs01event(plab, tA, tZ, a, ntp)
      endif
#elif defined SIBYLL
      if(ActiveMdl .eq. 'sibyll') then
         call sibyllevent(plab, tA, tZ, a, ntp)
      endif
#else
      if(ActiveMdl .eq. 'qgsjet2' ) then
         call cxsecQGS(plab, tA, xs)
      endif
      call chAcol(plab, tA, tZ,  a, ntp)
#endif
      do j = 1, ntp
!               boost to  target mooving system
         call cibst1(j, tg,  a(j), a(j))
      enddo
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
!

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
      real*8  p, wx, wy, wz 
      real  x1, y1, z1, x2, y2, z2
      character*100  tracefile

      write(tracefile, *) TraceDir(1:klena(TraceDir))//'/trace', nev
      call kseblk(tracefile, ' ', leng)
      call copenfw(TraceDev, tracefile(1:leng), icon)
      if(icon .ne. 0) then
         call cerrorMsg('tracefile could not be opened',0)
      endif
      do j = 1, ntp
         i = indx(j)
         p= sqrt( a(i)%fm%p(1)**2 + a(i)%fm%p(2)**2
     *        +      a(i)%fm%p(3)**2 )
         wx = a(i)%fm%p(1)/p               
         wy = a(i)%fm%p(2)/p               
         wz = a(i)%fm%p(3)/p 
         if(xyz .eq. 0) then
            x1 = 0.
            y1 = 0.
            z1 = 0.
         else
            x1 = xpos
            y1 = ypos
            z1 = zpos
         endif
         x2 = x1 + wx*trackl
         y2 = y1 + wy*trackl
         z2 = z1 + wz*trackl
         write(TraceDev,'(3g14.5, i3, g14.4, i3, i2)') 
     *      x1, y1, z1,
     *      a(i)%code,  a(i)%fm%p(4) - a(i)%mass, a(i)%charge,
     *      0
         write(TraceDev, '(3g14.5, i3, g14.4, i3, g14.4)' )
     *      x2, y2, z2,
     *      a(i)%code,  a(i)%fm%p(4) - a(i)%mass, a(i)%charge,
     *      trackl 
         write(TraceDev, *) 
         write(TraceDev, *) 
      enddo
      close(TraceDev)
      end

