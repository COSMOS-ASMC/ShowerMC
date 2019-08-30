      subroutine epmuPrsmpP(media, Emu, prob, path)
      implicit none
#include "Zmedia.h"

       type(epmedia):: media    !input. media
      real*8 Emu             ! input.  muon total energy in GeV
      real*8 prob            ! output. muon pair creation prob. /X0
      real*8 path            ! output. sampled path in r.l
                             ! if Emu < media.cnst.muPrEmin,  prob=0 
                             !  and path becomes big

      real*8 u, ale

      if(Emu .le. media%cnst%muPrEmin) then
         prob = 0.
      elseif(Emu .le. media%cnst%muPrEmax1) then
         ale = log10(Emu)
         call kintp3(media%tbl%MuPrTX, 
     *   1,  media%cnst%muPrTXT, media%cnst%muPrLEmin,
     *   media%cnst%muPrdETX, ale, prob)
      else
!          prob is const for  v > vmin; Emu > muPrEmax1         
         prob = media%tbl%MuPrTX(media%cnst%muPrTXT)
      endif
      if(prob .gt. 0.) then
         call rndc(u)
         path =- log(u)/prob
      else
         path = 1.d30
      endif
      end
      
      subroutine epmuPrsmpE(media, Emu, Epair)
      implicit none
#include "Zmedia.h"
!     ****************
      integer PairId(maxErgNodeForMuPair), mupairused
      common /muPairIDcom/ PairId, mupairused
!     **************

       type(epmedia):: media  ! input.  media
      real*8 Emu           ! input. muon total energy in GeV
      real*8 Epair         ! output. sampled energy loss of muon
                           !        (Epair/Emu > media.cnst.muPrVmin)
      integer,save::first=0 

      real*8  a, b, v
      if(Emu <= media%mu%muPairErg(maxErgNodeForMuPair)  .and.
     *   Emu >= media%mu%muPairErg(1) ) then
!         if( media.mu.NoOfErgNodeForMuPair == 0 ) then
         if( first == 0 ) then
!                 init for csampAF
            call epMuPaircsampAF0(media)
!            media.mu.NoOfErgNodeForMuPair =maxErgNodeForMuPair
            mupairused = maxErgNodeForMuPair
            first = 1
         endif
         call epMuPaircsampAF(media, Emu, v)
         Epair = Emu*v
      elseif( Emu > media%mu%muPairErg(maxErgNodeForMuPair) ) then
         a = media%mu%pa*Emu**(-0.40) + media%mu%ra
         b = media%mu%pb*Emu**media%mu%qb  +  0.92

         call epmuPrsmp0(a, b, media%cnst%muPrVmin, v)
         Epair = Emu*v
      else
         Epair = 0.
      endif
      end


!     **********************************
      subroutine epmuPrsmp0(a, b, vc, v)
!     **********************************
      implicit none
!       sample v from  dv/v**(2-b)/(1+av**b)**2
!       This is a good approximation for pair spectrum
!       from muon 
      real*8 a     ! input.
      real*8 b     ! input.
      real*8 vc    ! input.  minimum of v
      real*8 v     ! output. sampled v


      real*8 u, temp


      logical ok

      ok = .false.
      do while (.not. ok)
         if(b .ne. 1.0d0) then
            temp = vc**(b-1)
            call rndc(u)
            v = ((1-temp)*u + temp)**(1./(b-1))
            call rndc(u)
            ok =(u .lt. ((1.+a*temp*vc)/(1.+a*v**b))**2) 
         else
            call rndc(u)
            v = vc**u
            call rndc(u)
            ok =(u .lt. ((1.+a*vc)/(1.+a*v))**2) 
         endif
      enddo

      end

      subroutine  epMuPaircsampAF0(mediain)
!           initialize  pair creation by muon for low energies
!        ( 3 GeV ~ 1.5TeV ),  which is difficult by normal sampling
!        method. We use csampAF so we must make a ds/dv table 
!        at several energies and use that table for general Emu
!        in   epMuPaircsampAF.
!
      use modcsampAF
      implicit none
#include "Zmedia.h"
#include "ZmuBPNgene.h"

       type(epmedia):: mediain    !input

!     ****************
      integer PairId(maxErgNodeForMuPair), mupairused
      common /muPairIDcom/ PairId, mupairused
!     **************
      integer,save:: first=0  ! this is not common

      integer noOfv
      parameter(noOfv=100)

      real(8):: xmin, x, xmax, tprob, f, epmuPrS,epmuPairVmn
      real(8):: v(noOfv), dsdv(noOfv)
      integer i, j

      if( first == 0 ) then
         xmax = 1.
!          we have to inform media to ZmuBPNgene
!          which is used in ds/dv calculation
         media = mediain 
         do i = 1, maxErgNodeForMuPair
            Emu = media%mu%muPairErg(i)
            xmin = epmuPairVmn(Emu)
            xmin = max(xmin, media%cnst%muPrVmin)
            call eptotcmuP(xmin, xmax, tprob)
            j = 0
            x = xmin
            do while ( x< xmax) 
               f = epmuPrS(x)/tprob
               if( f > 0. ) then
                  j = j + 1
                  if(j > noOfv ) then
                     write(0,*) ' too large # of v'
                     write(0,*) ' in epMuPaircsampAF0 '
                     write(0,*) ' increase noOfv there'
                     stop
                  endif
                  v(j) = x
                  dsdv(j) = f/x
               elseif( x > 0.98 ) then
                  exit
               endif
!              1.e-4       1.e-3  3e-3  1.e-2   1.e-1     1.0
!
!             default: 1e-4 to 3e-3: ~60
!                      3e-3 to 0.8:  ~25
!                      0.8  to 1:    ~10
               if( x > 0.8 ) then
                  x = x + 0.02d0
               elseif( x > 3.d-3) then
                  x = x * 10.d0**0.1
               else
                  x = x * 10.d0**0.025
               endif
            enddo
!            iform (v,ds/dv) to csampAF
!            call csampAF0byArray( v, dsdv, j,  muPairId(i))
            call csampAF0byArray( v, dsdv, j,  PairId(i))
         enddo
      endif
      first = 1
!      mediain = media   ! not needed; only ...Id(i) may recive
!       the information.
      end
      subroutine epMuPaircsampAF(media, Emu, v)
      use modcsampAF
      implicit none
#include "Zmedia.h"
!     ****************
      integer PairId(maxErgNodeForMuPair), mupairused
      common /muPairIDcom/ PairId, mupairused
!     **************

       type(epmedia)::  media
      real(8),intent(in):: Emu ! muon energy in GeV
      real(8),intent(out):: v  ! sampled fractonal energy
                               ! of electron pair: Epair/Emu
      real(8)::logE, logE1, logE2
      real(8)::u, probE1, vmin, epmuPairVmn
      integer i


      vmin = epmuPairVmn(Emu)
      vmin = max(vmin, media%cnst%muPrVmin)

!      do i = 1, media.mu.NoOfErgNodeForMuPair-1
      do i = 1, mupairused - 1
         if( Emu >= media%mu%muPairErg(i) .and.
     *       Emu <= media%mu%muPairErg(i+1)  ) then
            logE = log(Emu)
            logE1 = media%mu%logmuPairErg(i) 
            logE2 = media%mu%logmuPairErg(i+1) 
            probE1 = (logE2-logE)/(logE2-logE1)

            v = 0.
            do while( v< vmin)
               call rndc(u)
               if(u <= probE1) then
!                  call csampAF(media.mu.muPairId(i), v)
                  call csampAF(PairId(i), v)
               else
!                  call csampAF(media.mu.muPairId(i+1), v)
                  call csampAF(PairId(i+1), v)
               endif
            enddo
            return  ! ************
         endif
      enddo
      write(0,*) ' Emu=',Emu, ' for epMuPaircsampAF '
      write(0,*) ' strange ;media=',media%name
      stop
      end

               
      
                 
