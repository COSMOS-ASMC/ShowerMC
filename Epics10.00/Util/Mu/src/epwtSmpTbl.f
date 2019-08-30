!
!    For example, if we want to express the atmosphere by a number of
!    spherical shell, we need a lot of 'Air' with different density.
!    Except for small changes due to the LPM effect (at energies < 10^17 eV)
!    many constants are the same and sampling table can be used commonly.
!    So we may define a few Air with different density, and many other 
!    Air with different density can be treated simply changing some of the
!    constants which depend on the density. Such constants are listed
!    below.
!
!    density dependent quantity:
!   
!        media.X0
!        media.wp
!        media.nd
!        media.gtocm
!        media.mbtoPcm
!        media.mbtoPcm2
!        media.sh.c
!        media.sh.x0
!        media.sh.x1
!        media.sh.sa
!        media.sh.xa
!        media.dEdxatp3m
!        cnst.BremEeminLPM
!        cnst.BrEgminH
!        cnst.BrEe1H
!        cnst.BrLEe1H
!        cnst.BrEe2H
!        cnst.PairEgmaxL
!        cnst.PrEg1H
!        cnst.PrLEg1H
!        cnst.PrEg2H
!

!     *******************************
      subroutine epwtmedia(media)
!       print some of basic media information(tbl, cnst are not
!      included)
!        should be called after epGetEffZA, epExpot
!
      implicit none
#include "Zmedia.h"
       type(epmedia):: media


!
!      to avoid that the output by DEC fortran is put "\n" for
!      a long write(*,*), we use  write(*,'(a)')
!

      write(*,'(a,f7.2,a)')
     *   ' Excitation (Ionization) potential =',
     *   media%I*1.d9, ' eV'
      write(*,'(4(a,1p,g13.5))')
     *   ' radiation length=', media%X0, ' cm', media%X0g, ' g/cm^2'

      write(*, '(4(a, 1p, g13.5))')
     *  ' <Aeff>=',  media%Aeff,  ' <Zeff=>',  media%Zeff, 
     *           ' <Z2/A>= ', media%Z2byAeff

      write(*,  '(4(a, 1p, g13.5))')
     *  ' plasma energy=', media%wp,  ' GeV.  No density=',
     *    media%nd, ' /cm^3'

      write(*, '(4(a, 1p, g13.5))')
     * ' <Z2eff=>', media%Z2eff, ' <Z/A>=', media%ZbyAeff
      write(*,'(a,1p,g13.5)') 'Migdal s1=',  media%s1
      write(*,*) ' Conversion factors:'
      write(*,*) ' g/cm^2 to cm=', media%gtocm

      write(*, '(4(a, 1p, g13.5))')
     *  ' mb to /(g/cm2)=', media%mbtoPgrm,
     *  ' mb to /cm=', media%mbtoPcm, ' mb to /X0=', media%mbtoPX0

      write(*,'(a)') ' If a comound is specified as an Atom, use the'
      write(*,'(a,1p,3g13.5)') ' following instead: ', 
     * media%mbtoPgrm2, media%mbtoPcm2, media%mbtoPX02

      write(*,'(a, 1p,g13.5)')
     *  ' pi x Re**2 * N* Z /A *X0g  = 0.15 Z/A*X0g=', media%basearea

      write(*, '(a, 1p,3g13.5)')
     * ' media%A, Z, Z^2=',media%A, media%Z, media%Z2
      write(*, *) 
      write(*, '(a)') ' **** Sternheimers consts'
      write(*, '(a,1p,3g13.5)')
     * ' a, b, c =', media%sh%a, media%sh%b, media%sh%c

      write(*,'(a, 1p,2g13.5)') 
     *    ' x0, x1 =', media%sh%x0, media%sh%x1
      write(*, *) ' sa=', media%sh%sa

      write(*,'(4(a, 1p,g13.5))')
     * ' minimum dEdx_restricted=',media%dEdxatp3m,
     *    ' GeV/(g/cm^2)  with Tcut=',media%sh%tcut, ' GeV'
      write(*, *) 
      write(*, '(a)') ' **** photoelectric consts'
      write(*, '(a,1p,3g13.5)')
     *  'b0, b1, b2=', media%pe%b0, media%pe%b1, media%pe%b2

      write(*, '(a, 1p, 3g13.5)')
     * 'fa, a, p=', media%pe%fa, media%pe%a, media%pe%p
      write(*, '(a, 1p,2g13.5)') 'l, ek=', media%pe%l, media%pe%ek
      write(*, *)
!          separator; when having read basic media later, skip next to
!          this separator (*-format is no usable for compaq, dec)
      write(*, '(a)') '#-#-#-#-#-#-#-#-#-#-#-#-'
      end

!     ******************************* 
      subroutine epwt1dTbl(com, erg,  tbl, size, name)
!           print  total cross section
      implicit none
      integer size
      real*8 tbl(size), erg(size)
      character(*),intent(in):: com, name

      integer i

      write(*, '(a)' ) trim(name)
      write(*, '(a)' ) trim(com)
      do i = 1, size
         write(*,'(1p,2g13.5)') erg(i),  tbl(i)
      enddo
      write(*,*)
      end
!     ******************************************
      subroutine epwt2dTbl(com, tbl, sizeu,  sizee)
!  
      implicit none
      integer sizeu,sizee
      real*8 tbl(sizeu,  sizee)
      character(*),intent(in):: com 

      integer i, iu

      write(*, '(a)') trim(com)

      do i = 1, sizee
         write(*, '(7f11.5)') ( tbl(iu, i), iu = 1, sizeu )
         write(*, *)
      enddo
      write(*,*)
      end
!     ****************************
      subroutine epwtBrCnstS(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst

      write(*, '(a)') 
     * ' constant for lower Seltzer Brems sampling table '
      write(*, '(1p, 5g15.7)')
     *  cnst%BrEeminS, cnst%BrEgminS, cnst%BrLEeminS,
     *  cnst%BrEemaxS, 
     *  cnst%BrUminSA, cnst%BrUmaxSA, cnst%BrTXTS,
     *  cnst%BrdUSA, cnst%BrdETXS, cnst%BrdES, cnst%BrUminSB,
     *  cnst%BrUmaxSB, cnst%BrdUSB,  cnst%BrES,
     *  cnst%BrUszSA, cnst%BrUszSB, cnst%how,
     *  cnst%NormS, cnst%NormPS, cnst%NormCS, cnst%NormSH
!         last 5 are special for brems normalization consts
!         for how=-1 at table creation time. (how->HowNormBrems)


      write(*, '(a)') 
     * ' constant for upper Seltzer Brems sampling table '
      write(*, '(1p, 5g15.7)')
     *  cnst%BrEeminS2, cnst%BrEgminS2, cnst%BrLEeminS2,
     *  cnst%BrEemaxS2, 
     *  cnst%BrUminSA2, cnst%BrUmaxSA2, cnst%BrTXTS2,
     *  cnst%BrdUSA2, cnst%BrdETXS2, cnst%BrdES2, cnst%BrUminSB2,
     *  cnst%BrUmaxSB2, cnst%BrdUSB2,  cnst%BrES2,
     *  cnst%BrUszSA2, cnst%BrUszSB2

      write(*, *)

      end

!     ****************************
      subroutine epwtBrCnst(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst

      write(*,'(a)') 
     * ' constant for Brems sampling table in p.s region'
      write(*,'(1p,5g15.7)') cnst%CompScrE,
     *  cnst%BremEgmin, cnst%BremEemin, cnst%BremLEemin,
     *  cnst%BremEeminLPM, cnst%BrScrE,
     *  cnst%BremUminLA, cnst%BremUmaxLA, cnst%BremTXTL,
     *  cnst%BremdULA, cnst%BremdETXL, cnst%BremdEL, cnst%BremUminLB,
     *  cnst%BremUmaxLB, cnst%BremdULB,  cnst%BremEsize,
     *  cnst%BremUszLA, cnst%BremUszLB

      write(*, *)
      end
!     ****************************
      subroutine epwtBrCnstH(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst


      write(*, '(a)') 
     * ' constant for Brems sampling table at high energies'
      write(*,'(1p,5g15.7)')
     *   cnst%BrEgminH, cnst%BrEe1H, cnst%BrLEe1H,
     *   cnst%BrneH, cnst%BrdU1H, cnst%BrdEH,
     *   cnst%BrEe2H, cnst%BrdU1H,  cnst%BrU1H,
     *   cnst%BrU2H,  cnst%Brnu1H,  cnst%BrU3H,
     *   cnst%BrU4H,  cnst%Brnu2H,  cnst%BrdVU2H,
     *   cnst%BrdU2H, cnst%BrneH2,  cnst%BrdEH2,
     *   cnst%BrEe2H2, cnst%BrPow

      write(*, *)

      end
!     ****************************
      subroutine epwtPrCnst(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst


      write(*,'(a)') 
     * ' constant for Pair  sampling table at low energies'
      write(*,'(1p, 5g15.7)')  cnst%PairEgmin, cnst%PairLEgmin,
     * cnst%PairNonSc, cnst%PrScrE,
     * cnst%PairEgmaxL,  cnst%PairTXTL, cnst%PairEsize,
     * cnst%PairUminLA,  cnst%PairUmaxLA, cnst%PairUszLA,
     * cnst%PairdULA, cnst%PairdETXL,  cnst%PairUminLB,
     * cnst%PairUmaxLB, cnst%PairUszLB, cnst%PairdULB,
     * cnst%PairdELA, cnst%PairdELB

      write(*, *)

      end
!     ****************************
      subroutine epwtPrCnstH(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst

      write(*,'(a)') 'constant used in Pair cre. sampling with LPM'
      write(*, '(1p,5g15.7)')
     *  cnst%PrEg1H, cnst%PrneH, cnst%PrdU1H, cnst%PrdEH,
     *  cnst%PrU1H,  cnst%PrU2H, cnst%Prnu1H, cnst%PrLEg1H,
     *  cnst%PrEg2H

      write(*, *)
      end

!     ****************************
      subroutine epwtmuNCnst(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst

      integer klena

      write(*,'(a)') 'constant used in muon Nuc. int. sampling'
      write(*,'(1p,5g15.7)')
     *   cnst%muNVmin,  cnst%muNdETX, cnst%muNdE,  cnst%muNEmin,
     *   cnst%muNEmax,  cnst%muNdU,   cnst%muNUsize, cnst%muNEsize,
     *   cnst%muNTXT, cnst%muNEmax1

      write(*, *)
      end
!     ****************************
      subroutine epwtmuPrCnst(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst
      
      write(*, '(a)') 'constant used in muon pair creation sampling'
      write(*, '(1p,5g15.7)')
     *  cnst%muPrVmin, cnst%muPrdETX, cnst%muPrdE, cnst%muPrEmin,
     *  cnst%muPrEmax, cnst%muPrdU, cnst%muPrUsize, cnst%muPrEsize,
     *  cnst%muPrTXT, cnst%muPrEmax1

      write(*, *)
      end
!     ****************************
      subroutine epwtmuBrCnst(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst
      
      write(*,'(a)') 'constant used in muon Brems sampling'
      write(*, '(1p, 5g15.7)' )
     *   cnst%muBrVmin, cnst%muBrdETX, cnst%muBrdE, cnst%muBrEmin,
     *   cnst%muBrEmax,  cnst%muBrdU, cnst%muBrUsize, cnst%muBrEsize,
     *   cnst%muBrTXT, cnst%muBrEmax1

      write(*, *)
      end
