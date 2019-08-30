!     This is a modified version of epwtSmpTbl.f; 
!   The output format is data statement usable in Fortran
!   program.  (For making special sampling routine dedicated
!      to a particular media; used for Air in Cosmos)
!  Note:  subroutine names are the same as used in
!         epwtSmpTbl;
!    but  argument is little bit diff. 
!
!To convert declaration part in epwtSmpTbl, awk -f nextone ..
!      was used.
!{for(i=1; i<=NF;i++){
! gsub(/,/,"",$i);
!
! if(length($i) > 1) {
!	 x= $i;
!	 gsub(/cnst./,"", x);
!	 print "      write(*, *)'      real*8 ", x, "'";
!  }
! }
!}
!
! To make data statement, next awk is used
!
!{for(i=1; i<=NF;i++){
! gsub(/,/,"",$i);
! if(length($i) > 1) {
!	 x= $i;
!	 gsub(/cnst./,"", x);
!	 print "      write(*, *)'      data ", x,"/',",$i,",'/'";
!  }
! }
!}
!
!
!
!     *******************************
      subroutine epwtmedia(media)
      implicit none
!       print some of basic media information(tbl, cnst are not
!      included)
!        should be called after epGetEffZA, epExpot
!

#include "Zmedia.h"
       type(epmedia):: media

!
      write(*, *) '      real*8  Ipot '
      write(*, *) '      real*8  X0'
      write(*, *) '      real*8  X0g'
      write(*, *) '      real*8  Aeff, Zeff, Z2byAeff'
      write(*, *) '      real*8  Z2eff, ZbyAeff'
      write(*, *) '      real*8  wp, nd'
      write(*, *) '      real*8  s1'
      write(*, *) '      real*8  gtocm'
      write(*, *) '      real*8  mbtoPgrm'
      write(*, *) '      real*8  mbtoPcm'
      write(*, *) '      real*8  mbtoPX0'

      write(*, *)
     *      'c      Excitation (Ionization) potential(GeV)'
      write(*, *)'      data Ipot/', media%I, '/'

      write(*, *) 
     *      'c        radiation length; cm and g/cm^2'
      write(*, *)'      data X0/', media%X0, '/,',
     *           '  X0g/',  media%X0g, '/'


      write(*,*) 'c        <Aeff>,   <Zeff>,  <Z2/A> '
      write(*,*) '      data Aeff/',  media%Aeff,  '/,  Zeff/',
     *           media%Zeff, '/'
      write(*,*) '      data  Z2byAeff/', media%Z2byAeff, '/'

      write(*,*) 'c        <Z2eff> and <Z/A>='
      write(*,*) '      data Z2eff/', media%Z2eff, 
     *           '/,  ZbyAeff/', media%ZbyAeff, '/'

      write(*,*) 'c           plasma energy GeV; No density (/cm^3)'
      write(*,*) '      data wp/', media%wp,
     *           '/,  nd/', media%nd, '/'


      write(*,*) "c            Migdal's s1"
      write(*,*) '      data   s1/',  media%s1, '/'

      write(*,*) 'c             Conversion factors:'
      write(*,*) 'c      g/cm^2 to cm='
      write(*,*) '      data   gtocm/', media%gtocm,'/'

      write(*,*) 'c       mb to /(g/cm2); mb to /cm='
      write(*,*) '      data   mbtoPgrm/', media%mbtoPgrm,'/' 
      write(*,*) '      data   mbtoPcm/', media%mbtoPcm, '/'

      write(*,*)'c          mb to /X0='
      write(*,*)'       data  mbtoPX0/', media%mbtoPX0, '/'

      write(*,*)
     *       'c      If a comound is specified as an Atom, use the'
      write(*,*)
     *       'c      following instead: ', 
      write(*,*) 'c   ',
     * media%mbtoPgrm2, media%mbtoPcm2, media%mbtoPX02

      write(*,*)
      write(*,*) 'c     pi x Re**2 * N* Z /A *X0g  = 0.15 Z/A*X0g=',
     *           media%basearea
!
!   ************************************************************
!                   next consts cannot be fixed without fixing
!            knockon minimum energy.
!   ************************************************************
!      write(*, *) 
!      write(*, *) 'c   **** Sternheimers consts'
!      write(*, *) 'c    a, b, c =', media.sh.a, media.sh.b, media.sh.c
!      write(*, *) 'c     x0, x1 =', media.sh.x0, media.sh.x1
!      write(*, *) 'c     sa, xa =', media.sh.sa, media.sh.xa
!      write(*, *) 'c     minimum dEdx_restricted=',media.dEdxatp3m,
!      write(*, *) 'c     GeV/(g/cm^2)  with Tcut=',media.sh.tcut, ' GeV'
!      write(*, *) 
      write(*, *) 'c   **** photoelectric consts'
      write(*, *)
     *  'c   b0, b1, b2=', media%pe%b0, media%pe%b1, media%pe%b2
      write(*, *) 'c   fa, a, p=', media%pe%fa, media%pe%a, media%pe%p
      write(*, *) 'c   l   ek=', media%pe%l, media%pe%ek
      write(*, *)
!          separator; when having read basic media later, skip next to
!          this separator
      write(*, *)'c         #-#-#-#-#-#-#-#-#-#-#-#-'
      end

!     ******************************* 
      subroutine epwt1dTbl(com, name,  tbl, size)
!           print  total cross section
      implicit none
      integer size
      real*8 tbl(size)
      character*(*) name, com
      integer klena



      write(*, *)'c   ',  com(1:klena(com)) 

      call kmkDataStm1(tbl, size,  name(1:klena(name)), 
     *       'g12.5', 14)
      write(*,*)
      end
!     ******************************************
      subroutine epwt2dTbl(com, name, tbl, sizeu,  sizee)
!  
      implicit none
      integer sizeu,sizee
      real*8 tbl(sizeu,  sizee)
      character*(*) com
      character*(*) name !  variable name for the array.
      integer  klena

      write(*, *)
      write(*, *)'c  ', com(1:klena(com))

      call kmkDataStm2a(tbl, sizeu, sizee, name(1:klena(name)),
     *         'f11.5', 11)

!      do i = 1, sizee
!         call kmkDataStm2b(tbl(1, i), 
!         write(*, '(7f11.5)') ( tbl(iu, i), iu = 1, sizeu )
!         write(*, *)
!      enddo
      write(*,*)
      end
!     ****************************
      subroutine epwtBrCnstS(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst
      write(*,*)
      write(*, *) 
     * 'c    consts for Seltzer Brems sampling table '

      write(*, *)'      real*8  BrEeminS '
      write(*, *)'      real*8  BrEgminS '
      write(*, *)'      real*8  BrLEeminS '
      write(*, *)'      real*8  BrEemaxS '
      write(*, *)'      real*8  BrUminSA '
      write(*, *)'      real*8  BrUmaxSA '
      write(*, *)'      real*8  BrTXTS '
      write(*, *)'      real*8  BrdUSA '
      write(*, *)'      real*8  BrdETXS '
      write(*, *)'      real*8  BrdES '
      write(*, *)'      real*8  BrUminSB '
      write(*, *)'      real*8  BrUmaxSB '
      write(*, *)'      real*8  BrdUSB '
      write(*, *)'      real*8  BrES '
      write(*, *)'      real*8  BrUszSA '
      write(*, *)'      real*8  BrUszSB '
      write(*, *)
      write(*, *) 
     * 'c    consts for Seltzer Brems sampling table '

      write(*, *)'      data  BrEeminS /', cnst%BrEeminS ,'/'
      write(*, *)'      data  BrEgminS /', cnst%BrEgminS ,'/'
      write(*, *)'      data  BrLEeminS /', cnst%BrLEeminS ,'/'
      write(*, *)'      data  BrEemaxS /', cnst%BrEemaxS ,'/'
      write(*, *)'      data  BrUminSA /', cnst%BrUminSA ,'/'
      write(*, *)'      data  BrUmaxSA /', cnst%BrUmaxSA ,'/'
      write(*, *)'      data  BrTXTS /', cnst%BrTXTS ,'/'
      write(*, *)'      data  BrdUSA /', cnst%BrdUSA ,'/'
      write(*, *)'      data  BrdETXS /', cnst%BrdETXS ,'/'
      write(*, *)'      data  BrdES /', cnst%BrdES ,'/'
      write(*, *)'      data  BrUminSB /', cnst%BrUminSB ,'/'
      write(*, *)'      data  BrUmaxSB /', cnst%BrUmaxSB ,'/'
      write(*, *)'      data  BrdUSB /', cnst%BrdUSB ,'/'
      write(*, *)'      data  BrES /', cnst%BrES ,'/'
      write(*, *)'      data  BrUszSA /', cnst%BrUszSA ,'/'
      write(*, *)'      data  BrUszSB /', cnst%BrUszSB ,'/'

      write(*, *)
      end

!     ****************************
      subroutine epwtBrCnst(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst

      write(*, *) 
     * 'c    consts for Brems sampling table at low energies'
      write(*, *)'      real*8  CompScrE '
      write(*, *)'      real*8  BremEgmin '
      write(*, *)'      real*8  BremEemin '
      write(*, *)'      real*8  BremLEemin '
      write(*, *)'      real*8  BremEemaxL '
      write(*, *)'      real*8  BrScrE '
      write(*, *)'      real*8  BremUminLA '
      write(*, *)'      real*8  BremUmaxLA '
      write(*, *)'      real*8  BremTXTL '
      write(*, *)'      real*8  BremdULA '
      write(*, *)'      real*8  BremdETXL '
      write(*, *)'      real*8  BremdEL '
      write(*, *)'      real*8  BremUminLB '
      write(*, *)'      real*8  BremUmaxLB '
      write(*, *)'      real*8  BremdULB '
      write(*, *)'      real*8  BremEsize '
      write(*, *)'      real*8  BremUszLA '
      write(*, *)'      real*8  BremUszLB '
      write(*, *)
      write(*, *) 
     * 'c    consts for Brems sampling table at low energies'

      write(*, *)'      data  CompScrE /', cnst%CompScrE ,'/'
      write(*, *)'      data  BremEgmin /', cnst%BremEgmin ,'/'
      write(*, *)'      data  BremEemin /', cnst%BremEemin ,'/'
      write(*, *)'      data  BremLEemin /', cnst%BremLEemin ,'/'
      write(*, *)'      data  BremEemaxL /', cnst%BremEemaxL ,'/'
      write(*, *)'      data  BrScrE /', cnst%BrScrE ,'/'
      write(*, *)'      data  BremUminLA /', cnst%BremUminLA ,'/'
      write(*, *)'      data  BremUmaxLA /', cnst%BremUmaxLA ,'/'
      write(*, *)'      data  BremTXTL /', cnst%BremTXTL ,'/'
      write(*, *)'      data  BremdULA /', cnst%BremdULA ,'/'
      write(*, *)'      data  BremdETXL /', cnst%BremdETXL ,'/'
      write(*, *)'      data  BremdEL /', cnst%BremdEL ,'/'
      write(*, *)'      data  BremUminLB /', cnst%BremUminLB ,'/'
      write(*, *)'      data  BremUmaxLB /', cnst%BremUmaxLB ,'/'
      write(*, *)'      data  BremdULB /', cnst%BremdULB ,'/'
      write(*, *)'      data  BremEsize /', cnst%BremEsize ,'/'
      write(*, *)'      data  BremUszLA /', cnst%BremUszLA ,'/'
      write(*, *)'      data  BremUszLB /', cnst%BremUszLB ,'/'
      write(*, *) 
      end
!     ****************************
      subroutine epwtBrCnstH(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst
      write(*, *) 
     * 'c   consts for Brems sampling table at high energies'

      write(*, *)'      real*8  BrEgminH '
      write(*, *)'      real*8  BrEe1H '
      write(*, *)'      real*8  BrLEe1H '
      write(*, *)'      real*8  BrneH '
      write(*, *)'      real*8  BrdU1H '
      write(*, *)'      real*8  BrdEH '
      write(*, *)'      real*8  BrEe2H '
!      write(*, *)'      real*8  BrdU1H '  ! doubly def.
      write(*, *)'      real*8  BrU1H '
      write(*, *)'      real*8  BrU2H '
      write(*, *)'      real*8  Brnu1H '
      write(*, *)'      real*8  BrU3H '
      write(*, *)'      real*8  BrU4H '
      write(*, *)'      real*8  Brnu2H '
      write(*, *)'      real*8  BrdVU2H '
      write(*, *)'      real*8  BrdU2H '
      write(*, *)'      real*8  BrneH2 '
      write(*, *)'      real*8  BrdEH2 '
      write(*, *)'      real*8  BrEe2H2 '
      write(*, *)'      real*8  BrPow '

      write(*, *) 
     * 'c   consts for Brems sampling table at high energies'
      write(*, *)'      data  BrEgminH /', cnst%BrEgminH ,'/'
      write(*, *)'      data  BrEe1H /', cnst%BrEe1H ,'/'
      write(*, *)'      data  BrLEe1H /', cnst%BrLEe1H ,'/'
      write(*, *)'      data  BrneH /', cnst%BrneH ,'/'
      write(*, *)'      data  BrdU1H /', cnst%BrdU1H ,'/'
      write(*, *)'      data  BrdEH /', cnst%BrdEH ,'/'
      write(*, *)'      data  BrEe2H /', cnst%BrEe2H ,'/'
!      write(*, *)'      data  BrdU1H /', cnst.BrdU1H ,'/'  ! doubly def.
      write(*, *)'      data  BrU1H /', cnst%BrU1H ,'/'
      write(*, *)'      data  BrU2H /', cnst%BrU2H ,'/'
      write(*, *)'      data  Brnu1H /', cnst%Brnu1H ,'/'
      write(*, *)'      data  BrU3H /', cnst%BrU3H ,'/'
      write(*, *)'      data  BrU4H /', cnst%BrU4H ,'/'
      write(*, *)'      data  Brnu2H /', cnst%Brnu2H ,'/'
      write(*, *)'      data  BrdVU2H /', cnst%BrdVU2H ,'/'
      write(*, *)'      data  BrdU2H /', cnst%BrdU2H ,'/'
      write(*, *)'      data  BrneH2 /', cnst%BrneH2 ,'/'
      write(*, *)'      data  BrdEH2 /', cnst%BrdEH2 ,'/'
      write(*, *)'      data  BrEe2H2 /', cnst%BrEe2H2 ,'/'
      write(*, *)'      data  BrPow /', cnst%BrPow ,'/'
      write(*, *)
      end
!     ****************************
      subroutine epwtPrCnst(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst

      write(*, *) 
     * 'c  consts for Pair  sampling table at low energies'

      write(*, *)'      real*8  PairEgmin '
      write(*, *)'      real*8  PairLEgmin '
      write(*, *)'      real*8  PairNonSc '
      write(*, *)'      real*8  PrScrE '
      write(*, *)'      real*8  PairEgmaxL '
      write(*, *)'      real*8  PairTXTL '
      write(*, *)'      real*8  PairEsiz '
      write(*, *)'      real*8  PairUminLA '
      write(*, *)'      real*8  PairUmaxLA '
      write(*, *)'      real*8  PairUszLA '
      write(*, *)'      real*8  PairdULA '
      write(*, *)'      real*8  PairdETXL '
      write(*, *)'      real*8  PairUminLB '
      write(*, *)'      real*8  PairUmaxLB '
      write(*, *)'      real*8  PairUszLB '
      write(*, *)'      real*8  PairdUL '
      write(*, *)'      real*8  PairdELA '
      write(*, *)'      real*8  PairdELB '

      write(*, *) 
     * 'c  consts for Pair  sampling table at low energies'

      write(*, *)'      data  PairEgmin /', cnst%PairEgmin ,'/'
      write(*, *)'      data  PairLEgmin /', cnst%PairLEgmin ,'/'
      write(*, *)'      data  PairNonSc /', cnst%PairNonSc ,'/'
      write(*, *)'      data  PrScrE /', cnst%PrScrE ,'/'
      write(*, *)'      data  PairEgmaxL /', cnst%PairEgmaxL ,'/'
      write(*, *)'      data  PairTXTL /', cnst%PairTXTL ,'/'
      write(*, *)'      data  PairEsize /', cnst%PairEsize ,'/'
      write(*, *)'      data  PairUminLA /', cnst%PairUminLA ,'/'
      write(*, *)'      data  PairUmaxLA /', cnst%PairUmaxLA ,'/'
      write(*, *)'      data  PairUszLA /', cnst%PairUszLA ,'/'
      write(*, *)'      data  PairdULA /', cnst%PairdULA ,'/'
      write(*, *)'      data  PairdETXL /', cnst%PairdETXL ,'/'
      write(*, *)'      data  PairUminLB /', cnst%PairUminLB ,'/'
      write(*, *)'      data  PairUmaxLB /', cnst%PairUmaxLB ,'/'
      write(*, *)'      data  PairUszLB /', cnst%PairUszLB ,'/'
      write(*, *)'      data  PairdULB /', cnst%PairdULB ,'/'
      write(*, *)'      data  PairdELA /', cnst%PairdELA ,'/'
      write(*, *)'      data  PairdELB /', cnst%PairdELB ,'/'
      write(*, *)

      end
!     ****************************
      subroutine epwtPrCnstH(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst
      
      write(*, *) 'c   consts  used in Pair cre. sampling with LPM'
      write(*, *)'      real*8  PrEg1H '
      write(*, *)'      real*8  PrneH '
      write(*, *)'      real*8  PrdU1H '
      write(*, *)'      real*8  PrdEH '
      write(*, *)'      real*8  PrU1H '
      write(*, *)'      real*8  PrU2H '
      write(*, *)'      real*8  Prnu1H '
      write(*, *)'      real*8  PrLEg1H '
      write(*, *)'      real*8  PrEg2H '
      write(*, *)
      write(*, *)'c    consts used in Pair cre. sampling with LPM'
      write(*, *)'      data  PrEg1H /', cnst%PrEg1H ,'/'
      write(*, *)'      data  PrneH /', cnst%PrneH ,'/'
      write(*, *)'      data  PrdU1H /', cnst%PrdU1H ,'/'
      write(*, *)'      data  PrdEH /', cnst%PrdEH ,'/'
      write(*, *)'      data  PrU1H /', cnst%PrU1H ,'/'
      write(*, *)'      data  PrU2H /', cnst%PrU2H ,'/'
      write(*, *)'      data  Prnu1H /', cnst%Prnu1H ,'/'
      write(*, *)'      data  PrLEg1H /', cnst%PrLEg1H ,'/'
      write(*, *)'      data  PrEg2H /', cnst%PrEg2H ,'/'

      write(*, *)
      end

!     ****************************
      subroutine epwtmuNCnst(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst
      
      write(*, *)'c     consts used in muon Nuc. int. sampling'
      write(*, *)'      real*8  muNVmin '
      write(*, *)'      real*8  muNdETX '
      write(*, *)'      real*8  muNdE '
      write(*, *)'      real*8  muNEmin '
      write(*, *)'      real*8  muNEmax '
      write(*, *)'      real*8  muNdU '
      write(*, *)'      real*8  muNUsize '
      write(*, *)'      real*8  muNEsize '
      write(*, *)'      real*8  muNTXT '
      write(*, *)'      real*8  muNEmax1 '
      write(*, *)
      write(*, *)'c   consts used in muon Nuc. int. sampling'
      write(*, *)'      data  muNVmin /', cnst%muNVmin ,'/'
      write(*, *)'      data  muNdETX /', cnst%muNdETX ,'/'
      write(*, *)'      data  muNdE /', cnst%muNdE ,'/'
      write(*, *)'      data  muNEmin /', cnst%muNEmin ,'/'
      write(*, *)'      data  muNEmax /', cnst%muNEmax ,'/'
      write(*, *)'      data  muNdU /', cnst%muNdU ,'/'
      write(*, *)'      data  muNUsize /', cnst%muNUsize ,'/'
      write(*, *)'      data  muNEsize /', cnst%muNEsize ,'/'
      write(*, *)'      data  muNTXT /', cnst%muNTXT ,'/'
      write(*, *)'      data  muNEmax1 /', cnst%muNEmax1 ,'/'

      write(*, *)
      end
!     ****************************
      subroutine epwtmuPrCnst(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst
      
      write(*, *)'c    consts used in muon pair creation sampling'
      write(*, *)'      real*8  muPrVmin '
      write(*, *)'      real*8  muPrdETX '
      write(*, *)'      real*8  muPrdE '
      write(*, *)'      real*8  muPrEmin '
      write(*, *)'      real*8  muPrEmax '
      write(*, *)'      real*8  muPrdU '
      write(*, *)'      real*8  muPrUsize '
      write(*, *)'      real*8  muPrEsize '
      write(*, *)'      real*8  muPrTXT '
      write(*, *)'      real*8  muPrEmax1 '
      write(*, *)
      write(*, *)'c    consts used in muon pair creation sampling'
      write(*, *)'      data  muPrVmin /', cnst%muPrVmin ,'/'
      write(*, *)'      data  muPrdETX /', cnst%muPrdETX ,'/'
      write(*, *)'      data  muPrdE /', cnst%muPrdE ,'/'
      write(*, *)'      data  muPrEmin /', cnst%muPrEmin ,'/'
      write(*, *)'      data  muPrEmax /', cnst%muPrEmax ,'/'
      write(*, *)'      data  muPrdU /', cnst%muPrdU ,'/'
      write(*, *)'      data  muPrUsize /', cnst%muPrUsize ,'/'
      write(*, *)'      data  muPrEsize /', cnst%muPrEsize ,'/'
      write(*, *)'      data  muPrTXT /', cnst%muPrTXT ,'/'
      write(*, *)'      data  muPrEmax1 /', cnst%muPrEmax1 ,'/'
      write(*, *)
      end
!     ****************************
      subroutine epwtmuBrCnst(cnst)
      implicit none
#include  "ZbpSample.h"
       type(SmpCnst)::  cnst
      
      write(*, *)'c   consts used in muon Brems sampling'
      write(*, *)'      real*8  muBrVmin '
      write(*, *)'      real*8  muBrdETX '
      write(*, *)'      real*8  muBrdE '
      write(*, *)'      real*8  muBrEmin '
      write(*, *)'      real*8  muBrEmax '
      write(*, *)'      real*8  muBrdU '
      write(*, *)'      real*8  muBrUsize '
      write(*, *)'      real*8  muBrEsize '
      write(*, *)'      real*8  muBrTXT '
      write(*, *)'      real*8  muBrEmax1 '
      write(*, *)
      write(*, *)'c      consts used in muon Brems sampling'
      write(*, *)'      data  muBrVmin /', cnst%muBrVmin ,'/'
      write(*, *)'      data  muBrdETX /', cnst%muBrdETX ,'/'
      write(*, *)'      data  muBrdE /', cnst%muBrdE ,'/'
      write(*, *)'      data  muBrEmin /', cnst%muBrEmin ,'/'
      write(*, *)'      data  muBrEmax /', cnst%muBrEmax ,'/'
      write(*, *)'      data  muBrdU /', cnst%muBrdU ,'/'
      write(*, *)'      data  muBrUsize /', cnst%muBrUsize ,'/'
      write(*, *)'      data  muBrEsize /', cnst%muBrEsize ,'/'
      write(*, *)'      data  muBrTXT /', cnst%muBrTXT ,'/'
      write(*, *)'      data  muBrEmax1 /', cnst%muBrEmax1 ,'/'

      write(*, *)
      end


