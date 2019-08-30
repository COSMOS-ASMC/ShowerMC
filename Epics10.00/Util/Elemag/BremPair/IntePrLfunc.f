!     integrate  pair function, epPairLowE, to get total pair
!     creation cross-section  for the no screening case.
!     This is a special version: at all given energies, no  screening
!     function as given by Bethe-Heitler (with some correction)
!     is used.  Only a single Z value medium can be specified.
!     (must have been registered in Epics). e.g, H2, W, Fe etc
!     Usage:  make -f IntePrLfunc.mk.  Then, execute
!     ./a.out   or  ./a.out > output_file 
!     Then, the input parameters are requested as
!       Enter media name, Eg1, Eg2, step
!     You may enter, e.g,      
!         H2  2e-3 20000 0.25
!     to get cross-section from 2MeV to 20 TeV with log10 step 0.25
!     proton target. (Cross-section is in mb for H and not for H2)
!     No screening for all snergies mean thet we suppose
!     typically gamma + p where p is fully ionized so that
!     the outer electron contribution must  not be included.
!     The problems are:
!  1) The cross-section here is no screening one but includes
!     the outer electron contribution.  To get the contribution
!     from the nucleus only,  Z/(Z+1) may be multiplied. 
!  2) The cross-sction for low Z, (say, H=>1),  is too bad and is
!     not connected to the partial screening cross-section, although
!     for larger Z, (say, W=>74), the value is quite good. So
!     we must normalize the cross-section by referring the
!     cross-section obtained by normal procedure at, say, 20 MeV
!     (use ./IntePair.sh for this.)
! 3)  This program vs   Normal, x1/2       XCOM (N.Pair)
!     MeV  mb            mb                   /(g/cm2)   mb
!     20  3.26          7.72   3.86           1.967e-3  3.29
!    200  7.34          17.41  8.71           4.31e-3   7.21  
!     Using above two corrections, we may be able to get a
!     reliable cross-section.  For the proton case,
!     correction by  1) is x 1/2.  For 2), x 2.38. So the
!     net correction factor is x 1.19
!     The value of this program is 3.26mb so the  answer is
!     3.88 mb      
!  4) Another way is to use XCOM data; they are in Epics/Data/Media/
!     say,  H2.xcom.  The column, "N.Pair" is by nucleus and "E.Pair"
!     is by outer electrons.  For the normalizetion, we may use
!     the value for "N.Pair". At 20 MeV,  the  value   
!     1.9673e-3 (1/(g/cm2))--> /(5.974e-4)--> mb => 3.29 mb
!     So we may multiply 3.29/3.26= 1.01.   Previous gave
!     3.88 mb : 18 % over.        
!
! *** Next module is used in epPairLowNorm in epPairLowE.f.
!     If real one is used, more programs must be
!     included; they are all not needed in this application.
!     So BPPS is defined here as dummy to avoid such complicated
!     situation.      
      module  BPPS
      real(8)::  epPair
      end      module  BPPS
#include "ZepicsBD.h"
      module minterf
      real(8):: Z
      real(8):: Egmeinterf
      end module minterf
      program main
      use minterf
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZepTrackp.h"
!
      type(epmedia):: media

      real*8 Ee,  Eg1, Eg2, step, Eg
      character*80 file
      character*24 name

      real(4):: rhoc
      real*8 xs           ! input  Eg/me
      real(8):: xmin, xmax, tprob
      integer  icon, io
!      real(8):: Z    ! Atom's Z
      real*8 x        ! input  Ee/Eg.   me/Eg =< x <= 1.-me/Eg   
      io = 10
      write(0,*) 'Enter media name, Eg1, Eg2, step' 
      read(*, *)  name, Eg1, Eg2, step
      write(0,*) trim(name), Eg1,  Eg2, step
      call epgetRhoc(name, name, rhoc)
      
      file = "$EPICSTOP/Data/Media/"//trim(name) 
      call cerrorMsg(file,  1)
      call copenf(io, file, icon)
      call epReadTab(io, media)
      write(0,*) ' media%Z=',  media%Z
      Z =  media%Z
      media%rhoc = rhoc
      Eg =max(Eg1, masele*2.0d0*1.00001d0)
      do while( Eg < Eg2*1.000001d0 )
         Egmeinterf = Eg/masele
         xmax = 1. - masele/Eg
         xmin = masele/Eg
         call epPrNonScTX( xmin, xmax, tprob)
         write(*,'(1p,4g13.3)') Eg, tprob

         Eg = Eg*10.0d0**step
      enddo
      end
      subroutine epPrNonScTX(xmin, xmax, tx)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
      
      real*8  xmin    ! input.  Ee/Eg min.
      real*8  xmax    ! input.  Ee/Eg max.
      real*8  tx      ! output.  Integral of pair function from
!         x= xmin to xmax
      real(8),external:: interf
      real*8  ans1, ans2
      real*8  d,  vt
      
      d=(xmax-xmin)/20.d0
      vt=xmax-d
      call k16pGaussLeg(interf, xmin, vt, 16,  ans1)
      call k16pGaussLeg(interf, vt, xmax, 16,  ans2)
      tx = ans1+ans2
      end
      function interf(x) result(ans)
      use minterf
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
      real(8),external ::  epPairLowE
      real(8),intent(in):: x
      real(8)::ans
      ans = epPairLowE(Z, Egmeinterf, x) 
      end
      
      
      
      
      
