      implicit none
#include "Zmedia.h"
#include "Zmass.h"
!
!          test Moller scattering diff.  function.
!          as a function of recoil electron energy
!          at fixed electron energy
! 
       type(epmedia):: media
      integer i, io, icon
      real*8 Ee, prob,  w, path, v, vm
      real*8 xm , epmollertx,  g, func, epMoller
      character*130 file
      character(len=12):: name

      io = 10
      
      call  cerrorMsg(
     * "Enter  electron energy(T.E 1. GeV),"//
     * " minimu recoil energy(K.E 100e-6 GeV)"//  
     * " and media name( BGO )", 1)
      Ee = 1.
      w = 200.0d-6
      name = "BGO"
      read(*, *)  Ee, w, name
      call cerrorMsg(name,  1)
      file = '$EPICSTOP/Data/Media/'//trim(name)

      call copenf(io, file, icon)
      call epReadTab(io, media)
      close(io)
      call epGetEffZA(media)
      call epStern(media)  ! argm. is 1
      vm =  1.01*w/((Ee-masele)/2)
      v = 0.5
      call epmollerp(media, Ee, w,  prob, path)
      g = Ee/masele
      write(*,'(a)') "#  v=KEr/KEe, ds/dv, ds/dv/(g/cm2)"
      do while ( v .gt. vm )
         func = epMoller(g, v)  ! bare Moller function
         write(*,'(1p, 3g13.4)') v,  func,
     *       func*media%sh%a/((g-1)*masele)/(1.-1./g/g)
         v = v/10.**0.02
      enddo
      end

      real*8 function epMoller(g, v)
      implicit none
!         gives diff. Moller xs= ds/dv (v=Ek/T0)
!     Ek; recoil  electron kinetic energy
      real*8 v
      real*8 g  ! incident electron gamma factor
!     T0:      its kinetic energy
!            
!        if this is multiplied by C/m(g-1), prob/(g/cm2)
!       is obtained.  C=0.153..Z/A/beta^2
!

      epMoller=((g-1)/g)**2 +(1./v-(2*g-1)/g**2)/v +
     *  (1./(1.-v) - (2*g-1.)/g**2)/(1.-v)
      end
