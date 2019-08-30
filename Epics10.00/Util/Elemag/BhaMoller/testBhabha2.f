      implicit none
#include "Zmedia.h"
#include "Zmass.h"
!
!          test Bhabha scattering diff.  function.
!          as a function of recoil electron energy
!          at fixed positrotn energy
! 
       type(epmedia):: media
      integer i, io, icon
      real*8 Ee, prob,  w, path, v, vm
      real*8 g, func, epBhabha,epBhabhaRL, func2
      character*130 file
      character*12  name

      io = 10
      
      call  cerrorMsg(
     * "Enter  positron energy(0.1 T.E  GeV),"//
     * " minimu recoil energy(K.E 100d-6 GeV)"//  
     * " and media name('BGO')", 1)
      Ee = 0.1
      w = 100.0d-6
      name = "BGO"
      read(*, *)  Ee, w, name

      call cerrorMsg(name,  1)
      file = '$EPICSTOP/Data/Media/'//trim(name)
      call copenf(io, file, icon)
      call epReadTab(io, media)
      close(io)
      call epGetEffZA(media)
      call epStern(media)   ! now arg. is 1
      vm =  1.01*w/((Ee-masele))
      v = 1.0
      call epbhabhap(media, Ee, w,  prob, path)
      g = Ee/masele

      write(*,*) "# two functions of diff. authors: v=E+'/E+" 
      write(*,*)
     * "#  v     f1     f2     f1*Ek*C/beta2   f2*Ek*C/beta2 "
      do while ( v .gt. vm )
         func = epBhabha(g, v)  ! bare Bhabha function
         func2 = epBhabhaRL(g, v)  ! bare Bhabha function
         write(*,'(1p,5g12.4)') v,  func, func2,
     *       func*media%sh%a/((g-1)*masele)/(1.-1./g/g),
     *       func2*media%sh%a/((g-1)*masele)/(1.-1./g/g)
         v = v/10.**0.02
      enddo
      end

      real*8 function epBhabha(g, v)
      implicit none
!         gives diff. Bhabha xs= ds/dv (v=Ek/T0)
!       Expressions are used in many literatures
!     Ek; recoil  electron kinetic energy
      real*8 v
      real*8 g  ! incident positron gamma factor
!     T0:      its kinetic energy
!            
!        if this is multiplied by C/m(g-1), prob/(g/cm2)
!       is obtained.  C=0.153..Z/A/beta2
!
      real*8 B1, B2, B3, B4, y, beta2

      y = 1./(g+1.)
      B1 = 2.- y**2
      B4 =  (1.-2*y)**3
      B3 = B4+(1.-2*y)**2
      B2 = (1.-2*y)*(3.+y**2)
      beta2 = 1.- 1./g/g

      epBhabha=1./v**2 - beta2*(B1/v -B2 + B3*v -B4*v*v)
      end
      real*8 function epBhabhaRL(g, v)
      implicit none
!         gives diff. Bhabha xs= ds/dv (v=Ek/T0)
!      Expressions are based on Rohrlich & Carlso,
!        Phys.Rev.vol.93, No.1 (1954) p.38
!       Apparent formulat looks very much different from
!       the epBhabha, and not easy to prove that the both
!       are equivalent, numerical results coincide.
!
!     Ek; recoil  electron kinetic energy
      real*8 v
      real*8 g  ! incident positron gamma factor
!     T0:      its kinetic energy
!            
!        if this is multiplied by C/m(g-1), prob/(g/cm2)
!       is obtained.  C=0.153..Z/A/beta2
!
      real*8 B1, B2, B3, B4, y, beta2
      real*8 u
      y = 1./(g+1.)
      u = g-1.0
      beta2 = 1.- 1./g/g

      epBhabhaRL=1./v**2 - beta2/v + (u/g)**2/2
     * - u*y*( (g+2)/g/v - 2*beta2 +
     *   v*(u/g)**2 ) 
     * + (u*y)**2*(0.5 + 1./g + 1.5/g/g
     * - (u/g)**2*v*(1.-v))
      end
