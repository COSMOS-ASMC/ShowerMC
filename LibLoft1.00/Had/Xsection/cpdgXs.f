      module modpdgXs
      implicit none
      private 
!               if M below conflict in the user prog,
!               consider  alias like 
!              use  modpdgXs, Mpdg=>M
!               or 
!              use modpdgXs, only cSigmaT, csOfstu
      public :: cSigmaT, csOfstu, M

      real(8),parameter:: pi= 3.14159265358979323846d0
      real(8),parameter:: hbarc=197.327d-3  ! GeV%Fermi
      real(8),parameter:: M = 2.15   ! GeV
      real(8),parameter:: eta1 = 0.462, eta2 = 0.55
      real(8),parameter:: LB = pi*(hbarc/M)**2*10.0d0   ! mb 
      real(8),parameter:: delta = 0.003056, lambda = 1.630

      contains

      function csOfstu(p, Mp, Mt) result(ans)
!        s = (E+Mt)^2 - p**2 = 2*E*Mt + Mt**2 + Mp**2
      implicit none
      real(8),intent(in):: p   ! momentum  in GeV 
      real(8),intent(in):: Mp  ! proj. mass
      real(8),intent(in):: Mt  ! target mass
      real(8):: ans
      ans = (sqrt(p**2+Mp**2)*2 + Mt)*Mt + Mp**2
      end      function csOfstu

      function cSigmaT(a,b,p) result(ans) 
!         total cross-section by pdg (ed. 2012)
      implicit none
#include "Zmass.h"
      character(*),intent(in):: a  ! one of "p", "pb", "pi+","pi-"
                                   ! "pi0"
                              !      "K+","K-" "K0" "Sig", "g"
      character(*),intent(in):: b  ! one of "p", "n", "d"
      real(8),intent(in)::p  !  momentum of  a in GeV
      real(8):: ans       ! total cross section in mb.

      real(8):: LBt, sm, s

      real(8):: Z, Y1, Y2

      LBt = LB
      select case(a)
      case('p', 'pb')
         select case(b)
         case('p')
            Z = 34.71
            Y1 = 12.72
            Y2 = 7.35
            sm = (masp *2 + M)**2
            s = csOfstu(p,masp, masp)
         case('n')
            Z = 35.00
            Y1 = 12.19
            Y2 = 6.62
            sm = (masp + masn + M)**2
            s = csOfstu(p, masp, masn)
         case('d')   
            Z = 65.02
            Y1= 29.04
            Y2= 14.9
            sm = (masp + masd + M)**2
            s = csOfstu(p, masp, masd)
            LBt =LBt*lambda
         case default
            write(0,*)
     *       'b ', b, ' for csigma(a,b) is invalid'
            stop
         end select
      case('pi+', 'pi-', 'pi0')
         select case(b)
         case('p','n')
           Z = 19.02
           Y1 = 9.22
           Y2= 1.75
           sm = (maspic + masn + M)**2
           s = csOfstu(p, maspic, masp)
         case('d') 
           Z = 37.06
           Y1 = 18.28
           Y2 = 0.34
           sm = (maspic + masd + M)**2
           s = csOfstu(p, maspic, masd)
           LBt =LBt*lambda
         case default
            write(0,*)
     *       'b ', b, ' for csigma(a,b) is invalid'
            write(0,*) ' a is ',a
            stop
         end select  
      case('K+', 'K-', 'K0')
         select case(b)
         case('p')
            Z = 16.55
            Y1 =  4.02
            Y2 = 3.39
            sm = (maskc + masp + M)**2
            s = csOfstu(p, maskc, masp)
         case('n')
            Z = 16.49
            Y1 = 3.44
            Y2 = 1.82
            sm = (maskc + masn + M)**2
            s = csOfstu(p, maskc, masn)
         case('d')   
            Z = 32.34
            Y1 = 7.33
            Y2 = 5.59
            sm = (maskc + masd + M)**2
            s = csOfstu(p, maskc, masd)
            LBt =LBt*lambda
         case default
            write(0,*)
     *       'b ', b, ' for csigma(a,b) is invalid'
            write(0,*) ' a is ',a
            stop
         end select  
      case('Sig')
         select case(b)
         case('p')
            Z = 34.9
            Y1 = -55.0
            Y2 = -57.0
            sm = (massigmam + masp + M)**2
            s = csOfstu(p, massigmam, masp)
         case default
            write(0,*)
     *       'b ', b, ' for csigma(a,b) is invalid'
            write(0,*) ' a is ',a
            stop
         end select
      case('g')
         select case(b)
         case('p','n')
            Z = 34.71
            Y1 = 0.0128
            Y2= 0.
            sm = (masp+ masp + M)**2
            s = csOfstu(p, masp, masp)
            ans = (Z + LBt* (log(s/sm))**2)* delta
            sm = (masp + M)**2
            s = csOfstu(p, 0.d0, masp)
            ans = ans + Y1*(sm/s)**eta1
            return
         case('d')
            Z = 65.02
            Y1 = 0.0128
            Y2= 0.
            sm = (masp + masd + M)**2
            s = csOfstu(p, masp, masd)
            ans = (Z + lambda* LBt* (log(s/sm))**2 )* delta
            sm = ( masd + M)**2
            s = csOfstu(p, 0.d0, masd)
            ans = ans + Y1*(sm/s)**eta1
            return
         case('g')
            Z = 34.71
            Y1 = -0.034d-4
            Y2 = 0.
            sm = (masp + masp + M)**2
            s = csOfstu(p, masp, masp)
            ans = (Z +  LBt* (log(s/sm))**2 )* delta**2
            sm = ( masrho + M)**2
            s = csOfstu(p, 0.d0, masrho)
            ans = ans + Y1*(sm/s)**eta1
            return
         case default
            write(0,*)
     *       'b ', b, ' for csigma(a,b) is invalid'
            write(0,*) ' a is ',a
            stop
         end select
      case default
         write(0,*)
     *        'a ', a, ' for csigma(a,b) is invalid'
            write(0,*) ' b is ',b
            stop
      end select
      ans = Z + LBt* (log(s/sm))**2 + Y1*(sm/s)**eta1
      if( a == "pb" .or. a =="pi-" .or. a == "K-") then
         ans = ans + Y2*(sm/s)**eta2
      else
         ans = ans - Y2*(sm/s)**eta2
      endif
      end   function cSigmaT
      end   module modpdgXs

!      program  main
!      use modpdgXs
!      implicit none
!#include "Zmass.h"
!      character(4):: a = "g"
!      character(4):: b = "d"
!      real(8):: p,s
!      integer:: i
!      p = 5.
!      do while (p < 1.e8 ) 
!         s = csOfstu(p, 0.d0, masd)
!         write(*,*) p, cSigmaT(a, b, p), sqrt(s)
!         p = p*10.**0.1
!      enddo
!      end
