!    ******************************************************************
!    *                                                                *
!    *  prompt muon decay  mode is important.
!    *     others are roughly introduced to keep the energy conservation
!    *
!
!    ******************************************************************
! ******** note:  c2bdcy:   boost to lab is done inside
!                 cnbdcy:   boost must be done outside.
      subroutine cdDecay(pj, a, np)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"

      type(ptcl):: pj   ! input. demeson
      integer np                ! output. # of ptcls stored in a
      type(ptcl):: a(*)  ! output. to store produced ptcls

      real*8 u, sumbr


      call rndc(u)
      if(pj%charge .ne. 0) then
         if(u .lt. 0.093d0 ) then
!                   K0b + mu + neu    9.3%
            call cdDecay1( pj, a, np)
         elseif(u .lt. 0.132d0) then
!                   K + pi + mu + neumu 3.9 %  13.2
            call cdDecay2(pj,  a, np)
         else
            sumbr=0.4179  ! the total % below)/100. 
!                    K 2pi         9.22        9.22
!                    K0b + e  + neu   8.6 %   17.82
!                     K0s + 2pi    6.8        24.62
!                     K+3pi        6.0        30.62
!                     K + pi+ e+ neue  4.1 %  34.72  
!                     K0s + 3pi     3.02      37.74   
!                     K0L + pi      1.46      39.20
!                     K0s + pi      1.45      40.65
!                     3pi + pi0     1.14      41.79
            call rndc(u)
!            u=u*sumbr  ! 2013.  enu prob too big
            if(u .lt.  0.0922)  then
!                D+    K- pi+ p+         9.22        9.22
               call cdDecay3(pj, a, np)
            elseif(u .lt. 0.1782) then
!                    K0b  e+  neu   8.6 %   17.82
               call cdDecay4(pj, a, np)
            elseif(u .lt. 0.2462) then
!                     K0s + 2pi    6.8        24.62
               call cdDecay5( pj, a, np)
            elseif(u .lt. 0.3062) then
!                     K+3pi        6.0        30.62
               call cdDecay6(pj, a, np)
            elseif(u .lt. 0.3472) then
!                     K + pi+ e+ neue  4.1 %  34.72  
               call cdDecay7( pj, a, np)
            elseif( u .lt. 0.3774) then
!                     K0s + 3pi     3.02      37.74   
               call cdDecay8(pj, a, np)
            elseif( u .lt. 0.3920) then
!                     K0L + pi      1.46      39.20
               call cdDecay9(pj, a, np)
            elseif( u .lt. 0.4065) then
!                     K0s + pi      1.45      40.65
               call cdDecay10(pj, a, np)
            elseif( u < 0.4179) then
!                     3pi + pi0     1.14      41.79
               call cdDecay11(pj, a, np)
            else    ! 2013.Mar 2. neglect other
               np =0
            endif
         endif
      else
!            D0        
!                 K^- + mu^+ + neu     3.31   0.0331
!                 K^- + e^+  + neue    3.55   0.06863
!                 K*+ e^+ + neue       2.16   
!               treat as ==> K+ pi + e^+ + neue and combine with
!                 K^-+ pi0+e^+ +neue   1.6    
!                 K^0bar + pi^-+e^+ +neue 2.7  0.1332
!c                 K*+ mu^+ + neu       1.91
!               treat as K+pi+ mu^+ + neu      0.1523
!                 -----------------------
!                 K- pi+ pi-       13.9       0.2913
!                 K-+2pi+ pi-       8.08      0.3721
!                 K-2pi+ pi- pi0    4.2       0.4141
!              K- rho+=>K- pi+ pi0  10.8      0.5221

!                 K0s + pi+ pi- pi0 5.2       0.5741

! ............

         call rndc(u) 
         if(u .lt. 0.0331) then
!              D0    K-  mu+  neu  3.31
            call cdDecay20(pj, a, np)
         elseif( u < 0.06863) then
!                    K + e  + neue 3.55
            call cdDecay25( pj, a, np)
         elseif( u < 0.1332) then
            call cdDecay28(pj, a, np)
         elseif( u< 0.2523) then
!               treat as K+pi+ mu^+ + neu      0.1523
            call cdDecay29(pj, a, np)
         else
!            call rndc(u)
!            sumbr =0.3695
!            u=u*sumbr    ! avoid too large other branch
            if(u .lt. 0.2913) then
!                          K- pi+ pi-       13.9   13.9
               call cdDecay21(pj, a, np)
            elseif(u .lt. 0.3721 ) then
!                          K-+2pi+ pi-     8.08
               call cdDecay22( pj, a, np)
            elseif(u .lt. 0.4141) then
!                          K-2pi+ pi- pi0    4.2   23.5
               call cdDecay23( pj, a, np)
            elseif(u .lt. 0.5221) then
!                      K- pi+ pi0       10.8   27.39
               call cdDecay24( pj, a, np)
            elseif(u .lt. 0.5741 ) then
!                          K0s + pi+ pi- pi0  5.2
               call cdDecay27(pj, a, np)
            else ! neglect all others
               np = 0
            endif
         endif
      endif
      end
!     ********************************
      subroutine cdDecay1( pj, a, np)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
!        D+->  K0b + mu + neu    9.3%
      type(ptcl):: pj   ! input. demeson
      integer np                ! output. # of ptcls stored in a
      type(ptcl):: a(*)  ! output. to store produced ptcls

      integer muchg, nusubc, kchg, ksubc, icon, i, echg
      integer pichg

      real*8 u, w

      call rndc(u)
      if(u .lt. .50) then
         ksubc = k0s
      else
         ksubc = k0l
      endif
      muchg = pj%charge
      nusubc = muchg
!           neue
      call cmkptc(kneue, nusubc, 0, a(1))
      call cmkptc(kmuon, 0,  muchg,  a(2))
      call cmkptc(kkaon, ksubc, 0, a(3))
!            3  body pure phase space
      call cnbdcy(3, pj%mass, a,  0, w, icon)
      np = 3
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!     *******************************************
!        D+    K-  pi+ mu+  neumu 3.9 %  13.2
      entry cdDecay2(pj, a, np)

      kchg =-pj%charge
      nusubc = regptcl
      call cmkptc(kkaon, 0, kchg, a(1))
      call cmkptc(kpion, 0, -kchg, a(2))
      call cmkptc(kmuon, 0, -kchg, a(3))
      call cmkptc(kneumu, nusubc, 0, a(4))

      call cnbdcy(4, pj%mass, a,  0, w, icon)
      np = 4
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!     ***********************************
!            D+   K- pi+ p+         9.22        9.22
      entry cdDecay3(pj, a, np)


      kchg = -pj%charge
      call cmkptc(kkaon, 0, kchg, a(1))
      call cmkptc(kpion, 0, -kchg, a(2))
      call cmkptc(kpion, 0, -kchg, a(3))

      call cnbdcy(3, pj%mass, a,  0, w, icon)
      np = 3
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!*****************************************
!             K0b  e+  neu   8.6 %   17.82
      entry cdDecay4(pj, a, np)


      ksubc = k0s
      echg = 1
      call rndc(u)
      if(u .lt.0.5 ) then
         ksubc = k0l
         echg = -1
      endif
      nusubc = -echg
      call cmkptc(kkaon, ksubc, 0, a(1))
      call cmkptc(kelec, 0, echg, a(2))
      call cmkptc(kneue, nusubc, 0, a(3))

      call cnbdcy(3, pj%mass, a,  0, w, icon)
      np = 3
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!     *******************************
!            K0s pi+ pi0    6.8        24.62
      entry cdDecay5( pj, a, np)


      ksubc = k0s
      pichg = pj%charge
      call cmkptc(kkaon, ksubc, 0, a(1))
      call cmkptc(kpion, 0, pichg,  a(2))
      call cmkptc(kpion, 0, 0,  a(3))

      call cnbdcy(3, pj%mass, a,  0, w, icon)
      np = 3
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!     *************************************
!           K-pi+ pi+ pi0        6.0        30.62
      entry cdDecay6( pj, a, np)

      kchg = -pj%charge
      call cmkptc(kkaon, 0, kchg, a(1))
      call cmkptc(kpion, 0, -kchg, a(2))
      call cmkptc(kpion, 0, -kchg, a(3))
      call cmkptc(kpion, 0, 0, a(4))

      call cnbdcy(4, pj%mass, a,  0, w, icon)
      np = 4
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!     **********************************
!          K-  pi+ e+ neue  4.1 %  34.72  
      entry cdDecay7( pj, a, np)

      kchg =-pj%charge
      nusubc = -kchg
      call cmkptc(kkaon, 0, kchg, a(1))
      call cmkptc(kpion, 0, -kchg, a(2))
      call cmkptc(kelec, 0, -kchg, a(3))
      call cmkptc(kneue, nusubc, 0, a(4))

      call cnbdcy(4, pj%mass, a,  0, w, icon)
      np = 4
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!*****************************
!         K0s + pi+ pi+ pi-    3.02      37.74   
      entry cdDecay8(pj, a, np)

      ksubc = k0s
      pichg = pj%charge
      call cmkptc(kkaon, ksubc, 0, a(1))
      call cmkptc(kpion, 0, pichg, a(2))
      call cmkptc(kpion, 0, pichg, a(3))
      call cmkptc(kpion, 0, -pichg, a(4))


      call cnbdcy(4, pj%mass, a,  0, w, icon)
      np = 4
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!    ******************************
!           K0L + pi      1.46      39.20
      entry cdDecay9(pj, a, np)

      ksubc = k0l
      pichg = pj%charge
      call cmkptc(kkaon, ksubc, 0, a(1))
      call cmkptc(kpion, 0, pichg, a(2))
      call c2bdcy(pj, a(1), a(2))
      np =2
      return
!     ****************************
!           K0s + pi      1.46      39.20
      entry cdDecay10(pj, a, np)

      ksubc = k0s
      pichg = pj%charge
      call cmkptc(kkaon, ksubc, 0, a(1))
      call cmkptc(kpion, 0, pichg, a(2))
      call c2bdcy(pj, a(1), a(2))
      np =2
      return
!  ***************************
!                 3pi + pi0     1.14      41.79
      entry cdDecay11(pj, a, np)

      pichg = pj%charge
      call cmkptc(kpion, 0, pichg, a(1))
      call cmkptc(kpion, 0, pichg, a(2))
      call cmkptc(kpion, 0, -pichg, a(3))
      call cmkptc(kpion, 0, 0, a(4))

      call cnbdcy(4, pj%mass, a,  0, w, icon)
      np = 4
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!     **************************************
!              D0    K-  mu+  neu  3.31
      entry cdDecay20(pj, a, np)


      call cmkptc(kkaon, 0, -1, a(1))
      call cmkptc(kmuon, 1, 1, a(2))
      call cmkptc(kneumu, -1, 0, a(3))

      call cnbdcy(3, pj%mass, a,  0, w, icon)
      np = 3
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!********************************************
!     D0        K- pi+ pi-       13.9   13.9
      entry cdDecay21(pj, a, np)

      call cmkptc(kkaon, 0, -1, a(1))
      call cmkptc(kpion, 0, 1, a(2))
      call cmkptc(kpion, 0, -1, a(3))

      call cnbdcy(3, pj%mass, a,  0, w, icon)
      np = 3
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!****************************************************
!                 K-+2pi+ pi-       5.4   19.3
      entry cdDecay22( pj, a, np)


      call cmkptc(kkaon, 0, -1, a(1))
      call cmkptc(kpion, 0, 1, a(2))
      call cmkptc(kpion, 0, 1, a(3))
      call cmkptc(kpion, 0, -1, a(4))

      call cnbdcy(4, pj%mass, a,  0, w, icon)
      np = 4
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!      *********************************
!                  K-2pi+ pi- pi0    4.2   23.5
      entry cdDecay23( pj, a, np)


      call cmkptc(kkaon, 0, -1, a(1))
      call cmkptc(kpion, 0, 1, a(2))
      call cmkptc(kpion, 0, 1, a(3))
      call cmkptc(kpion, 0, -1, a(4))
      call cmkptc(kpion, 0, 0, a(5))

      call cnbdcy(5, pj%mass, a,  0, w, icon)
      np = 5
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!      ***************************
!            K- pi+ pi0       3.89   27.39
      entry cdDecay24( pj, a, np)

      call cmkptc(kkaon, 0, -1, a(1))
      call cmkptc(kpion, 0,  1, a(2))
      call cmkptc(kpion, 0,  0, a(3))
      call cnbdcy(3, pj%mass, a,  0, w, icon)
      np =3
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo

      return
!      *************************
!               K + e  + neue    3.55
      entry  cdDecay25( pj, a, np)
      if( pj%subcode == antip ) then
         kchg = -1
         echg = 1
      else
         kchg = 1
         echg = -1
      endif
      call cmkptc(kkaon, kchg, kchg, a(1))
      call cmkptc(kelec, echg, echg, a(2))
      call cmkptc(kneue, kchg, 0, a(3))

      call cnbdcy(3, pj%mass, a,  0, w, icon)
      np = 3
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
!**************************************
!                          K0s + pi+ pi-    2.99   33.96
      entry cdDecay26(pj, a, np)


      call cmkptc(kkaon, 0, -1, a(1))
      call cmkptc(kpion, 0, 1, a(2))
      call cmkptc(kpion, 0, -1, a(3))

      call cnbdcy(3, pj%mass, a,  0, w, icon)
      np = 3
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return
! *******************************
!              K0s + pi+ pi- pi0   5.2
      entry cdDecay27(pj, a, np)
      ksubc = k0s
      call cmkptc(kkaon, ksubc, 0, a(1))
      call cmkptc(kpion,  0, 1, a(2))
      call cmkptc(kpion,  0, -1, a(3))
      call cmkptc(kpion,  0, 0, a(4))

      call cnbdcy(4, pj%mass, a,  0, w, icon)
      np = 4
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
! *******************************
!           K + pi + e^+ +neue   
      entry cdDecay28(pj, a, np)
      call rndc(u)
      if( pj%subcode == antip ) then
         echg = 1
      else
         echg = -1
      endif
      if( u < 0.5 ) then
         kchg = 0
         pichg = -echg
      else
         kchg =  -echg
         pichg = 0
      endif
      call cmkptc(kkaon,  0, kchg, a(1))
      call cmkptc(kpion,  0, pichg, a(2))
      call cmkptc(kelec,  0, echg, a(3))
      call cmkptc(kneue, -echg, 0, a(4))

      call cnbdcy(4, pj%mass, a,  0, w, icon)
      np = 4
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      return

! *******************************
!         K+pi+ mu^+ + neu      0.1523         
      entry cdDecay29(pj, a, np)
      call rndc(u)
      if( pj%subcode == antip ) then
         muchg = 1
      else
         muchg = -1
      endif
      if( u < 0.5 ) then
         kchg = 0
         pichg = -muchg
      else
         kchg =  -muchg
         pichg = 0
      endif
      call cmkptc(kkaon,  0, kchg, a(1))
      call cmkptc(kpion,  0, pichg, a(2))
      call cmkptc(kmuon,  0, muchg, a(3))
      call cmkptc(kneumu, -muchg, 0, a(4))

      call cnbdcy(4, pj%mass, a,  0, w, icon)
      np = 4
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      end
