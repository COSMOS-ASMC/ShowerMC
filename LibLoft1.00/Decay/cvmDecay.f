!      ******************************
!        rho 0, +, - decay
       subroutine crhodc(vm, a, np)
!      ******************************
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
!       record /ptcl/ vm, a(*)
       type(ptcl),intent(in):: vm
       type(ptcl),intent(out):: a(*)

       integer,intent(out):: np

       integer::  icharge
       real(8):: u

       call rndc(u)
       if(vm%charge .eq. 0) then
          if(u < 4.55d-5) then
!               mu+ mu- channel
             call cmkptc(kmuon, 0, 1, a(1))
             call cmkptc(kmuon, 0, -1, a(2))
          else
!               pi+ pi-
             call cmkptc(kpion, 0, 1, a(1))
             call cmkptc(kpion, 0, -1, a(2))
          endif
       else
          !      rho+/-
          call cmkptc(kpion, 0, 0, a(1))
	  icharge = vm%charge
          call cmkptc(kpion, 0, icharge, a(2))
       endif
       call c2bdcy(vm, a(1), a(2))
       np = 2
       end
!      ******************************
       subroutine comgdc(vm, a, np)
!      ******************************
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
!       record /ptcl/ vm, a(*)
       type(ptcl),intent(in):: vm
       type(ptcl),intent(out):: a(*)
       
       integer,intent(out):: np
       

       real(8):: u
       integer:: i
       
!         omega--->pi+ pi- pi0   : 88.8 %-->89.2 (2016)
!                  pi0 gamma:       8.5  -->8.28 (//)  97.48
!                  pi+ pi- :        2.7  -->1.53 (//)  99.01
!                  mu+ mu- :        9.0e-5             99.01009
!                  pi0 mu+ mu-:    1.3e-4              99.01022
          call rndc(u)
          if(u .lt. .892d0) then
              call c3pidc(vm, a, np)
          elseif(u .lt. .9748d0) then
             call cmkptc(kpion, 0, 0, a(1))
             call cmkptc(kphoton, 0, 0, a(2))
             call c2bdcy(vm, a(1), a(2))
             np = 2
          elseif( u < 0.9901d0) then
!               pi+ pi-
             call cmkptc(kpion, 0, 1, a(1))
             call cmkptc(kpion, 0,-1, a(2))
             call c2bdcy(vm, a(1), a(2))
             np = 2
          elseif( u < 0.9901009d0) then
!               mu+ mu-
             call cmkptc(kmuon, 0, 1, a(1))
             call cmkptc(kmuon, 0, -1, a(2))
             call c2bdcy(vm, a(1), a(2))
             np = 2
          elseif( u <  0.9901022d0) then
!               mu+ mu- pi0
             call cmkptc(kmuon, 0, 1, a(1))
             call cmkptc(kmuon, 0, -1, a(2))
             call cmkptc(kpion, 0, 01, a(3))
             np = 3
             do i = 1, np
                call cibst1(i, vm, a(i), a(i))
             enddo
          else
              !  rest: regards 3pi
             call c3pidc(vm, a, np)
          endif
        end

!    * ****************************      
        subroutine c3pidc(vm, a, np)
!    * ****************************      
        implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
        type(ptcl),intent(in):: vm
        type(ptcl),intent(out):: a(*)
!        record /ptcl/ vm, a(*)
        integer,intent(out):: np

        integer i
        integer icon, ntry
        real(8):: w
        character*70  msg

        call cmkptc(kpion, 0, 1, a(1))
        call cmkptc(kpion, 0, 0, a(2))
        call cmkptc(kpion, 0, -1, a(3))
!             pi+ pi- pi0
        ntry = 0
        icon = 1
        do while( icon .ne. 0 .and. ntry < 200)
           call cnbdcy(3, vm%mass, a,  0, w, icon)
           ntry = ntry + 1
        enddo
        if(icon .ne. 0) then
           write(0, *) ' cnbdcy fails for omega 3 pi',
     *          'cms=',vm%mass, 'icon=',icon
           np=0
        else
           np=3
           do i = 1, np
              call cibst1(1, vm, a(i), a(i))
           enddo
!////////////
!           write(*,*) ' omega 3body '
!//////////////////
        endif
       end
       subroutine cphidc(vm, a, np)
        implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
        
!     record /ptcl/ vm, a(*)
        type(ptcl),intent(in):: vm
        type(ptcl),intent(out):: a(*)
        
        integer,intent(out):: np
        
        real(8)::  u
        integer ic, nx
!        record /ptcl/pw
        type(ptcl):: pw
        nx = 0
!            phai-->k+ k-      49.5 %--> 48.9 (2016)     %
!                   k0l k0s 34.4     --> 34.2 (//)       83.1 %
!                   rho pi      12.9 (rho+ pi-, rho- pi+, rho0 pi0)
!                 pi+ pi- pi0      1.9  %
!     sum of above two 15.32(2016) : so
!     we use               12.9-->12.9   1.9->2.42     96.0  %
!                                                      98.42 %
!                 eta + gamma            1.309 %       99.729 
!     mu+ mu-       2.87e-4     0.99729
!                              +0.000287 =   0.997577

          call rndc(u)
          if(u .lt. 0.495d0) then
             call cmkptc(kkaon, 0, 1, a(1))
             call cmkptc(kkaon, 0, -1, a(2))
             call c2bdcy(vm, a(1), a(2))
             np = 2
          elseif(u .lt. 0.831d0) then
             call cmkptc(kkaon, k0l, 0, a(1))
             call cmkptc(kkaon, k0s, 0, a(2))
             call c2bdcy(vm, a(1), a(2))
             np = 2
          elseif(u .lt. 0.960d0) then
             !      rho + pi
             call rndc(u)
             ic=int(3.*u)-1
             call cmkptc(kpion, 0, ic, a(1))
             call cmkptc(krho, 0, -ic, a(2))
             call c2bdcy(vm, a(1), a(2))
!                  rho is made to decay
             pw = a(2)
             call crhodc(pw, a(2), nx)
             np=nx+1
          elseif( u < 0.9842d0) then
             call c3pidc(vm, a, np)
          elseif( u < 0.99729d0) then
             !  eta gamma
             call cmkptc(keta, 0, 0, a(1))
             call cmkptc(kphoton, 0, 0, a(2))
             call c2bdcy(vm, a(1), a(2))
             np =2
          elseif( u <   0.997577d0 ) then
                ! mu+ mu-
             call cmkptc(kmuon, 0, 1, a(1))
             call cmkptc(kmuon, 0, -1, a(2))
             call c2bdcy(vm, a(1), a(2))
             np =2
          else
!  rest is assigned to main br
             call cmkptc(kkaon, 0, 1, a(1))
             call cmkptc(kkaon, 0, -1, a(2))
             call c2bdcy(vm, a(1), a(2))
             np = 2
          endif
        end
