!  Usage:  compile by make -f cgetXs.mk
!      verbose mode:  (may be target is simple)
!      ./a.out > result
!      file input mode:(may be target is complex)
!                  prepare a file like
!          4            /  1-> proton 2->pion 3->K  4->A
!          56 26        /  if 4, give A,Z of projectile else omit this line
!          14.0, 16.0, 40.0 /     A's for target. / is mandatory
!          78.1, 20.95, 0.94 /    portion of A's / is mandatory
!       and
!      ./a.out < thefile > result
!   **NOTE*** verbose mode / file input mode may be
!        used independently of the complexity of the target.
!
      module modgetXs
      implicit none
      integer,parameter::maxel=15   ! max # of elements
      integer::nel  !  actual @ of elements
!      real(8):: A(maxel)
      integer:: A(maxel)
      real(8):: portion(maxel)
      integer,save:: pjk
!      integer,parameter:: nel=3   ! # of elements
!      real(8):: A(nel)=(/14.0, 16.0, 40.0/)   ! mass # (integer value is used)
!      real(8):: portion(nel)=(/78.1, 20.95, 0.94/)    ! relative # of A's in unit volume

      real(8),parameter:: dE=0.025d0   ! 0.05 log 10 step
      integer,parameter:: nebin=400  ! # of  Energy bins (from  10^11 eV to 10^21 eV) log10
!                                     ! step 0.1
      real(8),parameter:: E1= 100.    ! from 100GeV

      contains
      subroutine computeXS(pj, xsa)
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
      type(ptcl):: pj  !  projectile . E is fixed here
      real(8),intent(out):: xsa(nebin)

      type(ptcl):: tg  ! input target particle in lab.
      integer:: i, j
      real(8):: E, sum, xs 
      integer TA
      E = E1
      do i = 1, nebin
         pj.fm.p(4) = E
         sum = 0.
         do  j = 1, nel
                                       ! charge is not use so rough
            call cmkptc(kgnuc, A(j), max(A(j)/2,1), tg)
            call cgetXsInterface(pj, tg, xs)
            sum = sum + xs*portion(j)
         enddo
         xsa(i) =  sum
         E = E*10.0d0**dE
      enddo
      end   subroutine computeXS
      end       module modgetXs


      program main
      use   modgetXs
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
      type(ptcl):: pj

      integer:: i, j
      real(8):: xs( nebin)
      real(8):: E, temp
      integer:: nela(1)
      integer:: pjA, pjZ, kpj
      
      write(0,*) 'Enter projectile '
      write(0,*) ' p->1; pi->2; K->3; A->4; 0->stop'
      read(*,*)  kpj

      if( kpj >= 1 .and. kpj <= 4) then
         if(kpj ==  4) then
            write(0,*) ' Enter proj. A, Z (intger)'
            read(*,*)  pjA, pjZ
         endif
      elseif( kpj == 0) then
         stop
      else
         write(0,*) ' input=',kpj, ' invalid'
         stop
      endif
      write(0,*)
     * "Enter target A (or A's for compound/mixture target) with",
     * "  / at last (<=15)"
      A(:) = 0.
      read(*,*) A(:)
      nela = minloc(A(:))
      nel = nela(1) -1
      if( nel == 0 ) stop
      if( nel == 1 ) then
         portion(nel) = 1.
      else
         write(0,*) "Enter relative portion of A's"
         read(*,*) portion(1:nel)
      endif

      temp = sum(portion(1:nel))
      portion(:) = portion(:)/temp

      call cgetXsIni
      ! make proton
      if( kpj == 1 ) then
         call cmkptc(knuc, -1, 1, pj)
      elseif( kpj == 2 ) then
         call cmkptc(kpion, -1, 1, pj)
      elseif( kpj == 3 ) then
         call cmkptc(kkaon, -1, 1, pj)
      elseif( kpj == 4) then
         call cmkptc(kgnuc, pjA, pjZ, pj)
      else
         stop 'strage '
      endif

      call computeXS(pj, xs)
!         output
      E= E1
      write(*,'(a,i2,i4,i3)') '# pj=',pj.code, pj.subcode, pj.charge
      write(*,'(a, 15(i4,f6.3))') '# tg=',(A(i), portion(i), i=1,nel)
      do i = 1, nebin
         write(*,'(1p, 4g13.4 )')  E, xs(i)
         E = E*10.d0** dE
      enddo
      end program main
