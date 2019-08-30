!  Usage:  compile by make -f cgetXs.mk
!      verbose mode:  (may be target is simple)
!      ./a.out > result
!      file input mode:(may be target is complex)
!                  prepare a file like
!          4            /  1-> proton 2->pion 3->K  4->A
!          56 26        /  if 4, give A,Z of projectile else omit this line
!          14.0, 16.0, 40.0 /     A's for target. "/" is mandatory
!          1.56 .42  .01
!       and
!      ./a.out < thefile > result
!   **NOTE*** verbose mode and file input mode may be
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

      real(8),parameter:: dE=0.01d0   !  log 10 step
      real(8),save:: E1= 0   ! starting Etotal( =>mass*1.1 )
      real(8),parameter:: E0max=1e12 ! total E max
      integer,save:: nebin=480  ! # of  Energy bins 


      integer,save:: pjA, pjZ
      integer,save:: targetnucchg=1  
      contains
      subroutine computeXS(pj, xsa)
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
      type(ptcl):: pj  !  projectile . E is fixed here
      real(8),intent(out):: xsa(nebin,3)

      type(ptcl):: tg  ! input target particle in lab.
      integer:: i, j
      real(8):: E, sum(3), xs(3) 
      integer:: TA
      integer:: code
      real(8):: temp
      real(8),external::cA2Z

      E = E1  


      do i = 1, nebin
         pj%fm%p(4) = E*pjA

         sum(:) = 0.
         do  j = 1, nel
            if( A(j) > 1 ) then
                                       ! charge is not used so rough
               temp = A(j)
               call cmkptc(kgnuc, A(j), int(cA2Z(temp)), tg)
            else
               call cmkptc(knuc, -1, targetnucchg, tg) 
            endif
            call cgetXsInterface2(pj, tg, xs)
            sum(:) = sum(:) + xs(:)*portion(j)
         enddo
         xsa(i,:) =  sum(:)
         E = E*10.0d0**dE
         if( E*pjA > E0max ) exit
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
      real(8),allocatable:: xs( :,:)
      real(8):: E, temp, p, KE
      integer:: nela(1)
      integer,save::  kpj, pjcharge, pjsubcode, charge

      
      pjA = 1
      
      write(0,*) 'Enter projectile '
      write(0,*)
     *  ' p/n/p~/n~->1; pi->2; K->3; A->4;g->5; 0->stop'
      read(*,*)  kpj

      if( kpj >= 1 .and. kpj <= 5) then
         if(kpj ==  4) then
            write(0,*) ' Enter proj. A, Z (intger)'
            read(*,*)  pjA, pjZ
         elseif( kpj /= 5 ) then
            write(0,*)
     *      ' Enter proj. subcode & charge: / ==> Default=(-1 1)'
            write(0,*) ' subcode for ptcl=', regptcl
            write(0,*) ' subcode for anti-ptcl=', antip
            pjsubcode= -1
            pjcharge = 1
            read(*,*) pjsubcode, pjcharge
            if( kpj == 3 .and.  pjcharge==0 ) then
               pjsubcode = 4  ! k0short (5-->long)
            endif   
         endif
      elseif( kpj == 0) then
         stop
      else
         write(0,*) ' input=',kpj, ' invalid'
         stop
      endif
      write(0,*)
     * "Enter target A (or A's for compound/mixture target) with",
     * "  / at last pos. (# of A's  <= 15)"
      A(:) = 0.
      read(*,*) A(:)
      nela = minloc(A(:))
      nel = nela(1) -1
      if( nel == 0 ) stop
      if( nel == 1 ) then
         portion(nel) = 1.
         if(A(1) == 1.) then
            write(0,*) "Enter charge of target nucleon"
            read(*,*)  targetnucchg
         endif
      else
         write(0,*) "Enter relative portion of A's"
         read(*,*) portion(1:nel)
      endif

      temp = sum(portion(1:nel))
      portion(:) = portion(:)/temp

      call cgetXsIni
      ! make proj.
      if( kpj == 1 ) then
         call cmkptc(knuc, pjsubcode, pjcharge, pj)
      elseif( kpj == 2 ) then
         call cmkptc(kpion, pjsubcode, pjcharge, pj)
      elseif( kpj == 3 ) then
         call cmkptc(kkaon, pjsubcode, pjcharge, pj)
      elseif( kpj == 4) then
         call cmkptc(kgnuc, pjA, pjZ, pj)
      elseif( kpj == 5 ) then
         call cmkptc(kphoton, -1, 0, pj)
      else
         stop 'strange '
      endif

      if(kpj /= 5 .and.  E1 == 0.) then
         E1=pj%mass*1.1
      else
         E1= 155.e-3  ! 153MeV is threshold
      endif


      
      nebin = log10( E0max/E1 )/dE
      allocate( xs(nebin,3) )

      call computeXS(pj, xs)
!         output

      write(*,'(a,i2,i4,i3)') '# pj=',pj%code, pj%subcode, pj%charge
      write(*,'(a, 15i6)')    '#  target=', A(1:nel)
      write(*,'(a, 15f6.3))') '# portion=', portion(1:nel)
      if( kpj == 4 .and. ( nel>1 .or.  A(1)> 1.) ) then
         write(*,'(a)')
     *   '# Et(GeV)  Sinel(mb) p(GeV)  KE(Gev)'
      else
         write(*,'(a)') 
     *   '# Et(GeV)  Sinel  Stot Sel (mb) p(GeV) KE(GeV)'
      endif
      E= E1
      do i = 1, nebin
         p = sqrt( (E*pjA)**2 - pj%mass**2 )
         KE = E*pjA- pj%mass
         if( kpj == 4  .and. (nel > 1 .or. A(1) > 1.0))  then
            write(*,'(1p, 4g13.4 )')  E*pjA, xs(i,1), p, KE
     *
         else
            write(*,'(1p, 6g13.4 )')
     *         E*pjA, xs(i,1), xs(i,2), xs(i,3), p, KE
         endif
         E = E*10.d0** dE
         if(E*pjA> E0max) exit
      enddo
      end program main
