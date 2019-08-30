      subroutine cDsDecay(pj,  a,  np)
!   Ds+/-  decay:
!     we consider decays with BR> 1% and  with neu.
!      no polarization is considered.
!     Ds+ -->  τ + neutau        5.48 
!              φ + e+ + neue            2.39 %
!            η +  e+ + neue+ η' +  e+ + neue  2. % this included
!                                                  in next  
!            η + e+ + neue     2.29 + 2.79 = 5.08 %

!  (PDG not understandable; e+ semileptonic means e+ neue+pi?)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"


       type(ptcl),intent(in):: pj         ! input. Ds
       type(ptcl),intent(out):: a(*) ! output. produced ptcls
       integer,intent(out):: np    ! output. no. of ptcls produced   
!
       integer,parameter:: nbr=4  ! # of branches
!    br:  ratio (%), the last one may be any. total is adjuested
!     to be  100 %       
       real(8),save::br(nbr)=(/5.48, 2.39, 5.08, 0./)
       real(8),save:: cbr(nbr)
       logical,save:: first = .true.
       integer:: i, icon, charge, subcode
       real(8):: u
       real(8),save:: w=0.
       integer,parameter:: jw=0
       
       integer,intent(out):: brinfo

       integer,save:: branch
       
       if( first )  then
          cbr(1) = br(1)
          do i = 2, nbr
             cbr(i) = cbr(i-1) + br(i)
          enddo
          cbr(:) = cbr(:)/100.
          cbr(nbr) = 1.0
          first = .false.
       endif

       if( pj%charge == 0 ) then
          write(0,*) 'cDsDecay is called for charge 0 '
          write(0,*)' code, subcode,charge=',
     *         pj%code, pj%subcode, pj%charge
          stop
       endif
       
       call rndc(u)
!        electron(or tau) charge 
       charge = pj%charge
!          neu subcode
       subcode =  -charge
       if( u < cbr(1) ) then
         !      5.48 % 
         !       neutau  
          call cmkptc(kneutau, subcode, 0, a(1))
         !       tau
          call cmkptc(ktau,  charge, charge, a(2))
          np = 2
          branch= 1

       elseif( u < cbr(2) ) then
!          φ + e+ + neue            2.39 %
!           neue or newebar
          call cmkptc( kneue, subcode, 0, a(1) )
!           e+ or e-
          call cmkptc( kelec, charge, charge, a(2) )
!     phi
          call cmkptc( kphi, 0,    0,   a(3) )
          np = 3
          branch = 2

       elseif( u< cbr(3)) then
!     η + e+ + neue     2.29 + 2.79 = 5.08 %
!         neue
          call cmkptc( kneue, subcode, 0, a(1) )
!     e+
          call cmkptc(kelec, charge, charge, a(2) )
!     eta
          call cmkptc(keta, 0,  0,  a(3))
          np = 3
          branch = 3

       else
!          neglect branches with no neu          
          np = 0
          branch = 4
          
       endif
       if( np > 0 ) then
          if( np == 2 ) then
             call c2bdcy( pj, a(1), a(2) )
             ! already boosted
          else
             call cnbdcy(np, pj%mass, a, jw,  w, icon)
             if( icon /= 0 ) then
                write(0,*) ' cDsDecay--> cnbdcy--> icon=',icon
                write(0,*) ' mass =', pj%mass, ' np=', np,
     *               ' branch =', branch
                write(0,*) ' a(1:np)%mass=', a(1:np)%mass
                stop
             endif
          !              boost to pj system
             do i = 1, np
                call cibst1(i, pj, a(i), a(i))
             enddo
          endif
      endif
      return
!      
      entry cDsDecayBr(brinfo)
      brinfo = branch
      end


!
