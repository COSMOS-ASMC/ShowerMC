      subroutine showXsec(IModel, proj, KE1, KE2, Estep)
      use modXsecMedia
      implicit none
#include "Zptcl.h"
#include "Zcode.h"

      character(*),intent(in):: IModel
      character(*),intent(in):: proj
      real(8),intent(in):: KE1, KE2, Estep

      real(8):: KE, xst, mfpi, mfpt, xsi, xs2, mfp

       type(ptcl)::   pj

      KE= KE1
      call setProj(proj, pj)
!!!!!!!!!
      write(0,*) ' name ', media(1)%name, NoOfMedia
!!!!!!!
      do while (KE < KE2 ) 
               ! p(1:3) will not be used
         if( pj%code == kgnuc ) then
            pj%fm%P(4) = pj%mass + KE*pj%subcode
         else
            pj%fm%P(4) = pj%mass + KE ! p(1:3) will not be used
         endif
         call cGetXsec(IModel, pj, media(1), xst, mfpt)
!         call cGetXsec2(IModel, pj, media(1), xs2, mfp)
!!         if( pj.code >= kpion .and. pj.code <= knuc ) then
!!            if( pj.fm.p(4) .lt.  4.1d0 ) then 
            call cGetXsec("other", pj, media(1), xsi, mfpi)
!!         endif

         mfpt = mfpt/10.d0   ! kg/m2 --> g/cm2
!!         mfp = mfp/10.d0
         mfpi = mfpi/10.d0

         write(*,'(1p, 5E13.4)') 
     *    KE, xst, mfpt, xsi, mfpi
          KE = KE*10.0d0**Estep
      enddo
      end      subroutine showXsec
      
      subroutine getInela(pj, md,  xso, mfp)
      use modXsecMedia 
      implicit none
#include "Zptcl.h"
       type(ptcl):: pj
      type(xsmedia),intent(inout):: md ! input/output 
      real(8),intent(out):: xso
      real(8),intent(out):: mfp

      real(8):: sumns,xs
      integer:: i

      sumns=0.
!!!!!!!!!
      write(0,*) ' md%noOfElem=', md%noOfElem
      write(0,*) ' md%mbtoPkgrm=',  md%mbtoPkgrm
!!!!!!!!!!!!

      do i =1, md%noOfElem
         call cinelx(pj, md%elem(i)%A,  md%elem(i)%Z, xs)
!!!!!!!!!
         write(0,*)' i,  md%elem(i)%A,  md%elem(i)%Z,md%elem(i)%No'
         write(0,*) i,  md%elem(i)%A,  md%elem(i)%Z,md%elem(i)%No
         write(0,*)
     *    ' xs =', xs, ' md%elem(i)%nsigma=',md%elem(i)%nsigma
!!!!!!!!!!!

         md%elem(i)%nsigma = xs*md%elem(i)%No
         sumns = sumns + md%elem(i)%nsigma
      enddo
      md%xs = sumns
      mfp = 1.0d0/( md%mbtoPkgrm * sumns)
      xso = sumns
      end subroutine getInela
      
      subroutine setProj(proj, pj)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      character(*),intent(in):: proj
       type(ptcl)::   pj
      integer:: A

      select case(proj)
          case('p') 
             call cmkptc(knuc, -1, 1, pj)
          case('n')   
             call cmkptc(knuc, -1, 0, pj)
          case('pbar')   
             call cmkptc(knuc, 1, -1, pj)
          case('nbar')   
             call cmkptc(knuc, 1, 0, pj)
          case('pi+')   
             call cmkptc(kpion, -1, 1, pj)
          case('pi-')   
             call cmkptc(kpion, 1, -1, pj)
          case('k+')   
             call cmkptc(kkaon, -1, 1, pj)
          case('k-')   
             call cmkptc(kpion,  1, -1, pj)
          case('k0')   
             call cmkptc(kpion, -1, 0, pj)
          case('A')                
             write(0,*) 'Enter mass # '
             read(*,*) A
             call cmkptc(kgnuc, A, max(int(A/2),1), pj) 
          case default 
             write(0,*) ' projectile : ', proj
             write(0,*) ' not permitted'
             stop
      end select
      end      subroutine setProj


      
