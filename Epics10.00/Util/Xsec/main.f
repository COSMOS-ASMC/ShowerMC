      program main 
      use modXsecMedia
      implicit none
#include "Zptcl.h"
#include "Zcode.h"

      character(len=12):: IModel="dpmjet3"

      real(8):: KE1, KE2, Estep=0.025
      real(8):: KE, xst, xsi,  mfpt, mfpi 
      character(len=4):: proj

       type(ptcl)::   pj

      call init(proj, KE1,KE2)
      call setProj(proj, pj)
      call showInfo(media(1), proj, pj)
      KE= KE1
      do 
               ! p(1:3) will not be used
         if( pj%code == kgnuc ) then
            pj%fm%P(4) = pj%mass + KE*pj%subcode
         else
            pj%fm%P(4) = pj%mass + KE ! p(1:3) will not be used
         endif
                !  in some model, xst is Total xs (mb)
         call cGetXsec(IModel, pj, media(1), xst, mfpt)
                !  get inelastic xs 
         call cGetXsec("other", pj, media(1), xsi, mfpi)
         mfpt = mfpt/10.d0   ! kg/m2 --> g/cm2
         mfpi = mfpi/10.d0
         write(*,'(1p, 5E11.3)') 
     *    KE, xst, mfpt, xsi, mfpi
         if(KE >= KE2) exit
         KE = KE*10.0d0**Estep
      enddo
      end program main
      
      subroutine showInfo(md, proj, pj)
      use modXsecMedia 
      implicit none
#include "Zptcl.h"

      type(xsmedia),intent(in):: md ! input
      character(*),intent(in):: proj
       type(ptcl):: pj   ! input

      integer:: i

      write(*,'(a, a, a, 3i4)')
     *  '# projectile=',proj, ' code,sub chg=', 
     *  pj%code, pj%subcode, pj%charge
      write(*,'(a,a)') '# target media ', md%name
      write(*,'(a,i3)') '# No of Elem=',  md%noOfElem
!      write(0,*) '# md%mbtoPkgrm=', md%mbtoPkgrm
      write(*,'(a)') '#  i   A    Z    No'
      do i =1, md%noOfElem
         write(*,'(a, i2,f7.2, f4.0, f7.3)')
     *   "# ", i, md%elem(i)%A, md%elem(i)%Z,md%elem(i)%No
      enddo
      write(*,*)
      write(*,*) '# In the table below '
      write(*,*) '# xs1 mfp1 may be Total xsec: depends on model'
      write(*,*) '# upto some energy: say, ~3GeV for dpmjet3'
      write(*,*) '# above that energy, it should be the same as xs2'
      write(*,*) '# '
      write(*,*) '# xs2 mfp2 is always inelastic xsec'
      write(*,*) '# *** important***'
      write(*,*) '# If media is molecule or mixture, xs(=s) may not'
      write(*,*) '# be correct but mfp(=L) should be correct: why?'
      write(*,*) '# If the number of media in /g is N,'
      write(*,*) '# N*L*s=1 but there is sometimes uncertanity of '
      write(*,*) '# N in the Base media file and what we do is'
      write(*,*) '# equivalent to showing  s/a as xsec,'
      write(*,*) '# where aN*L*s/a=1, a is some unknown number'
      write(*,*) '# '
      write(*,*)
     *  '#  KE(GeV)   xs1(mb)   mfp1(g/cm2)   xs2(mb)  mpf2(g/cm2)'
      write(*,*)'# -----------------------------------------------'
      end subroutine showInfo
      
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
             call cmkptc(kkaon,  1, -1, pj)
          case('k0')   
             call cmkptc(kkaon, -1, 0, pj)
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


      
