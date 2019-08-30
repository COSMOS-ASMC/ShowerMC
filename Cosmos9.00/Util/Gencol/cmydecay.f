!  activate next to see how the charged ptcl # changes by decay
!  #define SEEDETAIL
      subroutine cmydecay(Jdecay, tau,  a, nin, nout)
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
#include "Zprivate.h"

      integer,intent(in)::Jdecay(klast)  ! if Jdecay(i) !=  0
         !  the ptcl is decayed, if i is one of
         !  kbomega, kdmes, keta, kgzai,  klambda, klambdac, kpion, ksigma
         !  (kpion --> pi0)
      real(8),intent(in)::tau ! life time (sec), above which the partcles
                 ! are  regarded as stable and not made to decay.
                 ! Below this, decay is forced. For the decay product,
                 ! the same is applied.
            ! IF tau =0, random sampling is tried if the
            !  particle decay or not using correct life time and 
            !  distance dist. (see below)
            ! if decay   takes place, decay product is put 
      integer,intent(in):: nin  ! # of ptcls in a
      type (ptcl):: a(*)  ! input/output.  size of a must be large
                            ! say, > (nout - nin) *3 + nin
      integer,intent(out):: nout   !  ptcls in a. nout >= nin.


      type (ptcl):: b(maxn)  ! working
      integer i, j, code, m, n
      real(8):: pol, ctau
      real(8),save::dist=141.2 ! m from col. to LHCf
      real(8) u, gamma, decayl
      logical decay
!//////
#ifdef  SEEDETAIL
      integer nchgin, nchgout
      nchgin = 0
      do i = 1, nin
         if(a(i)%charge /= 0) nchgin= nchgin+1
      enddo
#endif
!//////////////
      n = 0
      do i = 1, nin
         code = a(i)%code
         j = 0
         if( Jdecay(code) /=  0 .or. tau == 0. ) then
            call cgetctau(a(i), ctau)
            if(tau == 0. ) then
!                random sample for decay
               gamma = a(i)%fm%p(4)/a(i)%mass
               call rndc(u) 
               decayl = -log(u)*ctau*gamma
               decay = (decayl < dist )
            else
               decay = (ctau/3.e8 < tau)
            endif
            if(.not. decay) then
               n = n + 1
               b(n) = a(i)
            else
               if( code == kpion ) then
                  if( a(i)%charge == 0) then
                     call cpi0Decay( a(i),  b(n+1), j)
                  endif
               elseif( code == kkaon ) then
                  if( a(i)%subcode == k0s ) then
                     call ckaonDecay(a(i), .false.,  b(n+1), j, pol)
                  endif
               elseif( code .eq. kdmes ) then
                  call cdDecay( a(i), b(n+1), j)
               elseif( code .eq. keta ) then
                  call cetaDecay( a(i), b(n+1), j)
               elseif( code .eq. kgzai ) then
                  call cgzaiDecay( a(i), b(n+1), j ) 
               elseif( code .eq. klambda ) then
                  call clambdaDcy( a(i), b(n+1), j ) 
               elseif( code .eq. klambdac ) then
                  call clambdacDcy( a(i), b(n+1), j ) 
               elseif( code .eq. ksigma ) then
                  call csigmaDecay( a(i), b(n+1), j ) 
               elseif( code .eq. kbomega ) then
                  call cbomegaDcy( a(i), b(n+1), j ) 
               else
                  write(*,*)  ' code =', code
                  call cerrorMsg('cmydecay error', 0)
               endif
               if(j > 0 ) then
!/////////////////////
#ifdef SEEDETAIL
                  write(0,*) ' j=',j, ' code=',code, ' chg=',a(i)%charge
#endif
!///////////////
                  n = n + j
               endif
            endif
         else
            n = n + 1
            b(n) = a(i)
         endif
      enddo
!/////////
#ifdef SEEDETAIL
      nchgout = 0
      do i = 1, n
         if(b(i)%charge /= 0) nchgout = nchgout+1
      enddo
      write(0,*) ' chgin out', nchgin, nchgout
#endif
!/////////////////
      do i = 1, n
         a(i) = b(i)
      enddo
      nout = n
      end
