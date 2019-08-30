      subroutine cdecay(tauc, a, nin, nout)
#include  "Zptcl.h"
#include  "Zlife.h"
      real*8  tauc          !  life time in sec.  If a particle  life time is 
                            !  < tauc, it is made  to decay.
      type(ptcl):: a(*)   !  input/output particles. 

      integer  nin          !   input.  number of particles in a
      integer  nout         !   output.  number of particles in a.
!                             if partilces in a have life time < tauc, 
!                     they are made to decay and again stored in a.
      type(ptcl):: tmp1(10), tmp2(10) ! temporary array

      integer i, j1, j2, np
      
      do i  = 1, nin
         tmp1(1) = a(i)
         j1 = 1
         j2 = 0
         do while (j1 .ne. j2)
            do j = 1, j1
               if(tmp1(j).code .eq. kpion .and. tmp1(j).charge .eq. 0) then
                  if(t0pi0 .lt. tauc) then
                     call cpi0Decay(tmp1(j), tmp2(j2+1), np)
                     j2 = j2 + np
                  endif
               elseif(tmp1(j).code .eq. kkaon 
     *               .and. tmp1(j).charge .eq. 0) then 
                  if(tmp1(j).subcode .eq. k0s) then
                     if(t0k0s .lt.  tauc) then
                        call ckShortDecay(tmp1(j), tmp2(j2 +1), np
                        j2 = j2 + 1
                     endif
                  elseif(tmp1(j).subcode .eq. k0l) then
                     if(t0k0l  .lt. tauc) then
                        call  ck
                  
                     
