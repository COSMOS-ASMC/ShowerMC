      subroutine howmuch(limit, E0, NN, cosz)
      implicit none
#include  "Zmaxdef.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Ztrackp.h"
#include "Zincidentp.h"
#include "Zprivate.h"
!  #include "Zprivate2.h"
                  
      real  E0  ! input primary total energy in GeV
      integer NN ! input primary nucleon number
      real  cosz ! input.  zenith in cos
      real  limit(4) ! input. approx max. number of ptcls
                  !  to be recorded in each area bin

      real depth, depthc, sum, Nx, age, r, depc
      integer idep, ir, j
      integer code
      real*8 cvh2thick


      if( .not. FreeC ) then
         depc = cvh2thick(HeightOfInj)*0.1 -55.4
      else
         depc = 0.
      endif

      do idep = 1, ansites
         j = indivdep(idep)
         depth = DepthList(j)*0.1  !  g/cm2
!
!            for very deeply penetrated shower,
!            depth correction needed.
         depthc =max(depth - depc, 50.0)
         do code = 1, 4
            write(0,*) ' code=',code,' depth=',depth,' g/cm^2', 
     *      ' depthc=',depthc
            call crecprob(depthc, code, limit(code),
     *      sngl(dfai), E0, NN,
     *      cosz,   nrbin,  recprob(1, code,  idep),
     *           nptcls(1, code, idep),  age,  sum, Nx)

            do ir= 1, nrbin
               if( recprob(ir, code,  idep) .gt. 1.e38) then
                   recprob(ir, code,  idep) = 1.e38
                endif
               r=0.01*10**((ir-1)*0.1)
               write(0, *) r, age, recprob(ir, code,  idep),
     *         nptcls(ir, code,  idep)
            enddo
         enddo
      enddo
      end
