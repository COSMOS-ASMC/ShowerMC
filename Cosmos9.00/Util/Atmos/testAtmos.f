#include "Zcondc.h"
      program main
      implicit none
      real(8),external:: cvh2thick,  cvh2temp, cvh2den,
     *     cvh2scaleh, cvh2denp, cvh2den2p, cthick2h, 
     *     cthick2den
      real(8):: height, depth, rho, temp, rhop, rhodp, sh
      real(8):: dep2h, dep2rho
      integer::model
      call creadParam(5)
      call ciniAtmos
      

#if ATMOSPHERE == 3
      call cNRLHeaderW0(6)
#else      
      call cqAtmosModel(model)
      write(*,*) '# Atmosphere model=', model
#endif
      
      write(*,'(a,a)')
     *   "#  H(m)   depth(g/cm2) T(K)  rho(kg/m3) SH(km) ",
     *   " rho'(SI)     rho''(SI)  dep->h  dep->rho"
      write(*,'(a)')
     *  "# ------------------------------------------------------"

      height = -500

      do  while( height < 500.d3)
         depth= cvh2thick( height )
         temp= cvh2temp( height )
         rho = cvh2den( height )
         sh = cvh2scaleh( height )
         rhop = cvh2denp( height )
         rhodp = cvh2den2p( height )
         dep2h = cthick2h( depth )
         dep2rho = cthick2den( depth )

         write(*,'(1p, E10.3, 1p, E11.3, 0p, f7.1, 1p, E10.3,
     *    -3p, f7.2, 1p, E12.3, E11.3, 1p, E10.2, E10.3)' )
     *     height, depth/10., temp, rho, sh, rhop, rhodp,
     *     dep2h, dep2rho
         if( height < 10005.d0 ) then
            height = height + 100.d0
         else
            height = height*10.d0**0.05d0
         endif
      enddo

      end

