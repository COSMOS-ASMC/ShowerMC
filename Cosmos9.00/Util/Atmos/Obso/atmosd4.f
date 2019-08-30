#include "BlockData/cblkGene.h"
      implicit none
#include "Zcondc.h"
#include "Zmanagerp.h"
#include "Zatmos.h"
#include "ZcosmosExt.h"
!
!            depth --> rho, height
!
      integer i

      real*8 cvh2den, cvh2denp, cvh2den2p, cvh2temp, cthick2h
      real*8 cvh2scaleh, cvh2thick, cthick2den, h

      real*8 sh
      real*8 d0, d, rho0, dep0
      
      real*8 rl  !  r.l in g/cm2
      data rho0/1.205d-3/ ! Epics Air default density g/cm3
!          d: obs. depth step/2
      parameter (rl=36.566717d0)
!         divide 1 rl into 1/nstep rl where the density is
!         taken to be constant. 
      integer ndep, nstep, detail
      parameter (detail=8, nstep=4*detail, d=rl/nstep,  ndep=nstep*25)
      real*8 h1(0:ndep), t(0:ndep), dep(0:ndep), rho(0:ndep), 
     *   rhoc(0:ndep)

      save 

#if ATMOSPHERE == 1
      call creadParam(5)
      call creadAtmosD
!
      call catmosCnst1
      call catmosCnst2
#elif ATMOSPHERE == 2
!           default
!          read segmented atmosphere data
      call creadAtmosD
!          manipulate data
      call catmosCnst1
#endif
!       
      write(0,*) 'Enter starting depth in g/cm2'
      read(*,*) d0
      dep0 = d0 
      write(0,'(a)')
     *   '#  depth(g/cm2)  height(m)   rho(g/cm3)   rho/rho0'
      do i = 0, ndep
         dep(i) = dep0
         t(i) = dep(i)*10.
         h1(i) = cthick2h(t(i))
         rho(i) = cthick2den(t(i))/1000.
         rhoc(i) = rho(i)/rho0
         write(0, '(1p5G13.5)')  dep(i), h1(i), rho(i), rhoc(i)
         dep0 = dep0 + d
      enddo
      write(*,'(a,f8.2,a,f8.2, a, i1, a)')
     *     "#  atmosphere from depth ",d0, " to ", dep(ndep), 
     *     " g/cm2 with step 2/", nstep,  "r.l"
      write(*,'(a)') '--------------------------------------------'
      do i = 0, ndep-2
         if(mod(i,2) .eq. 0) then
            if(mod(i,nstep) .ne. 0) then
               write(*,'(i3, a, f7.5, a, 1p, g12.5)') 
     *              i/2+1, " box Air*",rhoc(i+1), 
     *              " 0 = 0 / = = + = = ",
     *              ( h1(i)-h1(i+2) )*100.
            else
               if(i .eq. 0) then
                  write(*,'(i3, a, f7.5, a, 1p, g12.5)') 
     *                 i/2+1, " box Air*",rhoc(i+1), 
     *                 " 0 1 0 / 0 0 0 1.e6 1.e6 ",
     *                 ( h1(i)-h1(i+2) )*100.
               else
                  write(*,'(i3, a, f7.5, a, 1p, g12.5)') 
     *              i/2+1, " box Air*",rhoc(i+1), 
     *              " 1 = 0 / = = + = = ",
     *              ( h1(i)-h1(i+2) )*100.
               endif
            endif
         endif
      enddo
      write(*,'(a)')
     *    '------------------------------------------'
      write(*,*) 
     * " next can be used as DepthList definition in Cosmos"
      write(*,*) "HeightOfInj=", h1(0)
      write(*,'(" DepthList=")')
      write(*,'(10f9.2)')
     *    (dep(i)*10., i=nstep,ndep,nstep)
      end
