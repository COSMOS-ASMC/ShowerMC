      subroutine cinipipx
      implicit none
!         initialize leading particle x sampling for pi+p-->pi leading
!
#include "Zcinippxc.h"

      external cintdndx2

      integer  i, icon
      logical first/.true./
      save first

      real*8 s 
      real*8 x, x1, eps, ans
      real*8 xa(npitbl)
      real*8 dndx(npitbl)  ! unnormalized dn/dx of leading particle at
      data xa/0, 
     * .01, .02, .03, .04, .05, .06, .07, .08, .09, .1,
     * .11, .15, .20, .25, .30, .35, .40, .45, .50, .55,
     * .60, .65, .70, .725, .75, .77, .78, .79, .80,.82, 
     * .84, .85, .87, .88, .90, .91, .92, .93, .94, .95, 
     * .96, .97, .98, .99, 1.00/
      data dndx/0,
     *  31.,  38., 43., 48., 52., 53.5, 55.5, 58., 59., 60.,
     * 58.5,  58., 52, 44,  36,   30,  25, 21,  17, 14.5,
     * 12.5, 11, 9.6, 9.0, 8.8, 8.6, 8.5, 8.4, 8.4, 8.4,
     * 8.6, 8.8, 9.2, 9.5, 10.1, 11.5, 12.2, 14.2, 16., 21.,
     * 26.0, 33., 50, 80, 130./

       eps = 1.d-4
       if(first) then
!          move xa into xval; 
          do i = 1, npitbl
             xval(i) = xa(i)
          enddo
          first = .false.
       endif

!          integral of  dn/dx from 0 to 1. 
       call ktrpzIntT2(dndx, 1, npitbl, xa, 1,  0.d0, 1.d0, s)
!
!         make normalized table
!
       do i = 1, npitbl
          ndndx(i) = dndx(i)/s
       enddo
!         make table of intgral(0:x) of ndndx. for x = 0 to 1.0 step
!         0.01
        x = 0.
        dx = 0.01d0
        do i=2, n-1
           x = x + dx
           call ktrpzIntT2(ndndx, 1, npitbl, xa, 1, 0.d0,
     *     x, intendndx2(i) )
        enddo
        intendndx2(1) = 0.
        intendndx2(n) = 1.
!           solve  inte(0:x) of ndndx = u for u = 0 to 1.0 step 0.01
        x = dx
        uconst = 0.
        x1 = 0.
        do i = 2, n-1
           uconst = uconst + dx
           call kbinChop(cintdndx2, x1, 1.d0, x, eps, ans, icon)
           if(icon .ne. 0) then
              call cerrorMsg(
     *       'failed in making pi-leading particle sampling table', 0)
           endif
           pipsx(i) = ans        
           x = ans
           x1 = x
        enddo
        pipsx(1) = 0.
        pipsx(n) = 1.
      end
      real*8 function cintdndx2(x)
      implicit none
      real*8 x

#include "Zcinippxc.h"

      real*8 ans

      call ktrpzIntT2(ndndx, 1, npitbl, xval, 1,  0.d0,  x, ans)
      cintdndx2 = ans - uconst
      end

