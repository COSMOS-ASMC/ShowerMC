      subroutine csampAF0(iowk, filen, sampInfo)
      implicit none
#include "ZsampAF.h"
      integer iowk  ! input file logical number temporarily used
      character*(*)  filen  ! input. file name which contains (x,dn/dx)

      type(sampAF):: sampInfo

      integer icon
      integer i
      call copenf(iowk, filen, icon)
      if(icon .ne. 0 ) then
         write(0,*) ' error '
         write(0,*) 'file: ',filen
         write(0,*) ' could not be opened'
         stop
      endif
      i = 0
      call cskipComment(iowk, icon)
      if(icon .ne. 0 ) stop
      do while ( .true. )
         read(iowk, *, end=100) sampInfo.x(i+1), sampInfo.y(i+1)
         i  = i + 1
      enddo
 100  continue
      close(iowk)
      sampInfo.n =  i
      call ksampAF0(sampInfo.x, sampInfo.y, sampInfo.n, 
     *    sampInfo.coef, sampInfo.n, sampInfo.yi, 
     *    sampInfo.sum,  sampInfo.coef2)
      end
      subroutine csampAF(sampInfo, xs)
      implicit none
#include "ZsampAF.h"
      type(sampAF):: sampInfo
      real*8 xs
      call ksampAF(sampInfo.x, sampInfo.yi, sampInfo.n,
     *             sampInfo.coef2, sampInfo.n, xs)
      end
      subroutine csampAFIntp(sampInfo, xv, ans)
!         this is not for sampling but simply
!       get value of y at xf

      implicit none
#include "ZsampAF.h"
      type(sampAF):: sampInfo ! input obtained by csampAF0
      real*8 xv  ! input
      real*8 ans  ! output y at xv
      call  kcsplIntp(sampInfo.x, sampInfo.y, sampInfo.n,
     *   sampInfo.coef, sampInfo.n, xv, ans)
      end
      subroutine csampAFmax(sampInfo,  xmax, fmax)
!       find max position and value of given function
!      (approx value)
      implicit none
#include "ZsampAF.h"
      type(sampAF):: sampInfo ! input obtained by csampAF0
      real*8 xmax ! ouput. max position  in (x1,x2) ; approx value
      real*8 fmax ! outpu. max function value

      real*8 x1, x2
      real*8 x, dx, temp
      integer i

      x1 = sampInfo.x(1)
      x2 = sampInfo.x(sampInfo.n)
      x = x1
      dx = (x2-x1)/sampInfo.n/10.
      xmax = x1
      call csampAFIntp(sampInfo, xmax, fmax)
      call csampAFIntp(sampInfo, x2, temp)
      if( fmax .lt. temp ) then
         xmax = x2
         fmax = temp
      endif
      x = x + dx
      do while (x .lt. x2-dx/2 )
         call csampAFIntp(sampInfo, x, temp)
         if( fmax .lt. temp ) then
            xmax = x
            fmax = temp
         endif
         x = x + dx
      enddo
      end
