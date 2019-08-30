      subroutine cinippx
      implicit none
!         initialize leading particle x sampling.
!

#include "Zcinippxc.h"

      external cinteLdndx

      integer  i, icon
      real*8 s 
      real*8 x, x1, eps, ans
      real*8 dndx(n)  ! unnormalized dn/dx of leading particle at x= 0, 0.01,,

!                      to 1.0
          data ( dndx(i), i=  1,   n)/
     1  0,  0.178704,  0.330057,  0.439452, 0.518544,  0.576176,  
     2  0.6189, 0.651,  0.677, 0.699,  0.719,  0.737, 0.754,  0.770, 
     3  0.7855, 0.800,  0.814,  0.845,  0.869,  0.894,  0.843,  0.804,
     4  0.806, 0.9,  0.819,  0.878,  0.856,  0.831,  0.779,  0.811, 
     5  0.866, 0.831,  0.874,  0.804,  0.826,  0.781,  0.77,  0.758,
     6  0.777,  0.746,  0.738,  0.734,  0.71,  0.704,  0.74,  0.666,
     7  0.671,  0.669,  0.658,  0.653,  0.655,  0.631,  0.601,  0.583,
     8  0.568,  0.546,  0.516,  0.482,  0.485,  0.455,  0.479,  0.476,
     9  0.435,  0.419,  0.412,  0.348,  0.382,  0.347,  0.332,  0.337,
     a  0.3,  0.335,  0.314,  0.296,  0.267,  0.284,  0.277,  0.25, 
     b  0.217, 0.253,  0.249,  0.198,  0.197,  0.194,  0.209,  0.189,  
     c  0.215, 0.193,  0.208,  0.206,  0.212,  0.193,  0.199,  0.214, 
     d  0.248, 0.254,  0.307,  0.371,  0.574,  1.914,  4.043/

        dx = 0.01d0
        eps = 1.d-4

!          integral of  dn/dx from 0 to 1. 
        call ktrpzIntegT(dndx, 1, n, 0.d0, dx, 1.d0, s)
!
!         make normalized table
!
        do i = 1, n
           ndndx(i) = dndx(i)/s
        enddo
!         make table of intgral(0:x) of ndndx. for x = 0 to 1.0 step
!         0.01
        x = 0.
        do i=2, n 
           x = x + dx
           call ktrpzIntegT(ndndx, 1, n, 0.d0, dx, x, intendndx(i) )
        enddo
        intendndx(1) = 0.
        intendndx(n) = 1.
!           solve  inte(0:x) of ndndx = u for u = 0 to 1.0 step 0.01
        x = dx
        uconst = 0.
        x1 = 0.
        do i = 2, n-1
           uconst = uconst + dx
           call kbinChop(cinteLdndx, x1, 1.d0, x, eps, ans, icon)
           if(icon .ne. 0) then
              call cerrorMsg(
     *        'failed in making leading particle sampling table', 0)
           endif
           ppsx(i) = ans        
           x = ans
           x1 = x
        enddo
        ppsx(1) = 0.
        ppsx(n) = 1.
      end
      real*8 function cinteLdndx(x)
      implicit none
      real*8 x

#include "Zcinippxc.h"

      real*8 ans

      call ktrpzIntegT(ndndx, 1, n, 0.d0, dx, x, ans)
      cinteLdndx = ans - uconst
      end
      subroutine cinippxn
      implicit none
!         initialize leading particle x sampling.
! for p--> n case
!

#include "Zcinippxc.h"

      external cinteLdndxn

      integer  i, icon
      real*8 s 
      real*8 x, x1, eps, ans
      real*8 dndx(n)  ! unnormalized dn/dx of leading particle at x= 0, 0.01,,

!                      to 1.0
          data ( dndx(i), i=  1,   n)/
     1    0, 0.1161, 0.214, 0.285, 0.337, 0.374,0.402, 0.423,
     2    0.440, 0.4549, 0.467, 0.479, 0.490, 0.50, 0.510, 
     3    0.520, 0.529, 0.549, 0.564, 0.581, 0.58, 0.59, 0.59,
     4    0.59, 0.58, 0.58, 0.57, 0.564, 0.543, 0.532, 0.521,
     5    0.510, 0.505, 0.50, 0.50, 0.50, 0.50, 0.49, 0.48,
     6    0.47, 0.46, 0.45, 0.45, 0.44, 0.44,0.43, 0.41, 0.397,
     7    0.39, 0.38, 0.375, 0.37, 0.34, 0.33, 0.32, 0.315, 0.31,
     8    0.30, 0.295, 0.29, 0.28, 0.27, 0.25, 0.24, 0.23, 0.22,
     9    0.21, 0.20, 0.195, 0.19, 0.18, 0.172, 0.17, 0.16, 0.15,
     a    0.145, 0.14, 0.13, 0.121, 0.115, 0.102, 0.095, 0.091,
     b    0.0808, 0.076, 0.073, 0.070, 0.067, 0.062, 0.056, 0.054,
     c    0.047, 0.041, 0.037, 0.034, 0.029, 0.025, 0.020, 0.015,
     d    0.006, 0.001/


        dx = 0.01d0
        eps = 1.d-4

!          integral of  dn/dx from 0 to 1. 
        call ktrpzIntegT(dndx, 1, n, 0.d0, dx, 1.d0, s)
!
!         make normalized table
!
        do i = 1, n
           ndndxn(i) = dndx(i)/s
        enddo
!         make table of intgral(0:x) of ndndx. for x = 0 to 1.0 step
!         0.01
        x = 0.
        do i=2, n 
           x = x + dx
           call ktrpzIntegT(ndndxn, 1, n, 0.d0, dx, x, intendndxn(i) )
        enddo
        intendndxn(1) = 0.
        intendndxn(n) = 1.
!           solve  inte(0:x) of ndndx = u for u = 0 to 1.0 step 0.01
        x = dx
        uconst = 0.
        x1 = 0.
        do i = 2, n-1
           uconst = uconst + dx
           call kbinChop(cinteLdndxn, x1, 1.d0, x, eps, ans, icon)
           if(icon .ne. 0) then
              call cerrorMsg(
     *        'failed in making leading particle sampling table', 0)
           endif
           ppsxn(i) = ans        
           x = ans
           x1 = x
        enddo
        ppsxn(1) = 0.
        ppsxn(n) = 1.
      end
      real*8 function cinteLdndxn(x)
      implicit none
      real*8 x

#include "Zcinippxc.h"

      real*8 ans

      call ktrpzIntegT(ndndxn, 1, n, 0.d0, dx, x, ans)
      cinteLdndxn = ans - uconst
      end


