!      test kdblexpi1
!      implicit none
!     real*8 a/0./, ans, eps/1.d-10/
!     external func
!     real*8 func
!     call kdblexpi1(func, a, eps, ans)
!     write(*,*) ' ans=',ans, asin(1.d0)
!     end
!     real*8 function func(x)
!     real*8 x
!     func = 1./(1.+x**2)
!     end
!***************************************************
!      integral of a given function from a to inf.
!   by Mori & Moriguchi's double exponentil forumula.
!   This is adpted version of Maruzen lib.
!
!   /usage/  call kdblexpi1(func, a, eps, ans)
!    func: real*8. input. integrand function name.
!       a: real*8.  inut. the lower bound of the
!                 integral region.
!     eps: real*8. input. absolute error tolerance
!    ans : real*8. output. obtained integral.
!
!*************************************************

      subroutine kdblexpi1(func, a, eps, ans)
      implicit none
!----      include 'kdblexpi1.h'
#include "KKlib/kdblexpi1.h"

      real*8 func, a, eps, ans, epsv
      external func
      real*8 h
      real*8 vold, vnew, wm, wp
      real*8 eps0/1.0d-32/, epsq, epsm/0./, epsp/0./

      logical first/.true./
      save first
      integer is, ih, km, kp, nm, np, i, mstep
!
      if( first ) then
         call kdblexpi1aux
         first = .false.
      endif
         
      neval = 0

      if (abs(eps) .ge. eps0) then
        epsv = abs(eps)
      else
        epsv = eps0
      end if

      epsq = 0.2 * sqrt(epsv)

      h = half

      is = 2**npow
      ih = is

      l = 1

  101 continue

      km = 0
      kp = 0
      nm = 0
      np = 0

      vnew = 0

!     ---- initial step ----
!          integrate with mesh size = 0.5
!          and check decay of integrand
!
       do   i = is, nend, ih
        if (kp .le. 1) then
          wp = func(ap(i,l) + a) * bp(i,l)
          neval = neval + 1
          vnew = vnew + wp
          if (abs(wp) .le. epsv) then
            kp = kp + 1
            if (kp .ge. 2) then
              np = i - ih
              go to 111
            end if
          else
            kp = 0
          end if
        end if
       enddo
 111   continue
      if (l .le. 2) then
        if (np .eq. 0) then
          l = l + 1
          go to 101
        end if
      end if

       do   i = is, nend, ih

        if (km .le. 1) then
          wm = func(am(i,l) + a) * bm(i,l)
          neval = neval + 1
          vnew = vnew + wm
          if (abs(wm) .le. epsv) then
            km = km + 1
            if (km .ge. 2) then
              nm = i - ih
              go to 121
            end if
          else
            km = 0
          end if
        end if
       enddo
  121 continue

      vnew = vnew + func(a0(l) + a) * b0(l)
      neval = neval + 1

      if (nm .eq. 0) then
        nm = nend
        epsm = 0.2 * sqrt(abs(wm))
        write (*,
     *  '("kdblexpi1: Warning. slow decay on negative side ")')
      end if

      if (np .eq. 0) then
        np = nend
        epsp = 0.2 * sqrt(abs(wp))
        write (*,
     * '("Kdblexpi1: Warning. slow decay on positive side")')
      end if
!
      epsq = max(epsq, epsm, epsp)

!     ---- general step ----

      vold = h * vnew

       do   mstep = 1, npow
        vnew = 0.0
        ih = is
        is = is / 2
         do   i = is, nm, ih
          vnew = vnew
     $         + func(am(i,l) + a) * bm(i,l)
          neval = neval + 1
         enddo
         do   i = is, np, ih
          vnew = vnew
     $         + func(ap(i,l) + a) * bp(i,l)
          neval = neval + 1
         enddo

        vnew = (vold + h * vnew) * half
        if (abs(vnew - vold) .lt. epsq) then
!        ---- converged and return ----
          ans = vnew
          return
        endif
        h = h * half
        vold = vnew
      enddo
      write (*,
     * '("kdblexpi1: Warning. Insufficient mesh refinement")')
      ans = vnew
      end
!      ************************
      subroutine kdblexpi1aux
      implicit none
!----      include 'kdblexpi1.h'
#include "KKlib/kdblexpi1.h"
!
!     generate points and weights for double
!     exponential forula
!     over half infinite interval (a,infinity)        

!         l for de half infinite transformation       
!         l = 1  x = exp(0.5*t - exp(-t))            
!         l = 2  x = exp(t - exp(-t))                
!         l = 3  x = exp(2 * sinh t)                 

 
      real*8 h,eh, en, eni
      real*8 sh,ch
      integer n6
      parameter (n6 = 6)
      integer i

      npow = n6
      nend = 5 * 2**(npow+1)
      h = one / 2**(npow+1)
      eh = exp(h)

!     ---- de transformation x = exp(0.5*t-exp(-t)) ----

      a0(1) = exp(-one)
      b0(1) = 1.5d0 * a0(1)
      en = 1.0

       do   i = 1, nend
        en = en * eh
        eni = 1 / en
        sh = half * h * i
        ap(i,1) = exp(sh - eni)
        bp(i,1) = (half + eni) * ap(i,1)
        am(i,1) = exp(-sh - en)
        bm(i,1) = (half + en) * am(i,1)
       enddo

!     ---- de transformation x = exp(t-exp(-t)) ----
!                         l = 2

      a0(2) = exp(-one)
      b0(2) = 2.0 * a0(2)
      en = 1.0

       do   i = 1, nend
        en = en * eh
        eni = 1 / en
        sh = h * i
        ap(i,2) = exp(sh - eni)
        bp(i,2) = (1 + eni) * ap(i,2)
        am(i,2) = exp(-sh - en)
        bm(i,2) = (1 + en) * am(i,2)
       enddo

!     ---- de transformation x = exp(2*sinh t) ----
!                        l = 3

      a0(3) = 1.0
      b0(3) = 2.0
      en = 1.0

       do   i = 1, nend
        en = eh * en
        eni = 1 / en
        sh = en - eni
        ch = en + eni
        ap(i,3) = exp(sh)
        bp(i,3) = ch * ap(i,3)
        am(i,3) = 1 / ap(i,3)
        bm(i,3) = ch * am(i,3)
       enddo
      end
