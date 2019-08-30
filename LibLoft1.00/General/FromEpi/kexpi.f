!  test
!      implicit none
!      real*8 kexpi, x
!      integer i
!      x = 13.d0
!      do i = 1, 150
!         write(*, *) x, kexpi(x)
!         x = x + .01d0
!      enddo
!      end
!     ****************************************************************
!     *                                                              *
!     *    kexpi:  exponential integral                              *
!     *   keimlg:  kexpi(x)- log(x)                                  *
!     *   keiexp:  keimlg(x) * exp(-x)                               *
!     *                                                              *
!     ****************************************************************
!
!  /usage/  these are all functions.
!
!     definition:  kexpi(x)=  g  +  sum ( x**n/n/n' ) + log(x)
!
!                  where  g is euler const(0.577..)
!
!                  x may be <= 0, in keimlg or keiexp
!
!
!
      real*8 function kexpi(x)
      implicit none
      real*8 x
!
      real*8 kexpicore, s, g
      data g/0.577215664d0/

      s = kexpicore(x)
      if(x .lt. 14.) then
         kexpi=s + g + log(x)
      else
         kexpi=s*exp(x)
      endif
      end
!
!     ***********
      real*8 function keimlg(x)
!     ***********
!
!        ei(x)-log(x)
!
      implicit none
      real*8 x
      real*8 g, kexpicore, s
      data g/0.577215664d0/

      s = kexpicore(x) 
      if(x .lt. 14.) then
         keimlg = s + g
      else
         keimlg = exp(x)*s - log(x)         
      endif
      end

!     ***********
      real*8 function keiexp(x)
!     ***********
!
!     (ei(x)-log(x))*exp(-x)
!
      implicit none
      real*8 x
      real*8 g, kexpicore, s
      data g/0.577215664d0/

      s = kexpicore(x)
      if(x .lt. 14.) then
         keiexp = (s + g)*exp(-x)      
      else
         keiexp = s-exp(-x)* log(x)
      endif
      end
!     **********************************
      real*8 function kexpicore(x)
      implicit none
!
      real*8 x
      real*8 eps, s, tmp, fn, tmp2
!
      data eps/1.d-6/

      s=0.
      tmp=1.
      if(x .eq.  0. ) then
         kexpicore = 0.
      elseif(x .lt.14.) then
         fn=1.
 5       continue
         tmp=tmp*x/fn
         tmp2=tmp/fn
         s=s+tmp2
         if(abs(tmp2/s) .gt. eps) then
            fn=fn+1.
            go to 5
         endif
         kexpicore = s
      else
         fn=15.
 105     continue
         s=s+tmp
         fn=fn-1.
         if(fn .gt. 0.) then
            tmp=tmp*x/fn
            go to 105
         endif
         s=s/tmp/x
         kexpicore = s
      endif
      end         
