! This is a copy of kexpi.f from Epics.  The name is changed by
! attaching C in the last part of each name 
! test implicit none real*8
! kexpiC, x integer i x = 13.d0 do i = 1, 150 write(*, *) x, kexpiC(x) x =
! x + .01d0 enddo end
! **************************************************************** * * *
! kexpiC: exponential integral * * keimlg: kexpiC(x)- log(x) * * keiexp:
! keimlg(x) * exp(-x) * * *
! ****************************************************************
!
!  /usage/  these are all functions.
!
!     definition:  kexpiC(x)=  g  +  sum ( x**n/n/n' ) + log(x)
!
!                  where  g is euler const(0.577..)
!
!                  x may be <= 0, in keimlg or keiexp
!
!
!
      real*8 function kexpiC(x)
      implicit none
      real*8 x
!
      real*8 kexpicoreC, s, g
      data g/0.577215664d0/

      s = kexpicoreC(x)
      if(x .lt. 14.) then
         kexpiC=s + g + log(x)
      else
         kexpiC=s*exp(x)
      endif
      end
!
!     ***********
      real*8 function keimlgC(x)
!     ***********
!
!        ei(x)-log(x)
!
      implicit none
      real*8 x
      real*8 g, kexpicoreC, s
      data g/0.577215664d0/

      s = kexpicoreC(x) 
      if(x .lt. 14.) then
         keimlgC = s + g
      else
         keimlgC = exp(x)*s - log(x)         
      endif
      end

!     ***********
      real*8 function keiexpC(x)
!     ***********
!
!     (ei(x)-log(x))*exp(-x)
!
      implicit none
      real*8 x
      real*8 g, kexpicoreC, s
      data g/0.577215664d0/

      s = kexpicoreC(x)
      if(x .lt. 14.) then
         keiexpC = (s + g)*exp(-x)      
      else
         keiexpC = s-exp(-x)* log(x)
      endif
      end
!     **********************************
      real*8 function kexpicoreC(x)
      implicit none
!
      real*8 x
      real*8 eps, s, tmp, fn, tmp2
!
      data eps/1.d-6/

      s=0.
      tmp=1.
      if(x .eq.  0. ) then
         kexpicoreC = 0.
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
         kexpicoreC = s
      else
         fn=15.
 105     continue
         s=s+tmp
         fn=fn-1.
         if(fn.gt.0.) then
            tmp=tmp*x/fn
            go to 105
         endif
         s=s/tmp/x
         kexpicoreC = s
      endif
      end         
