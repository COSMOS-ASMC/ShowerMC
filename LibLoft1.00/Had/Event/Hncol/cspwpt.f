!c             test cspwpt
!       include '../KKlib/rnd.f'
!c------------------------------------
!       real*8 b, p, pt
!       integer i
!
!       b = 1.7d0
!       p = 3.2d0
!       do i = 1, 10000
!           call cspwpt(b, p, pt)
!           write(*, *) sngl(pt)
!       enddo
!       end
        subroutine cspwpt(b, p, pt)
!          sample pt from   x / (x+b)**p dx
!          solve  ( (1+x/b)(1-p) - (2-p) )(1+x/b)**(1-p) =-u
!          for b=3, p=15, <pt>=.5 gev.
!              b=3, p= 9, <pt>=1. gev.
!              b=6, p=15, <pt>= 1. gev
!     <pt>= 2*b/(p-3)
!              1872 msec / 50000  (b=3, p=9)
         real*8 b, p, pt
!
         real*8 eps0/1.d-3/
         real*8 u, x1, tmp1, tmp2, f, fp
         call rndc(u)
         x1=sqrt( 1. - u) /(p-1.)
!         *** until loop*** 
         do while (.true.)
             tmp1=(1.d0+x1)**(1.d0-p)
             tmp2= 1.d0+(p-1.d0)*x1
             f=tmp1* tmp2 -u
             fp=(p-1.d0)*tmp1  * (1.d0 - tmp2/(1.d0+x1))
             x=x1- f/fp
             if(x .lt. 1.d0) then
                eps=abs(x-x1)
             else
                eps=abs(( x-x1)/x )
             endif
             x1=x
            if(eps .lt. eps0) goto 100
         enddo
  100    continue
         pt=x1*b
         end
