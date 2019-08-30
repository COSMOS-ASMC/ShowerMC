!     **************************************************************
!     *
!     * csampPolAng: sample a polar angle from a distribution
!     *         (a+b*cos(t)+c*cos(t)**2 ) dcos(t)
!     *
!     **************************************************************
       subroutine csampPolAng(a, b, c, cost, icon)
       implicit none
       real*8 a, b, c, cost
       integer icon

       real*8 eps/1.e-3/, u, xn, tmpc, x, tmp
       integer nc
!
           call rndc(u)
           xn=2*u-1.
           tmpc= (2*a+2./3.*c)*u
           nc=0
!           *** until loop*** 
           do while (.true.)
               x=xn
               xn=x -
     *         ((x**3+1.)*c/3 + (x**2-1.)*b/2 +a*(x+1.)
     *          - tmpc )   / ( (c*x + b)*x + a )
                tmp=abs(x-xn)
               nc=nc+1
               if(nc .gt. 15 .or.  abs(xn) .gt.  1.) then
                   icon=1
                   goto 900
               endif
           if         (tmp .lt. eps)
     *                        goto 100
           enddo
  100      continue
           icon=0
           cost=x
  900      continue
        end
