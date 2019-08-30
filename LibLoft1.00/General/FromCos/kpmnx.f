!          test kpmnx
!
!      include  "kgamma.f"
!      real*4 x, y, sint
!      integer m, n
!      real*4 kpmnxn
!c      
!      read(*,*) m, n
!      
!      do x=-1.d0, 1.0000001d0, 0.01d0
!          y = min(x, 1.d0)
!          sint =sqrt(1.d0- y**2)
!          write(*,*) m, n, y, kpmnxn(m, n, sint, y)
!      enddo
!      end
!
!        kpmnx: P(n, m, x) = m-th  derivative
!               of Legendre Pn(x) times
!              (1-x**2)**(m/2)
!
!          Use recurrence relation:
!      (n-m) P(n, m, x) = (2n-1)xP(n-1,m, x) - (n-1+m)P(n-2, m, x)
!    
!       P(0, 0, x) = 1, P(n, m, x) = 0 (m> n)
!       P(1, 0, x) = 2x
!       P(1, 1, x) = 2(1-x**2)**(1/2)
!       P(1, m, x) = 0 ( m>=2 ) 
!
      real*4 function krmnx(m, n, x)
      implicit none
      integer m  !  input  >=0
      integer n  !  input. >=0
      real*4  x  !  input. argument |x| <= 1
      

!         m-th derivative of  Pn(x)
!     where Pn(x) = 1/(2^2 n!) (x^2-1)^n
!
      real*4  r0m, r1m, rn, rn1, rn2, krnnx
      real*4 kpnx

      integer i

      if(n .lt. 0 .or. m  .lt. 0) then
         write(0, *) ' error input to krmnx; n, m=', n, m
         stop 9999
      endif
      if(m .eq. 0) then
!             m=0
         krmnx = kpnx(n, x)
      elseif(m .gt. n) then
!         krmnx = 0.d0
         krmnx = 0.
      elseif(n .ge. 1) then
         if(m .eq. 1) then
!            r1m = 1.d0
            r1m = 1.0
            r0m = 0.
         else
!              m > 1
            r1m = 0.
            r0m = 0.
         endif
         
         if( n .ge. 2) then
            rn2 = r0m
            rn1 = r1m
            do i = 2, n
               if(i .eq. m) then
                  rn = krnnx(i, x)
               else
!                  (n-m) r(n, m, x) = (2n-1)xr(n-1,m, x) - (n-1+m)r(n-2, m, x)
                  rn = ( (2*i-1)*x*rn1 - (i-1+m)*rn2 ) /(i-m)
               endif
               rn2  = rn1
               rn1 = rn
            enddo
            krmnx = rn
         else
            krmnx = r1m
         endif
      else
!          n=0, m=0
         krmnx = 1
      endif
      end
!     **********************************************
      real*4 function krnnx(n,  x)
      implicit none
!          kpmnx for n=m case. 
      integer n  ! input
      real*4  x  ! input.
!
!        r(n,m) = r(n-2,m) + (2n-1)*r(n-1, m-1)
!        since  n=m, r(n-2,m) is always 0.
!    So r(n) = (2n-1)*r(n-1)
!        r(0) = 1. 
!
      real*4 r
      integer i


      if(n .lt. 0) then
         write(0,*) 'error input to krnnx; n=',n
         stop 999
      else
!         r = 1.d0
         r = 1.0
         do i = 1, n
            r = (2*i-1)*r
         enddo
         krnnx = r
      endif
      end
!-------------------------------------------------------
!c            testing kpnx, kpmnx, kdpmnx, kdpnx, kdpmnxn, kpmnxn
!c         kpnx:  Legendre polynomial
!c         kpmnx:  Legendre polynomial p(m, n, x)
!c        kdpmnx:  d p(m, n, x)/dx 
!c         kdpnx:   m-th derivative of legendre function, kpnx
!c             =       p(m,n,x)/(1-x**2)**(m/2)
!c       kdpmnxn:  sqrt( em* (n-m)!/ (n+m)! ) * kdpmnx
!c        kpmnxn:  sqrt( em* (n-m)!/(n+m)!  ) * kpmnx
!      implicit none
!      integer m, n
!      real*4 x, y, y1, y2, z, kpmnx
!      real*8 kgamma
!      real*4 kpnx, kdpmnx, kdpmnxn
!      integer mv
!      read(*, *) mv
!      do  m=mv, mv
!c
!         do  n=1, 8
!            do  x=-.996, 1., .004
!c               y = kpmnx(m, n, x)
!c               y1 = kgamma(dble(n-m)+1.d0)
!c               y2 = kgamma(dble(n+m)+1.d0)
!c               z =sqrt((2*n+1)/2.d0*y1/y2)*y
!c                z = kpnx(n, x)
!c                 z = kdpmnx(m, n, x)
!                z = kdpmnxn(m, n, x)
!                write(*, *) sngl(x), sngl(z)
!            enddo
!            write(*,*)
!         enddo
!      enddo
!      end
!*********************************************************************
      real*4 function kpnx(n, x)
      implicit none
      real*4 x
      integer n
!           Legendre polinomial  (n>=0)
      real*4 pim, pimm, pi
      integer nc/0/, i
!
           if(n .eq. 0) then
!              kpnx=1.d0
              kpnx=1.0
           elseif(n .eq. 1) then
              kpnx=x
           elseif(n .ge. 2) then
              pim=x
!              pimm=1.d0
              pimm=1.0
              do   i=2, n
                  pi=( (2*i-1)*x *pim - (i-1)*pimm )/i
                  pimm=pim
                  pim=pi
              enddo
              kpnx=pi
           else
              if(nc .lt. 20) then
                 write(*,*) ' n=',n,' invalid for kpnx'
                 nc=nc+1
              endif
!              kpnx=1.d50
              kpnx=1.e35
           endif
       end
! **************************************
       real*4 function kpmnx(m, n, sint, x)
       implicit none
       real*4 x, sint
       integer m, n
!
       real*4 krmnx, sintm
       if(m .eq. 0) then
          sintm = 1.
       else
          sintm = sint**m
       endif
       kpmnx = sintm * krmnx(m, n, x)
       end
!      ********************************
       real*4 function kpmnxn(m, n, sint, x)
       implicit none
       real*4 x, sint
       integer m, n
!               sqrt(em*(n-m)!/(n+m)!)* kpmnx(m, n, x )
       real*4 kpnorm, kpmnx
       kpmnxn=kpnorm(m, n)* kpmnx(m, n, sint, x)
       end
!      *********************
       real*4 function kpmnxisin(m, n, sint, x)
       implicit none
       integer m, n   ! m >= 1 
       real*4 x      
       real*4 sint   !  sqrt(1.-x**2)
!     
       real*4 krmnx, sintm
       
       if(m .eq. 1) then
          sintm = 1.
       else
          sintm = sint**(m-1)
       endif
       kpmnxisin = sintm * krmnx(m, n, x)
       end
!       ********************
       real*4 function kpmnxisinn(m, n, sint, x)
       implicit none
       integer m, n  !  m >=1
       real*4 x
       real*4 sint  !  sint = sqrt(1-x**2)
       real*4 kpmnxisin, kpnorm
       kpmnxisinn = kpnorm(m, n)*kpmnxisin(m, n, sint, x)
       end
!        **********************
       real*4 function kdpmnxsin(m, n, sint, x)
       implicit none
       real*4 x
       real*4 sint
       integer m, n

       real*4 krmnx, sintm
       if(m .gt. 0) then
          if(m .eq. 1) then
             sintm = 1.
          else
             sintm = sint**(m-1)
          endif
          kdpmnxsin= (-m*x*krmnx(m, n, x) + sint*sint*krmnx(m+1, n, x))
     *              *sintm
       else
          kdpmnxsin=  sint*krmnx(1, n, x)
       endif
       end
       real*4 function kdpmnxsinn(m, n, sint, x)
       implicit none
       real*4 x
       real*4 sint
       integer m, n
!          
       real*4 kpnorm, kdpmnxsin
       kdpmnxsinn = kpnorm(m, n)* kdpmnxsin(m,  n, sint, x)
       end
! *************
       real*4 function kpnorm(m, n)
       implicit none
       integer m, n
       real*4 pnorms, em
       real*8 kgamma
       save pnorms
!
       integer ms/-1/, ns/-1/
       logical first


       integer nn, i
       parameter (nn=30)
       real*4 fact(0:nn), dv, dn
       save first, fact, ms, ns
       data first /.true./

       if(first) then
          fact(0) = 1.
          do i= 1, nn
             fact(i) = i* fact(i-1)
          enddo
          first = .false.
       endif
       if(m .eq. ms .and. n .eq. ns)then
       elseif(m .gt. n) then
          pnorms=0.0
       else
          ms = m
          ns = n 
          if(m .eq. 0)then
             em=1.0
          else
             em=2.0
          endif
          if(n-m  .le. nn) then
             dn = fact(n-m)
          else
             dn = kgamma(dble(n-m+1))
          endif
          if(n+m .le.  nn) then
             dv = fact(n+m)
          else
             dv = kgamma(dble(n+m+1))
          endif
          pnorms=sqrt(em* dn /dv)
       endif
       kpnorm=pnorms
       end



