!       include 'rnd.f'
!--------------------------------
!c           to test sampling from the negative binomial dist.
!       implicit none
!       real*8 k, aven
!       integer  n, i
!       aven = 10.
!       k = 1.5 
!       do i=1,10000
!           call knbino(k, aven, n)
!           write(*, *) float(n), float(n)/sngl(aven)
!       enddo
!      end
       subroutine knbino(k, avn, n)
!         samples n from the negative binomial distrubtion of
!         p(n)= (n+k-1)] / n]/(k-1)] ( a/(1+a) )** n / (1+a)**k
!         where a=avn/k.   k should be real*8
!         this is for unfixed k and avn.  If either or both of them are
!         fixed, user anothr method for faster sampling.
!
       implicit none
       real*8 k, avn
       integer n
!
       real*8 a, pn, ar, u, sum
       integer  nmax
!
        a=avn/k
        pn= 1.d0/(1.d0+a)**k
        ar=a/(1.d0+a)
        call rndc(u)
        sum=pn
        nmax=avn*15.d0
        n=0
        do   while (sum .lt. u .and. n .lt. nmax)
            n=n+1
            pn= pn* (n+k-1)/n * ar
            sum=sum+pn
        enddo
       end
