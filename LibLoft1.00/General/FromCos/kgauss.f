!         test kgauss  on sun fortran (new version)
!      
!       real*8 g1
!       do 100 i=1, 25000 
!          call kgauss(0.d0, 1.d0, g1)
!         write(*, *) sngl( g1 )
!  100  continue
!      end
      subroutine kgauss(m, v, g1)
!           kagauss and kgauss2:
!           Generate a Gaussian random variable with a given mean
!           and variance.
!        m:   real*8 input.   mean of the distribution.
!        v:   real*8 input.   variance of the distribution.
!       g1:   real*8 output.  a gaussian random variable.
!
!       g2:   real*8 output.  another gaussian random variable.
!             available only via kgauss2
!
!        subroutine needed:  rndc.
!
!          Since Box-Muler method with a switch to generate one
!          variable at a time efficently is not usable for
!          Cosmos which uses the skeleton-fleshing method,
!          this version generate two varibles at a time and
!          returns only one (kgauss), or two (kgauss2).
!          (This one is faster than Butcher's method).
!        100000 variable generation: by 90 MHz Pentium and Absoft Fortran.
!         kgauss needs  (1.29 sec) and 
!         kgauss2 needs (0.69 sec )
!         Butcher's method needs 1.69 sec.
!   
!          
!          To generate 25000 variables, 0.85 sec is needed on Sparc 2.
!           
        implicit none
        real*8 m, v, g1, g2
!
!        logical sw/.true./
        logical more
        real*8 u1, u2, r
        real*8  temp
!       save sw, temp
!        save u2
        logical entry2
!
        entry2 =.false.
        goto 10

!       ***************
        entry  kgauss2(m, v, g1, g2)
!       ***************
        entry2 = .true.

!       integer nt
!             -----------------
!        if(sw) then
!               counter for the loop at the test time.
!           nt=0
!               nt is distributed like
!               exp(-1.57nt) ( nt=1, 2, ...).  The average is 1.75.
 10      continue
            more=.true.
            do while (more)
!              nt=nt+1
               call rndc(u1)
               call rndc(u2)
               u1=u1*2 - 1.0d0
               u2=u2*2 - 1.0d0
               r=u1**2 + u2**2
               more= r .gt. 1.0d0
            enddo
!           write(*,*) nt
            temp=sqrt(-2*log(r)/r)
            g1=u1*temp*v + m
!           sw=.false.
            if(entry2) then
!        else    
!                g1=u2*temp*v + m
               g2 = u2*temp*v + m
            endif   
!            sw=.true.
!        endif    
      end    
