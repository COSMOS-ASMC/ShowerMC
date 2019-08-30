!c          for test of kbetar generator
!          implicit none
!          real*8 a, b, x
!          integer i
!          call cerrorMsg('Enter a, b', 1)
!          read(*, *) a, b
!          do  i=1, 100000
!             call kbetar(a, b, x)
!             write(*, *) sngl(x)
!          enddo
!          end
!- -------------------------- for drawing curve ---------------
!
!          anrm=gamma(a)*gamma(b)/gamma(a+b)
!
!             f=      x**(a-1.)* (1.-x)**(b-1.)/anrm
!       **************************
        subroutine kbetar(a, b, x)
!       **************************
        implicit none
        real*8 a  ! input. see below
        real*8 b  ! input. see below
        real*8 x  ! output. sampled random variable

!           generate a random variable x following beta function kernel
!           i.e.  density function is x**(a-1)* (1-x)**(b-1)dx
!
!          both of a and b cannot be < 1.0
!              if a and b unbalance largely, rejection rate will
!              become large: (a=2, b=5--->1/6.6 is accepted
!                             a=2.1,b=1.2-->30 % rejection)
!        To sample  100000 variables, 4.226 sec is needed.
!       (wiht a=1.5, b=2.5 with facom m780)
!       this is sometimes faster than ksbeta for a class of a,b
!       but has limitations.

        real*8 ai, bm, u, bi, am

        if(a .lt. 1. and. b .lt. 1.) then
           write(*,*) ' kbetar cannot accept a<1,b<1:'
           write(*,*) ' use, ksbeta in such case'
           write(*,*) ' a=',a, ' b=',b
           stop
        endif
        if( (a .ge. b .and. b .ge. 1.)  .or.  a .lt. 1.) then
           ai=1./a

           bm=b-1.
!             *** until loop*** 
           do while (.true.)
              call rndc(u)

              x=u**ai
              call rndc(u)
              if         ( u .le. (1.-x)**bm)
     *             goto 100
           enddo
 100       continue
        else
           bi=1./b
           
           am = a-1.
!            *** until loop*** 
           do while (.true.)
              call rndc(u)
              x=1.-u**bi
              call rndc(u)
              if         ( u .le. x**am)
     *             goto 200
           enddo
 200       continue
        endif
        end
