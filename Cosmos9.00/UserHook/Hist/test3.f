      program main
      use modHistogram
      use modHistogram1
      use modHistogram3
      implicit none
      integer ns, nrbin, nebin, fno
      parameter (ns=3, nrbin=10, nebin=8)
      type(histogram3)  h(ns)
      save h
      real*8 x, y, z, fx, R, tm
      integer i, j
!       assume  Energy spectrum  E**-adE
!               R spectrum       R**-bdR
!               T spectrum       Texp(-T/c)dT
!     where
!        a=2log(R)-(s-1).  (  E>1 )
!        b=1-(s-1)          R<30   ( R>1  )
!          3+(s-1)          R>30
!        c= 0.1 + 10.0/E + R**2/300.+(s-1)
! 
! We take T-spectrum, regarding R and E as parameters.
! for s=0.8, 1.0, 1.2
! Let's see T for R= 1 to 500 with log10 step 0.5 (bin=6)
!                 E= 1 to  30 with log10 step 0.5 (bin=3)
!        
      real*8 coef(2),  pw(2), node(3), xp(1), ca(2)
      real*8 a, b, c, u, E, s(ns)
      character*24 dirstr, key
      data s/0.8, 1.0, 1.2/ 

      call kwhistso(1)
      fno=-6
      do i = 1, ns 
         call kwhisti3(h(i),
     *     1.0, 0.2, nrbin, b'00001',   !  R
     *     1.0, 0.2, nebin, b'01001',    !  E count  overflow 
     *     0.1, 0.1, 50,    b'00011' )   ! T
         call kwhistai3(h(i),
     *   "Time spectrum",
     *   "ret", "ptcls", .true., 0.,
     *   "r", "m", "E", "GeV", "T", "ns")
         call kwhistc3(h(i))
         pw(1)  = 1.d0+(s(i)-1.)
         pw(2)  = 3.d0+ (s(i)-1.)
         node(1) = 1.d0
         node(2) = 30.d0
         node(3) = 1.d5
         coef(1) = 1.d0
         coef(2) = 1./30.d0

         write(dirstr,'("s",i1,"/")') i
         call kwhistdir3(h(i), dirstr)
         write(key, '("s=",f5.2," r:")') s(i)
         call kwhistid3(h(i), key)
         call ksampPwX(coef, pw, node, 2, xp, ca)
         do j = 1, 300000
            call ksampPw(j,  ca, pw, node, 2,  R, fx)
            a = -2.0 - log(R) + (s(i)-1.0)
            call rndc(u)
            E = u**(1./(a+1.))
            c = 0.1 +  10./E + R**2/300. + (s(i)-1.)
            call ksgmis(1, c, tm)
            call  kwhist3(h(i), sngl(R), sngl(E), sngl(tm), 1.0)
         enddo
         call kwhists3( h(i), 0. )
         call kwhiststep3(h(i), 2, 2)
         call kwhistp3( h(i), fno)
      enddo
      end program main
      
