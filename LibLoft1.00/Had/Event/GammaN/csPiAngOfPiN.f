!      implicit none
!      real*8 ecm, cs
!      integer cn, cpi, i
!      read(*,*)  ecm
!      cn=1
!      cpi=0
!      do i=1, 10000
!         call csPiAngOfPiN(ecm, cn, cpi, cs)
!         write(*, *) cs
!      enddo
!      end
!        decay angle of a pion from pi-n system is sampled
!    call csPiAngOfPiN(cmse, cn, cpi, cs)
!  ------------------------------
!  cmsein: mass of pi-n system in GeV
!      cn: charge of n after decay
!     cpi: charge of pi after decay
!      cs: sampled pion polare angle (cos)
!
!
      subroutine csPiAngOfPiN(cmsein, cn, cpi, cs)
      implicit none
      
      real*8  cmsein, cs
      integer cn, cpi
      
!
      real*8  cmse, a, b, c
      integer icon
!      
!        a+,b+,c+:              for (gamma + p)---> (n,pi+)
!        a0,c0 (b0=0 or b+)     for (gamma + p)---> (p,pi0)
!    for (n,pi0) use  0
!    for (p,pi-) use  +
!               a+(cms gev)   + for (gamma + p)---> (n,pi+)
!       xmin=    1.239     xmax=    1.777
      real*8  ap, bp, cp, a0, c0
      real*8 x
!
       ap(x)=(( 523.8081    *x-2523.137    )*x+ 3994.917    )*x-2066.6
     *08
!               b+(cms)
!       xmin=    1.239     xmax=    1.777
       bp(x)=((( 551.7981    *x-3704.184    )*x+ 9214.102    )*x-10054.2
     *2    )*x+ 4057.190
!               c+(cms)
!       xmin=    1.239     xmax=    1.777
       cp(x)=((-337.0378    *x+ 1605.070    )*x-2511.868    )*x+ 1287.63
     *  5
!                a0(cms)
!       xmin=    1.239     xmax=    1.777
         a0(x)=((( 5566.031    *x-33254.71    )*x+ 73843.62    )*x-72204
     *.50    )*x+ 26239.76
!                c0(cms)
!       xmin=    1.239     xmax=    1.777
         c0(x)=(((-5099.891    *x+ 30526.40    )*x-67929.69    )*x+ 6658
     *6.12    )*x-24266.87
!
         cmse=min(cmsein, 1.8d0)
         cmse=max(cmse, 1.239d0)
         if(cn .eq. 0 .and. cpi .eq.  1) then
             a=ap(cmse)
             b=bp(cmse)
             c=cp(cmse)
         elseif(cn .eq. 1 .and. cpi .eq. 0) then
             a=a0(cmse)
             b=0.
             c=c0(cmse)
         elseif(cn .eq. 1 .and. cpi .eq. -1) then
             a=ap(cmse)
             b=bp(cmse)
             c=cp(cmse)
         else
             a=a0(cmse)
             b=0.
             c=c0(cmse)
         endif
!         *** until loop*** 
         do while (.true.)
            call csampPolAng(a, b, c, cs, icon)
            if(icon .eq. 0) goto 30
         enddo
   30    continue
      end
