!   ifort  CompMolFunc.f90   -L$LIBLOFT/lib/PCLinuxIFC64 -lloft

!    compute Moliere's F0, F1, F2  using epMolfunc.f90
! which is
! 2pi  f(teta) teta dteta = ( F0(x)  + F1(x)/B   F2(x)/B^2) dx
!                 F0= 2exp(-x)
!                 F1 = epMoliere1(x)
!                 F2 = epScotD2(x) (don't use epMoliere2)
!   x = teta^2 (teta reduced angle; rad^2)
!  
!
 program main
 implicit none
 real(8):: teta, x
 real(8):: epMolfunc1, epScotD2
 real(8):: f0, f1, f2, B
 integer,parameter::nB=10
 real(8):: Ba(nB)
 integer:: i
Ba(:)=0.
write(0,*) "Enter list of B's with last / (<=10)"
read(*,*) Ba



do i = 1, nB
   B=Ba(i)
   if(B == 0.) stop
   teta = 0.
   do while( teta<14.)
      x = teta**2
      f0 = 2*exp(-x)
      f1 = epMolfunc1(x)
      f2 = epScotD2(x)
      write(*,'(1p,7g14.5)') B, teta, x, f0, f1/B, f2/B**2, f0+f1/B+f2/B/B
      teta = teta + 0.01d0
   enddo
enddo
end program main
  
