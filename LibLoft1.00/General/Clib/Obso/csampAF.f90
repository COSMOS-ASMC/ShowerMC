subroutine csampAF0(iowk, filen, sampInfo)
  use modSampAF
  implicit none
  integer,intent(in)::iowk  !  file logical number temporarily used
  character(len=*),intent(in)::filen  ! file name which contains (x,dn/dx)
  type(sampAF),intent(out)::sampInfo  ! distribution info is stored

  integer::icon
  integer::i
  call copenf(iowk, filen, icon)
  if(icon \= 0 ) then
     write(0,*) ' error '
     write(0,*) 'file: ',filen
     write(0,*) ' could not be opened'
     stop
  endif
  i = 0
  call cskipComment(iowk, icon)
  if(icon .ne. 0 ) stop
  do while ( .true. )
     read(iowk, *, end=100) sampInfo%x(i+1), sampInfo%y(i+1)
     i  = i + 1
  enddo
100 continue
  close(iowk)
  sampInfo%n =  i
  call ksampAF0(sampInfo%x, sampInfo%y, sampInfo%n, &
         sampInfo%coef, sampInfo%n, sampInfo%yi,  &
         sampInfo%sum,  sampInfo%coef2)
end subroutine csampAF0


subroutine csampAF(sampInfo, xs)
  use modSampAF
  implicit none
  type(sampAF),intent(in)::sampInfo
  real(8),intent(out):: xs
  call ksampAF(sampInfo%x, sampInfo%yi, sampInfo%n, &
               sampInfo%coef2, sampInfo%n, xs)
end subroutine csampAF
   
subroutine csampAFIntp(sampInfo, xv, ans)
!         this is not for sampling but simply
!       get value of y at xv
  use modSampAF
  implicit none
  type(sampAF),intent(in):: sampInfo ! input obtained by csampAF0
  real(8),intent(in)::  xv  ! input
  real(8),intent(out):: ans  ! y at xv
  call  kcsplIntp(sampInfo%x, sampInfo%y, sampInfo%n, &
        sampInfo%coef, sampInfo%n, xv, ans)
end subroutine csampAFIntp

subroutine csampAFmax(sampInfo,  xmax, fmax)
!       find max position and value of given function
!      (approx value)
  use modSampAF
  implicit none
  type(sampAF),intent(in):: sampInfo ! input obtained by csampAF0
  real(8),intent(out):: xmax !  max position  in (x1,x2) ; approx value
  real(8),intent(out):: fmax !  max function value

  real(8):: x1, x2
  real(8):: x, dx, temp
  integer:: i

  x1 = sampInfo%x(1)
  x2 = sampInfo%x(sampInfo%n)
  x = x1
  dx = (x2-x1)/sampInfo%n/10.
  xmax = x1
  call csampAFIntp(sampInfo, xmax, fmax)
  call csampAFIntp(sampInfo, x2, temp)
  if( fmax < temp ) then
     xmax = x2
     fmax = temp
  endif
  x = x + dx
  do while (x .lt. x2-dx/2 )
     call csampAFIntp(sampInfo, x, temp)
     if( fmax <  temp ) then
        xmax = x
        fmax = temp
     endif
     x = x + dx
  enddo
end subroutine csampAFmax
