!  Rtoutines for sampling a random variable  related to Brems angle
!  but  the sampled variable could be > pi at this level.

! epSmpBARF_m1: sample x from  p* x/( (x/q)**4 + 1) dx ~   x/(x**4 + q**4)  dx
! epSmpBARF_m2: sample x from  r* x/( (x/s)**2 + 1)**2 dx ~ x/(x**2+s**2)**2 dx
! epSmpBARF_m: sample x form  m = m1 +  m2
! epBARF_m1:   compute m1
! epBARF_m2:   compute m2
! epBARF_m:    compute m
module modBremAngRelF
  !  This is used only when we want to know which of m1  or m2 is used
  !  for sampling a variable in epSamp_m. (and to use pi, hpi)
  implicit none

  real(8),parameter::hpi=asin(1.0d0), pi=hpi*2
  integer,save:: midx

end module modBremAngRelF

subroutine epSmpBARF_m1(q, x)
  use  modBremAngRelF
  implicit none
  ! m1(x)  =  p*x/((x/q)**4 + 1)
  !      ~  x/(x**4 + q**4)
  !
  real(8),intent(in):: q  ! >0 ( could be negative but same as positive)
  real(8),intent(out):: x  ! sampled variable.   0.
      
  real(8):: u   ! uniform random variable.

  call rndc(u)
  x = sqrt(tan(u*hpi)) * q
end subroutine epSmpBARF_m1

subroutine epSmpBARF_m2(s,x)
  use  modBremAngRelF
  implicit none
  !  m2(x) =r*x/((x/s)**2 + 1)**2
  ! ~ x/(x**2+s**2)**2
  real(8),intent(in):: s  ! (>0 ;same as q)
  real(8),intent(out):: x ! sampled variable
  

  real(8):: u
      
  call rndc(u)
  x = sqrt(1.0d0/u-1.0d0) * s
end subroutine epSmpBARF_m2

subroutine epSmpBARF_m(p, q, r, s, x)
  use  modBremAngRelF
  implicit none
  real(8),intent(in):: p, q, r, s  ! input parametes.
     ! m(x) =   m1(p, q,x)  +  m2(r,s,x)
  real(8),intent(out):: x ! sampled varible
  ! Depending of q, s's parameter meaning.  this may be a reduced angle
  ! so that teta = x/g  (where g = Ee/Me;electorn gamma factor) is
  ! a real angle in teta (if < pi) .
  ! 


  real(8):: u,  s1, s2, ratio


  s1 = p*q**2* pi/4.0d0  ! area of m1
  s2 = r*s**2/2          ! area of m2 
  ratio = s1/(s1+s2)

  call rndc(u)
  if( u <  ratio ) then
     call epSmpBARF_m1(q,x)
     midx=1
  else
     call epSmpBARF_m2(s,x)
     midx=2
  endif
end subroutine epSmpBARF_m

subroutine epqSmpBARF(idx)
  use  modBremAngRelF
  !     query for  function index, idx ( i or 2 for m1 or m2)
  implicit none
  integer,intent(out):: idx
  idx = midx
end subroutine epqSmpBARF

function epBARF_m1(p,q,x)  result(ans)
  implicit none
  real(8),intent(in)::p, q, x
  real(8):: ans

  ans = p*x/((x/q)**4 + 1)
end function epBARF_m1

function epBARF_m2(r,s,x)  result(ans)
  implicit none
  real(8),intent(in)::r, s, x
  real(8):: ans

  ans = r*x/((x/s)**2 + 1)**2
end function epBARF_m2

function epBARF_m(p,q,r,s,x)  result(ans)
  implicit none
  real(8),intent(in)::p, q,  r, s, x
  real(8):: ans

  real(8),external:: epBARF_m1, epBARF_m2

  ans  = epBARF_m1(p,q,x) +  epBARF_m2(r,s,x)
end function epBARF_m





