
! generic gamma function
!  for Mac Gfortran, this maxnfac cannot be put in gammafunc
!  so  separate module is used.
module modgamma
  integer,parameter:: maxnfac=170  ! max N for N!
end module modgamma

module gammafunc
  interface kgamma
!         real(4) in/out  For real(8) in/out use
!                             kgamma(x) without "use gammafunc"
     function krgamma(x) result(ans)
       real(4),intent(in)::x
       real(4)::ans
     end function krgamma
!    complex(8) in/out 
     function kcgamma(z) result(ans)
       complex(8),intent(in)::z
       complex(8)::ans
     end function kcgamma
!         integer in / real(8) out
     function kigamma(n) result(ans)
       integer,intent(in):: n
       real(8):: ans
     end function kigamma
  end interface kgamma
end module gammafunc  

function kcgamma(z) result(ans)
  use modgamma
  implicit none
! Gamma func. by the  Lanczos method
! ref:    http://www.numericana.com/answer/info/godfrey.htm
!         relative accuracy <  10^-13.
  complex(8),intent(in):: z
  complex(8):: ans

  real(8),parameter::pi=3.14159265358979323846d0
  real(8),parameter::twopi=pi*2
  real(8),parameter::sqrt2pi=2.50662827463100050241d0
  complex(8):: zz
  real(8),parameter::c(11) =(/1.000000000000000174663d0, &
       5716.400188274341379136d0, &
       -14815.30426768413909044d0, &
       14291.49277657478554025d0, &
       -6348.160217641458813289d0, &
       1301.608286058321874105d0, &
       -108.1767053514369634679d0, &
       2.605696505611755827729d0, &
       -0.7423452510201416151527d-2, &
       0.5384136432509564062961d-7, &
       -0.4023533141268236372067d-8/)
  integer,parameter::g=9

  complex(8)::t, s, ss
  integer::k, m
  real(8),external:: kfactorial

  
  zz = z
  if( real(zz) >= 1.d0  .and. aimag(zz) == 0. ) then
     m = int(real(zz)) 
     if( m == real(zz) .and. m <= maxnfac ) then
        ans = cmplx( kfactorial(m-1), 0.d0, kind=8)
        return !   !!!!!!
     endif
  endif
  
  if(real(zz) < 0.0 ) then
     zz = -z
  endif
  t = zz + g

  s=0.0
  do k = g+2,2,-1
     s=s+c(k)/t
     t=t-1
  enddo

  s=s+c(1)
  ss=(zz+g-0.5)
  s=log(s*sqrt2pi) + (zz-0.5)*log(ss)-ss

  ans = exp(s)

  if( real(z) < 0. ) then
     ans = -pi/(zz*ans*sin(pi*zz))
  endif
end function kcgamma


function kfactorial(n) result(ans)
  use modgamma
  implicit none
  integer,intent(in):: n  ! one of 0,1,2...maxn
  real(8)::ans
  real(8)::fac(0:maxnfac)=0.
    
  integer:: i
  if( fac(0) == 0.d0) then
     fac(0) = 1.
     do i = 1, maxnfac
        fac(i) = fac(i-1)*i
        !          write(0,*) i, fac(i)
     enddo
  endif
  if(n >= 0) then
     if( n<= maxnfac) then
        ans = fac(n)
     else
        ans = fac(maxnfac)   
     endif
  else
     write(0,*) ' input n for kfactorial is ',n, ' invalid'
     stop
  endif
end function kfactorial

!     *************************************************************
!     *                                                           *
!     * krgamma: gamma function in real domain. 
!     *         single precision accuracy.
!     *   this is same as kgamma in KKlib but real(4) is used
!     *************************************************************
!
!  Usage:  y=krgamma(x).   x
!
!     Computes gamma(x) with 6 significant digit.  gamma(x)=factorial of
!     (x-1).  gamma(1)=gamma(2)=1 .........................................
!
!
!
!
function krgamma(x) result(ans)
  implicit none
  real(4),intent(in)::x
  real(4):: ans
!
  real*8 z, f, t, temp
!
  real(4),parameter::pi=3.141592653
!
  if( abs(x)  > 15.0) then
     z=x
     if( z <= 0 ) then
        f=pi/sin(pi*z)
        z=1.d0 - z
     endif
     temp = 2.506628274d0*exp(-z)*z**(z-0.5d0)*  &
          ((3.47222222222d-3/z+8.3333333133d-2)/z+1.d0)
     if(x <  0.0) then
        temp=f/temp
     endif
  else
     f=1.
     z=x
     do while (z >  3.0)
        z=z-1.
        f=f*z
     enddo
     do while (z < 2.0) 
        f=f*z
        z=z+1.
     enddo
     z=z-2.0
     t = (((((1.08298598d-2*z - 3.42705226d-3)*z + 7.7549276d-2)*z &
       + 8.01782477d-2)*z + 4.12102903d-1)*z +4.22766368d-1)* z +  &
       1.0000002d0
     if(x < 2.0) f=1.0d0/f
     temp=t*f
  endif
  ans = temp
end function krgamma

function kigamma(n) result(ans)
  implicit none
  integer,intent(in):: n
  real(8):: ans
  real(8),external:: kfactorial

  ans = kfactorial(n-1)  ! error check is in kfactorial
end function kigamma
