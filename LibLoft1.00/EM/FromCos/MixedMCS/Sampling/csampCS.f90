subroutine cfixHardMuc(Ein, A, muc, uc)
!  Prepare for  single coulomb scattering anglesampling.
!  fix A,muc,uc. At sampling
!  only mu>muc is sampled. muc is fixed  by getMuc
  use modMCS0
  use modcMCS
  implicit none

  real(8),intent(in):: Ein !  eV. e KE energy
  real(8),intent(out):: A  !  reduced sampling parameter
  real(8),intent(out):: muc ! mu> muc is hard scatterign
  real(8),intent(out):: uc  ! corresponding uniform random
           ! number is uc. for sampling, u> uc must be used. 
  real(8):: logE, A0,  v, error

  logE = log(Ein)
  call kcsplIntp(logKEele, TPXSnow%A0, nEneg, TPXSnow%coefA0,&
       nEneg-1, logE, A0)
  A = A0

  call cgetMuc(Ein, muc)

  if( muc>=1.0d0) then  ! never happen
     uc = 1.d0
  else
!         To sample  only mu> muc, v> vc=muc(A+1)/(A+muc)
!         so we must get vc and corresponding uc 

     call cgetUc(MCSnow%sampleV,  Ein, A, muc, uc)

  endif
!  write(0,*) ' Ein=',Ein, ' A=',A
!  write(0,*) ' uc =',uc, ' muc=',muc
end subroutine cfixHardMuc


subroutine csampCSPolA(Ein, A, uc, mu, cosa ) 
  use modcMCS  
  implicit none
  real(8),intent(in):: Ein !  eV. e KE energy
  real(8),intent(in):: A   !  reduced pdf parameter
  real(8),intent(in):: uc  !  uniform random number lower limit
             ! uc=0--> all angle DCS sampling
  real(8),intent(out):: mu  ! sampled angle.= (1-cosa)/2
  real(8),intent(out):: cosa  ! 
  real(8)::  u,  v, error

  call rndc(u) 
  u = (1.0-uc)*u + uc

!          for u ~ 1.0-eps,  eps <~ 10^-4, mu becomes bit
!          larger than what should be.   this cannot be
!          mproved by putting 3 or 4  instead of 2 
!          after nEneg.

  call  kpolintp2(uarray, 1, 0.d0, KEele, 1, 0.d0, &
       MCSnow%sampleV,  usize, usize, nEneg, 2, 2, u, Ein, v,  error)
  mu =min( A*v/(A+1-v), 1.0d0)  ! avoid mu= 1.0000000322 etc 
  cosa =1- 2*mu
end subroutine csampCSPolA

subroutine cgetUc(sampleV, Ein, A, muc, uc)
  use modcTPXS
  use modMCS0
  implicit none
  real(8),intent(in):: sampleV(usize, nEneg)
  real(8),intent(in):: Ein ! e-/e+ energy in eV.
  real(8),intent(in):: A   !  A0 for Ein
  real(8),intent(in):: muc ! mu>= muc is to be sampled
  real(8),intent(out):: uc ! for that, u>= uc must be
                     ! used.
  integer::ie1, ie2
  real(8):: vc, uc1, uc2, w1, w2, error, temp
  if( muc > 0. ) then
     vc = muc*(A + 1) /(A+muc)
     call kdwhereis(Ein, nEneg, KEele, 1,  ie1)
     if( ie1 == nEneg) then
        ie2 = ie1 
        ie1 = ie2 -1
     else
        ie2 = ie1 + 1
     endif

!     write(0,*) ' muc=',muc, ' A=',A, ' ie1, ie2=', ie1, ie2
!     write(0,*) 'Ein=',Ein, ' vc=', vc

     call kpolintpLogxyFE(sampleV(1, ie1), 1, uarray, 1, usize, 3, &
       2,  vc, uc1, error)


     call kpolintpLogxyFE(sampleV(1, ie2), 1, uarray, 1, usize, 3, &
       2,  vc, uc2, error)
     temp = log(KEele(ie2)/KEele(ie1))
     w1 = log(KEele(ie2)/Ein)/temp
     w2 = log(Ein/KEele(ie1))/temp
     uc = w1*uc1 +w2*uc2
  else
     uc = 0.
  endif
end subroutine cgetUc

subroutine cu2vRedPDF(sampleV, Ein, u, v)
  use modcTPXS

  use modMCS0
  implicit none
  real(8),intent(in):: sampleV(usize, nEneg)
  real(8),intent(in):: Ein ! e-/e+ energy in eV.
!  real(8),intent(in):: A   !  A0 for Ein
  real(8),intent(in):: u   ! 
  real(8),intent(out):: v !

  integer::ie1, ie2
  real(8):: vc, v1, v2, w1, w2, error, temp
  call kdwhereis(Ein, nEneg, KEele, 1,  ie1)
  if( ie1 == nEneg) then
     ie2 = ie1 
     ie1 = ie2 -1
  else
     ie2 = ie1 + 1
  endif
  call kpolintpLogxyFE(uarray, 1, sampleV(1, ie1), 1, usize, 3, &
       1,  u,  v1, error)
  call kpolintpLogxyFE(uarray, 1, sampleV(1, ie2), 1, usize, 3, &
       1,  u,  v2, error)
  temp = log(KEele(ie2)/KEele(ie1))
  w1 = log(KEele(ie2)/Ein)/temp
  w2 = log(Ein/KEele(ie1))/temp
  v = w1*v1 +w2*v2
end subroutine cu2vRedPDF
