! this is diff. from the one in mkSampTbl (cgetLambdah)
!  Ein may be not on the grid.
subroutine cgetLamh(Ein, lh)
!  get Lambda(h) which is max(lambda(el), c1*lambda1) 
  use modcMCS
  implicit none
  real(8),intent(in):: Ein ! electron energy in eV
                        ! 100<Ein<=1e9
  real(8),intent(out):: lh !  g/cm2

  real(8):: logE, error, y
  
  if(Ein>= 100. .and. Ein<=1.e9 ) then
     logE = log(Ein)
  else
     write(0,*) 'Ein =',Ein, ' in valid for cgetLambdah'
     stop
  endif

  call kpolintpFE(logKEele,1, MCSnow%loglambdah, 1, nEneg,3, &
       logE, y, error)
  lh = exp(y)
end subroutine cgetLamh

subroutine cgetTPMFP(Ein, ls1, ls2)
  use modcMCS
  implicit none
  real(8),intent(in):: Ein ! electron energy in eV
                        ! 100<Ein<=1e9
  real(8),intent(out):: ls1, ls2 !  g/cm2
         ! soft scattering transport mfp 1 and 2.
         ! soft is integration from mu=0 to muc
  real(8):: logE, error, y
  integer:: ie
  
  if(Ein>= MCSnow%minNon0mucE  .and. Ein<=1.e9 ) then
     logE = log(Ein)
  else
     write(0,*) 'Ein =',Ein, ' in valid forcgetTPMFP'
     stop
  endif

  ie = MCSnow%minNon0mucEindex
  call kpolintpFE(logKEele(ie), 1, MCSnow%loglambdas1(ie),1,  &
          nEneg-ie+1, 3,  logE, y, error)
  ls1 = exp(y)
  call kpolintpFE(logKEele(ie), 1, MCSnow%loglambdas2(ie),1, &
       nEneg-ie+1, 3,  logE, y, error)
  ls2 = exp(y)
end subroutine cgetTPMFP

subroutine cgetAveMu12(Ein, s, avemu, avemu2)
  use modcMCS
  implicit none
  real(8),intent(in):: Ein ! electron energy in eV
                        ! 100<Ein<=1e9
  real(8),intent(in):: s ! distacne between two 
                    ! hard elestic scatt. g/cm2

  real(8),intent(out):: avemu, avemu2
         ! soft scattering  average angle avemu = <mu>
         !  avemu2 = <mu2>
  real(8):: ls1, ls2

  call cgetTPMFP(Ein, ls1, ls2)

  avemu = (1.0d0 - exp(-s/ls1))/2
  avemu2 = avemu - (1- exp(-s/ls2) )/6.d0
end subroutine cgetAveMu12
