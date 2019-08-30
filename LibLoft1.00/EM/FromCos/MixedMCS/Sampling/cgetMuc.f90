subroutine cgetMuc( Ein, muc)
!  get muc above which scattering is regarded as hard
!  This is diff. from the one in mkSampTbl; (Ein: need not
!  be on grid). If muc is 0, all angle is target and
!  detailed simulaton is intended
  use modcMCS
  implicit none
  real(8),intent(in):: Ein ! KE energy of e. in eV

  real(8),intent(out):: muc ! predefined critical angle .
         ! mu> muc is hard scattering region.

  real(8):: logE, y, error
  integer::ie

  
  logE = log(Ein)
  if(Ein >= MCSnow%minNon0mucE) then
     ie = MCSnow%minNon0mucEindex
     call kpolintpFE(logKEele(ie), 1, MCSnow%logmuc(ie), 1,  &
       nEneg-ie+1, 3,  logE, y, error)
     muc = exp(y)
  else
     muc = 0.
  endif
end subroutine cgetMuc
