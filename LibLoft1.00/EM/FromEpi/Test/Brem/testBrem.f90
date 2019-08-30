program main


  implicit none
  real(8):: Ek
  real(8),parameter::  me=0.511d-3

  
  real(8)::  teta,  Eg, Ee, Z 
  logical,save :: pr=.true.
  integer:: i, n, ab

  Z= 26.

  write(0,*) ' Enter n ab  Ek, Eg <Ek, in MeV  print'
  read(*,*) n, ab,  Ek, Eg, pr

  Ek = Ek*1.d-3
  Eg = Eg*1.d-3
  
  Ee =(Ek + me)

  if(pr) then
     do i = 1, n
        if(ab == 1) then
           call epBremAng(Ee, me, Eg, Z, teta)
        elseif( ab == 2) then
           call epSmpBremAng2BN(Ee, Eg, teta)
        else
           call epSmpBremAngTsaiFE(Ee, Eg, teta)
        endif
        write(*,*) teta
     enddo
  else
     do i = 1, n
        if(ab == 1) then
           call epBremAng(Ee, me, Eg, Z, teta)
        elseif( ab == 2) then
           call epSmpBremAng2BN(Ee, Eg, teta)
        else
           call epSmpBremAngTsaiFE(Ee, Eg, teta)
        endif
     enddo
  endif
  write(0,*) i, teta  
end program main
