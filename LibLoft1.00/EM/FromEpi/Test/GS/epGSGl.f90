module GSmscat
  implicit none
  real(8):: GScoef1 ! sum 2pik^2N t 
  real(8):: GScoef2 ! sum 2pik^2N t * loggzai
  real(8):: loggzai ! see below
  real(8):: GSk
!  N # /cm3.  t length in cm
!   k: Eq.17  = Ze^2/mc^2 /gbeta^2
!             = Z e^2/(hbar c)( hbar c/mc^2)/gbeta^2
!             = Z/137 (hbar c/mec^2)(me/m)^2/gbeta^2
! for   k^2--> Z^2 should be Z(Z+1) for electron
!              Z^2 for heavy
!            hbar c/mec^2 = 197.3 x 10^-13/0.511   cm
!  loggzai = log( 1.10a/lambda bar)
!            a= a0/Z^(1/3). 
!         however,  Bethe: a=0.885a0/Z^(1/3)
!         a0; Bohr redius = hbarc/mec^2  hbarc/e^2
!                          
!         lambda bar/a0 = me/m /137 /gbeta
!         N= N0/A rho ni; ( rho t = g/cm2 );  ni = Ni/sum(Ni)
!         A= sum niAi; 

subroutine epGSGl(l,  Gl)
!  compute Gl; Eq.37 of Goudsmit and Saunderson
!         Phys.Rev. Vol.57.(1940)
!   Gl is the average of angle
  implicit  none
  integer,intent(in):: l  ! l in Eq. 37 of GS
  

  real(8),intent(out):: Gl   
  integer::i
  real(8)::sum
  sum = 0.d0
  do i = l,2,-1   ! if l=1, sum = 0.
     sum = sum + 1.d0/i 
  enddo
     
  Gl =exp(( -GScoef2 + GScoef1*sum)*l*(l+1) )
end subroutine epGSGl


