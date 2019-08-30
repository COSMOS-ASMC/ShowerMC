        subroutine cMWconstForMedia(md, pm, aTPXS)
        use modTPXS
        implicit none
#include "Zmedia.h"    
        type(epmedia),intent(in):: md
        integer,intent(in):: pm   ! >0 for e+, < 0 for e-
        type(TPXSconst),intent(inout):: aTPXS  !TPXSconst for media.

        type(TPXSconst),allocatable::zTPXS(:) ! for each Z in media
                                ! allocated here and deallocated 

        integer:: n  ! # of elements in a medium

        integer,parameter::io=11  ! temporaray disk logcial dev. #
                          ! should be the same as TempDeV
        integer:: i, j,  Z, icon, nE
        real(8)::MWavemu, MWavemu2, avemu, avemu2, temp, EK, S0, S1, S2

        n = md%noOfElem
        allocate(zTPXS(n))

        do i = 1, n
           Z = md%elem(i)%Z
                 ! create e-/e+ TPXS file
           call cTPXS_Z2file(Z, pm, filename)
           call copenf(io, filename, icon) 
           if(icon /= 0 ) then
              write(0,*) ' file=',trim(filename), ' cannot be opened '
              stop
           endif
           call cReadTPXS(io, pm, zTPXS(i))  ! for each Z set S0,S1,S2
        enddo
        nE = nEneg
        if(pm > 0 ) nE = nEpos
!      weigth is normalzied so that the absolute value of S0,S1,S2
!    for Air will be smaller than correct one but we deon't use
!    the abosulte values  (we use ratio or N*lamda*S=1 where S here is
!    S(true)/c, and N=N(true)*c so  
!    N*c lambdaa S/c  =1   gives correct lambda
!    write(0,*) ' weight'
!    write(0,*) md%elem(1:n)%NO    ! 
!!!
        aTPXS%n = nE
        do i = 1, nE
           temp = 0.
           do j = 1, n
              temp = temp + md%No(j) * zTPXS(j)%S0(i)
           enddo
           aTPXS%S0(i) = temp
           temp = 0.
           do j = 1, n
              temp = temp + md%No(j) * zTPXS(j)%S1(i)
           enddo
           aTPXS%S1(i) = temp
           temp = 0.
           do j = 1, n
              temp = temp + md%No(j) * zTPXS(j)%S2(i) 
           enddo
           aTPXS%S2(i) = temp
        enddo
        call cPrepIntpTPXS( aTPXS)
        
        do i = 1, nE
           EK = KEele(i)
           call cTPXS(aTPXS, 1, EK, S0, S1, S2)  ! S0... not in log
           call  cMWscatFixAB(S0, S1, S2, 
     *       MWavemu, MWavemu2,  
     *       aTPXS%mu(i),  aTPXS%mu2(i),  
     *       aTPXS%A0(i), aTPXS%A(i), aTPXS%B(i))
        enddo

        deallocate(zTPXS)
!          A0,A,B  cubic spline coef.
!       !!!
!                   A0 etc not in log
        call kcsplCoef(logKEele, aTPXS%A0, nE, aTPXS%coefA0, nE-1)

        call kcsplCoef(logKEele, aTPXS%A, nE, aTPXS%coefA, nE-1)

        call kcsplCoef(logKEele, aTPXS%B, nE, aTPXS%coefB, nE-1)
!            <mu>, <mu2>   cubic spline coef.
        aTPXS%logmu(:) =log( aTPXS%mu(:))
        call kcsplCoef(logKEele, aTPXS%logmu, nE, aTPXS%coefmu, nE-1)
     
        aTPXS%logmu2(:) =log( aTPXS%mu2(:))
        call kcsplCoef(logKEele, aTPXS%logmu2, nE, aTPXS%coefmu2, nE-1)

      end subroutine cMWConstForMedia

      subroutine cMWscatGetAB(name, Ein, mu, mu2, A0, A, Bp) 
        use modTPXS
        implicit none
        type(TPXSconst),intent(in):: name
        real(8),intent(in):: Ein !  e+ or e- energy in eV
        real(8),intent(out):: mu, mu2 !  <mu>, <mu2> at this Ein.
         real(8),intent(out)::A0, A, Bp  
          ! Bp is negative in region II. make is >0 when use it.
        integer:: n
        real(8):: logKE
        logKE =log( Ein )
        n = name%n
        call kcsplIntp(logKEele, name%logmu, n, name%coefmu,
     *        n-1, logKE, mu)
        mu = exp(mu)

        call kcsplIntp(logKEele, name%logmu2, n, name%coefmu2,
     *        n-1,logKE, mu2)
        mu2 = exp(mu2)
        call kcsplIntp(logKEele, name%A0, n, name%coefA0,
     *        n-1, logKE, A0)
        call kcsplIntp(logKEele, name%A, n, name%coefA, n-1, logKE, A)
        call kcsplIntp(logKEele, name%B, n, name%coefB, n-1, logKE, Bp)
      end subroutine cMWscatGetAB
