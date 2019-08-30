!            test cbremLMPXsec  (95/08/19)
!         
!      implicit none      
!      integer j, i
!      real*8  rho, e, xs 
!
!      do  j=1, 8
!          rho=3.16217 * 10.**(-(j-1)*0.75)
!          do  i=1, 16
!             e=2.e6 * 10.**( (i-1)/2.0)
!             call cbremLMPXsec(e, rho, xs)
!            write(*, *) e, xs
!        enddo
!         write(*, *) 
!      enddo
!      end
!     ********************** brems x-seciton with Landau effect
!                            for the atmosphere. ****************
!            brems total x-section/r.l for air when the landau effect is
!       taken into account at very high energy.
!       gamma energy cut is Eg/Ee = 10^-4
!       This one uses the fact that the emissin probability of 
!       fractional energy v by electron of energy E at air density
!       rho  f(v, E, rho)dv scales as f(v, E*rho)dv
!       The integral of f from vmin=10^-4 to 1 is approximated as
!       a funciton of X = E*rho.
!       
        subroutine cbremLPMXsec(e, rhoin, xs)
        implicit none
!
!  input.   e: electron/positron energy. in GeV. must be  1.e6<e<1.e13
!           i.e.,  between 10**15 to 10**22 eV.
!  input. rhoin: air density in kg/m**3. must be in 1.3 to 3.162x10^-6
!         (approx from  0km to 90 km above sea level.
!
!  output.   xs: brems probability /r.l (number of occurence of brems)
!
!       Note: Input out side of the mensioned range may be safely
!             accepted.
!

        real*8   e, rhoin, xs
        real*8   x, xxx(6), xlog
        integer i
          data ( xxx(i), i=  1,   6)/
     1       2.6807884    , -0.12901903    ,  0.27612910E-01,
     2  -0.30710641E-02,  0.86706897E-04, -0.67901808E-06               
     * /       
        x = e * rhoin * 1.e-3  ! GeV g/cm3
        if(x .lt. 30.) then
           xs = 11.447
        elseif(x .lt. 1.e10) then
           xs = 0.
           xlog = log(x)
           do i = 6, 2, -1
              xs = (xs + xxx(i))* xlog 
           enddo
           xs = xs  + xxx(1)
           xs = exp(xs)
       else
           xs =0.04198*(x/1.e10)**(-0.5)
       endif
       end
