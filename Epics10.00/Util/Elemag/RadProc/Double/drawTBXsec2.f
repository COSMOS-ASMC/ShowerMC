!          this is to use the fact f(v, Ee, rho)dv = f(v, Ee * rho)dv
!        Let  Ee * rho be X; since we consider
!        Ee=10^6 GeV to 10^13 GeV
!        rho = 10^-3 to 10^-9 g/cm3.  
!        the range of X is 10^-3 to 10^10 GeV g/cm3.
!       In this calculation, we fix the energy to be Ee = 10^10 GeV (10^19 eV)
!       and variate the rho so that X is in the range mentioned above.
!       
      implicit none
      real*8 rho,  z, a,  ans
      real*8 e, X
      common /landu1/ e
      real*8 v(300)
      integer nx
      external fff

      real*8 vmin
 

      z = 7.25
      a = 14.5


      nx = 0
      X = 10.**(-3)
      e = 1.e19/1.e9
      do while (X .lt. 1.3e10)
         nx = nx + 1
         rho = X/e
         call zpart( z, a,  rho)
!         vmin=max(1.d-3/e, 1.d-5)
         call k16pGaussLeg(fff, 0.d0, 1.d0, 16, ans)
         write(*,*) sngl(X), sngl(ans)
         v(nx) = log10(ans)
         X = X * 10.**0.1
!         call kmkDataStm2b(v, ne, nrho, 'bremxs', 'f9.4',9)
      enddo
      end
      real*8 function fff(v)
      real*8 v



      real*8 fbrem

      fff = v *  fbrem(v)
      end
