      subroutine crecprob(depthin, code, limit, dfai, E0in, NN, cosz, 
     *  nr, recprob, nptcl, age,  sum, Nx)
      implicit none
!
!        compute probabilty that we accept observed particles
!        which fall in at Moliere unit r and r+dr with dfai
!        azimuthal angle range; the probabilty is determined 
!        so that the maximumk number may not exceed 'limit'
!        in each region.  
!        This is based on 10^18 eV shower. 
!
      real depthin   ! input g/cm^2
      integer code ! input gamma-1, e-2, mu-3, had-4
      real limit   ! input max number of particls to be recorded
                   ! in actual simulation run.
      integer nr   ! input number of r bin; 42. is standard
                   ! r is digitized from r1=0.01 (m.u) with 0.1 log10
                   ! bin.
      real dfai    ! input  fai angle (deg) bin. 30 is standard
      real E0in      ! input  primary total energy in GeV.
      integer NN   ! input  primary nucleon number.
      real cosz    ! input  primary zenith angle in cosine
      real recprob(nr) ! ouput. probability to record particles
                       ! at r (at depth) in dfai region. if >1.
                       ! accept always
      real nptcl(nr)   ! output. average number of particles falling
                       ! on the given band region.
      real age        ! output. given depth is at age

      real sum  !  output sum of particles
      real Nx    !  output number of particles at depth (Nx ~ sum)

      integer idep, ir
      real r1, r2,  A,  s, r
      real sizex, asdensity, rho,  n
      real pi, fita, fitb, fitc
      real  eps, depth, E0, fac
      real cogdep, dep2cogdep, cogdep2s
!
      pi=3.141592
      sum = 0.
      if(NN .ge. 1) then
         E0=E0in/NN
         fac = NN
      else
!         to be updateded
         E0 = E0in
         fac = 1
      endif 
      depth = depthin/cosz
      cogdep = dep2cogdep(E0in, NN, depth)
      s = cogdep2s(E0in, NN, cogdep)
      Nx = sizex(E0in, NN, code, s)
      r1 = 0.01/10**0.05
      do ir=1, nr
         r2=r1*10.**0.1
         A = pi*(r2**2- r1**2)*dfai/360.
         r = sqrt(r1*r2)
         rho=asdensity(E0, code, s, r)
         n= Nx*rho*A
         nptcl(ir)=n
         sum = sum + n
         eps =limit/n
         recprob(ir)= eps
         r1= r2
      enddo
      age = s
      end
