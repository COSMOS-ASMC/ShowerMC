module modTPXS
    implicit none
!    next 3:    used when reading TPXS
    character(70):: line
    character(200):: filename
    character(1):: field(1)
!
    integer,parameter:: nEneg=106  ! # of energy grid  for e-
    integer,parameter:: nEpos=106  ! # of energy grid  for e+
    real(8),parameter:: eps= 1.d-5
      !  share energy for e- and e+
    real(8),save:: KEele(nEneg) = 0. ! nEneg >= nEpos assumed
    real(8),save:: logKEele(nEneg) = 0. ! nEneg >= nEpos assumed
    real(8),save::maxKEpos=50e-3   ! K.E GeV (100MeV) abvoe this
               !  e+ is same as e-. not used yet

    type TPXSconst
!                  total xs, 1st 2nd tpxs 
       real(8):: S0(nEneg), S1(nEneg), S2(nEneg)
       real(8):: logS0(nEneg), logS1(nEneg), logS2(nEneg)
!       cubic spline coef. for logS0 ....
       real(8):: coefS0(nEneg-1, 3), coefS1(nEneg-1,3), coefS2(nEneg-1,3) 
       !  <mu> <mu2>
       real(8):: mu(nEneg), mu2(nEneg)
       real(8):: logmu(nEneg), logmu2(nEneg)
!           for  log value coef.
       real(8):: coefmu(nEneg-1, 3), coefmu2(nEneg-1,3)
!         used for MW model
       real(8):: A0(nEneg), A(nEneg), B(nEneg)
       real(8):: coefA0(nEneg-1, 2), coefA(nEneg-1,3), coefB(nEneg-1,3) 
       integer::n=0   
    end type TPXSconst
  end module modTPXS
