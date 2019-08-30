       integer NdepthCompM     ! max number of layers of comp. shower.
       integer Ncs
       parameter (NdepthCompM = 16, Ncs=420)

       integer NoOfCompS  !  number of component showers
       integer NdepthComp ! actual number of layers for component showers
       real*4  DepthComp(NdepthCompM)  ! depths where c.s sizes are given 
       real*4  ElecSizeComp(Ncs, NdepthCompM)   ! size of c.s at each depth
                                           ! in log10
       real*4  ElecAgeComp(Ncs, NdepthCompM)   ! age
       real*4  E0Comp     ! component shower primary energy.
       real*4  FirstColDep ! in kg/m**2. First collision depth

       common /ZASbyH/
     1 DepthComp, ElecSizeComp, ElecAgeComp,  E0Comp, FirstColDep,
     2 NoOfCompS, NdepthComp


