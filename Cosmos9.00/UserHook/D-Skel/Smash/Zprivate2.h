      integer Maxp, Maxob, MaxPtclPerCpu, MaxCPU
      parameter (Maxp=1200000, Maxob=60000)
!      for p10^17 eV Emin=6e3, Maxp~280000 Maxob~23000
!      parameter (MaxCPU =41960)   ! must be >= NCPU
      parameter (MaxCPU =1280)   ! must be >= NCPU
      parameter (MaxPtclPerCpu=Maxp/MaxCPU*4/3)
      integer basefn  ! file #; 21 22,....
      parameter (basefn=20)
      real*8  erg(Maxp), averg
      real*8  sumergi(MaxCPU), sumergw(MaxCPU),  cpupw(MaxCPU)
      integer nOnCpu(MaxCPU)
      integer idxlist(MaxPtclPerCpu, MaxCPU)
      integer idx(Maxp), idxlocal(MaxCPU)
      integer numba(MaxCPU) 
      integer ctc, Ncpu,  Mcpu, Margin
      type(track):: ct(Maxp)
      type(parent):: pp 
      type(ob):: oo(Maxob) 
      character*128 msg
      character*120 skelfile(MaxCPU)
      common /arrangec/ ct, pp, oo,  sumergi, sumergw, cpupw,
     * erg,  averg, idxlist, idx, numba, nOnCpu, ctc, Ncpu,
     * Mcpu, Margin
      common /arrangecc/ msg, skelfile


