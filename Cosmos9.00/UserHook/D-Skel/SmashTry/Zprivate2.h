      integer Maxp, Maxob, MaxPtclPerCpu, MaxCPU
      parameter (Maxp=10000000, Maxob=100000)
      parameter (MaxCPU = 1000)
      parameter (MaxPtclPerCpu=Maxp/MaxCPU*4/3)
      integer basefn  ! file #; 21 22,....
      parameter (basefn=20)
      real*8  erg(Maxp), averg
      real*8  sumergi(MaxCPU), sumergw(MaxCPU),  cpupw(MaxCPU)
      integer nOnCpu(MaxCPU)
      integer idxlist(MaxPtclPerCpu, MaxCPU)
      integer idx(Maxp), idxlocal(MaxCPU)
      integer numba(MaxCPU) 
      integer ctc, Ncpu
      type(track):: ct(Maxp)
      type(parent):: pp 
      type(ob):: oo(Maxob) 
      character*128 msg
      common /arrangec/ ct, pp, oo,  sumergi, sumergw, cpupw,
     * erg,  averg, idxlist, idx, numba, nOnCpu, ctc, Ncpu
      common /arrangecc/ msg


