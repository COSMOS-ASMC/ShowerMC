*cmz :  3.14/16 13/03/89  14.48.41  by  nick van eijndhoven (cern)
*-- author :
      function pmltpc(np,nm,nz,n,b,c)
c
c *** nve 03-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (14-sep-1987)
c
      rlnnpf=0.
      if(np.le.1) goto 2
      do 1 i=2,np
    1 rlnnpf=rlnnpf+log(i*1.)
    2 rlnnmf=0.
      if(nm.le.1) goto 4
      do 3 i=2,nm
    3 rlnnmf=rlnnmf+log(i*1.)
    4 rlnnzf=0.
      if(nz.le.1) goto 6
      do 5 i=2,nz
    5 rlnnzf=rlnnzf+log(i*1.)
    6 pmltpc=-(np-nm+nz+b)**2/(2*(c*n)**2)-rlnnpf-rlnnmf-rlnnzf
      if(pmltpc.lt.-50.) goto 10
      pmltpc=exp(pmltpc)
      return
   10 pmltpc=0.
      return
      end
