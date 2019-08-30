*cmz :  3.14/16 13/03/89  14.48.40  by  nick van eijndhoven (cern)
*-- author :
      subroutine normal(ran)
c
c *** nve 14-apr-1988 cern geneva ***
c
c origin : h.fesefeldt (27-oct-1983)
c
      dimension rndm(12)
      ran=-6.
      call grndm(rndm,12)
      do 1 i=1,12
      ran=ran+rndm(i)
 1    continue
      return
      end
