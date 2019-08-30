*cmz :  3.14/16 13/03/89  14.48.39  by  nick van eijndhoven (cern)
*-- author :
      function ranres(x)
c
c *** restricted random numbers to be used as argument in log etc... ***
c *** nve 13-jul-1988 cern geneva ***
c
c note : 0 < ranres < 1
      dimension rndm(1)
c
 1    continue
      call grndm(rndm,1)
      ranres=rndm(1)
      if ((ranres .le. 0.) .or. (ranres .ge. 1.)) go to 1
      return
      end
