*cmz :  3.14/16 13/03/89  14.48.42  by  nick van eijndhoven (cern)
*-- author :
      real function fctcos(t)
c
c *** nve 01-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (27-oct-1983)
c
      common/coscom/aa,bb,cc,dd,rr
c
c --- boundary limits for arguments of intrinsic functions ---
c --- xl denotes lower bound whereas xu denotes upper bound ---
      common /limits/ expxl,expxu
c
c
      double precision test1,test2
c
      test1=-bb*t*1.0d0
      if (test1 .gt. expxu) test1=expxu
      if (test1 .lt. expxl) test1=expxl
      test2=-dd*t*1.0d0
      if (test2 .gt. expxu) test2=expxu
      if (test2 .lt. expxl) test2=expxl
c
      fctcos=aa*exp(test1)+cc*exp(test2)-rr
c
      return
      end
