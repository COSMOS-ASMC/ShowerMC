*cmz :  3.14/16 28/09/90  10.09.25  by  nick van eijndhoven (cern)
*-- author :
      subroutine selfab(sprob)
c
c *** self-absorbtion in heavy molecules ***
c *** nve 16-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (11-oct-1987)
c
      common /vecuty/ pv(10,200)
c
c
      sprob=0.
      ekw=pv(4,200)-abs(pv(5,200))
      if(ekw.lt.5.) return
      alekw=log(ekw-4.)
      sprob=0.6*alekw
      if(sprob.gt.1.) sprob=1.
      return
      end
