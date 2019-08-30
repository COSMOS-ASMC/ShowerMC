*cmz :  3.14/16 27/09/90  12.11.12  by  nick van eijndhoven (cern)
*-- author :
      subroutine defs1(i,j,k)
c
c *** nve 16-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (15-jan-1984)
c
      common /vecuty/ pv(10,200)
c
c
      data pi/3.1415927/
      px=pv(1,i)
      py=pv(2,i)
      pz=pv(3,i)
      if(abs(pv(1,j)).lt.1.e-6.and.abs(pv(2,j)).lt.1.e-6) goto 1
      call lengtx(j,p)
      cost=pv(3,j)/p
      sint=sqrt(abs(1.-cost*cost))
      ph=pi/2.
      if(pv(2,j).lt.0.) ph=3.*pi/2.
      if(abs(pv(1,j)).gt.1.e-6) ph=atan2(pv(2,j),pv(1,j))
      cosp=cos(ph)
      sinp=sin(ph)
      pv(1,k)= cost*cosp*px-     sinp*py+sint*cosp*pz
      pv(2,k)= cost*sinp*px+     cosp*py+sint*sinp*pz
      pv(3,k)=-sint     *px             +cost     *pz
      return
    1 pv(1,k)=px
      pv(2,k)=py
      pv(3,k)=pz
c --- take the case of theta=pi into account (mr/nve 27-sep-1990) ---
      if (pv(3,j) .lt. 0.) pv(3,k)=-pz
      return
      end
