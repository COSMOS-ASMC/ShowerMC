c***********************************************************************
c...Include file of JAM
c***********************************************************************

C...Double precision.
      implicit double precision(a-h, o-z)

c...Maximum number of all particles.
c     parameter(mxv=800)   ! AGS Au+Au
c     parameter(mxv=2000)
      parameter(mxv=30000)
c     parameter(mxv=100000)

c...The event record.
      common/jamevnt1/r(5,mxv),p(5,mxv),v(5,mxv),k(11,mxv)
      common/jamevnt2/nv,nbary,nmeson
      save /jamevnt1/,/jamevnt2/

c...Parton properites or decay orientation.
      common/jamjet1/vq(10,mxv),kq(2,mxv)
      save /jamjet1/

c...Deuteron coalesence.
      common/jamevnt3/weid(mxv)
      save /jamevnt3/

c...Save two-body collisision information.
      common/jamsave1/rcp(5,6),pcp(5,6),vcp(5,6),kcp(11,6)
      common/jamsave2/vqcp(10,6),kqcp(2,6)
      save /jamsave1/,/jamsave2/

c...Max. collision number in each time step.
      parameter (mxcoll=500000)
c     parameter (mxcoll=50000)

c...Collision matrix.
      common/jamcol/coll(6,mxcoll),icoll(3,mxcoll),mentry
      save /jamcol/

c...Box by maru  cell.
      common/jamcell/rcell(3,27)
      save /jamcell/

c=======================================================================
