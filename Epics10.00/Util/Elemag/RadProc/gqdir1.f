!     ****************************************************************
!     *                                                              *
!     *  gqdir1:  integral of f(y)/sqrt(b-y)  from a to b ( b > a )  *
!     *  gqdir2:  integral of f(y)/sqrt(y-a)  from a to b ( b > a )  *
!     *                                                              *
!     ****************************************************************
!
!   /usage/
!                  call  gqdir1(f, a, b, m, ans)
!                  call  gqdir2(f, a, b, m, ans)
!
!    f:  function name with one argument
!    a:  lower boundary
!    b:  upper boundary
!    m:  m=1 specifies 10-point formula to be used
!        m=2 or other // 20 //
!  ans:  obtained integral value
!
!      *** note ***   form p.889 of abramobitz & segn
!                     if b < a, message printed, and ans=0 results.
!
!
!
      subroutine gqdir1(f, a, b, m, ans)
!
!
      dimension w(30), w2(20),  x(30), x2(20)
!
      equivalence (w(11), w2(1)),  (x(11), x2(1))
!
!                  data for 10 point formula
!
      integer i

      data (x(i),i=1,10)/
     *     .99414369,.948113606,.860343759,.739014906,.595435715,
     *    .442988686, .295882707, .167828348, .070758124, .013695586/
!
      data (w(i),i=1,10)/
     * .305506774,  .298345972, .284192218, .263377276, .236389064,
     *  .203860239,  .166553483,  .125344096,  .08120286, .03522804/
!
!                  data for 20 point formula
!
      data w2/
     a.1550118959, .1540796363, .1522207238, .1494463381, .1457731648,
     b.1412232947, .1358240916, .1296080269, .1226124849, .1148795382,
     c.1064556940, 9.739161528e-2, 8.77418164e-2, 7.7564336e-2,
     d 6.692039056e-2, 5.587401396e-2, 4.449169838e-2, 3.284211676e-2,
     e 2.099656906e-2, 9.042554196e-3/
!
      data x2/
     i.9984966997, .9865244886, .9628676424, .9280944057, .8830400419,
     j.8287867701, .7666377700, .6980858785, .6247777323, .5484742141,
     k.4710081557, .3942403130, .3200146689, .2501141389, .1862177425,
     l.1298602647, 8.23953675e-2, 4.49629902e-2, 1.846152e-2,
     m 3.521475e-3/
      j=0
!
   10 continue
      bma=b-a
      s=0.
      if(bma .le. 0.) goto 900
!
!         select n  points
!
      n=20
!
      if(m .eq. 1) n=10
      indx=n-10
!
       do   i=1,n
      if(j .ne. 0) goto 20
      y=bma*x(i+indx)+a
      goto 25
   20 continue
      y=-bma*x(i+indx)+b
   25 continue
      s=s+w(i+indx)*f(y)
       enddo
      ans=sqrt(bma)*s
      return
!
!
      entry gqdir2(f, a, b, m, ans)
!
      j=-1
      goto 10
!
  900 continue
      ans=0.
      if(bma .eq. 0.) return
      write(6,910) a,b
  910 format('0****err input to gqdir1(2): a,b=',2g9.3)
      return
      end
e(6,910) a,b
  910 format('0****err input to gqdir1(2): a,b=',2g9.3)
      return
      end
nd
