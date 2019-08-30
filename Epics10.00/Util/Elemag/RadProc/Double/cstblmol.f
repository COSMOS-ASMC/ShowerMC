!     ****************************************************************
!     *                                                              *
!     * cstblmol:  create sampling table for moliere scattering      *
!     *                                                              *
!     ****************************************************************
!
!
!                  set jpunch(^=0 --> output in disk)
!
!  **** note ****
!        no z-dependent part.  b is a function of energy and z but
!        is parameter here so that it must be given in sampling routine
!
!
!
!
!        rutbnr, gquadt,expi,fmol2t,fmoli,fmolr,mkdt needed
!
      real ubscat( 7,7),ubsca2(11,7)
      external fmolr
      common/cmol/cnstu,b
      data er/0.0001/
      jpunch=0
      dt=0.2
      b=4.5
      db=2.
       do   ib=1,7
      t=0.
       do   i=1,35
      x=t*t
      f=fmol(x)
      f2=fmol2t(x)
      f3=fmol3t(x)
      write(6,50) t,x,f,f2,f3
   50 format(5g15.4)
      t=t+dt
       enddo
      cnst=fmoli(100.)
      write(6,5) b,cnst
    5 format(1h0,2g15.5)
      u=0.6
      du=0.05
      x=0.9
      xmin=0.
       do   iu=1,7
      cnstu=cnst*u
      call rutbnr(fmolr,x,er,xmin,100.,10.,x,j)
      if(j .gt. 0) goto 180
      write(6,175) b,u,x
  175 format(3g15.5)
  180 continue
      xmin=amax1(x-0.1,0.)
      ubscat(iu,ib)=x
      u=u+du
       enddo
      dusq=sqrt(0.1)/10.
       do   iu=1,10
      u=1.0-(      float(11-iu)*dusq)**2
      cnstu=cnst*u
      call rutbnr(fmolr,x,er,xmin,100.,10.,x,j)
      if(j .gt.0) goto 280
      write(6,175) b,u,x
  280 continue
      xmin=amax1(x-0.1,0.)
      ubsca2(iu,ib)=alog(x)
       enddo
      ubsca2(11,ib)=alog(100.)
      b=b+db
       enddo
      call mkdt('ubscat',ubscat,1,49,0.,0,jpunch)
      call mkdt('ubsca2',ubsca2,1,77,0.,0,jpunch)
      stop
      end
      function fmol(x)
!
!        moliere function x=theta**2, including 3rrd ter.
!
      common/cmol/cnstu,b
      fmol=exp(-x)+(fmol2t(x)+fmol3t(x)/b)/b
      return
      end
      function fmol2t(x)
      dimension t1(11),t2(11),t3(5)
      data t1/ 2.4929, 2.0694, 1.0488, -0.0044, -0.6068, -0.6359,
     1  -0.3086, 0.0525, 0.2423, 0.2386, 0.1316/,
     2  t2/1.053, 0.475, -0.775, -1.483, -1.676, -1.448, -1.008,
     3 -0.566, -0.222, 0.005, 0.1375/
      data t3/0.1375, 0.2521, 0.2602, 0.2456, 0.2264/
!        t1:  f2 from o to 2 step 0.2
!        t2: f2/2*theta**4 from 2 to 4 step 0.2
!        t3: f2/2*thta**4 from 4 to 6 step 0.5
!
      expx=exp(-x)
      fmol2t=(x-1.)*eiexp(x)-1.+expx+expx
      return
!
!     ************
      entry fmol3t(x)
!     ************
!
!        3rd term/2 of moliere expansion.  x=theta**2
      t=sqrt(x)
      if(t .gt. 2.) goto 100
!
      fmol2t=fintp4(t1,0.,0.2,11,t)/2.
      return
!
  100 continue
      if(t .gt.4.) goto 110
      tmp=fintp4(t2,2.,0.2,11,t)
  105 continue
      fmol2t=tmp/t**4
      return
!
  110 continue
      if(t .gt.5.5) goto 120
      tmp=fintp4(t3,4.,0.5,5,t)
      goto 105
!
  120 continue
      tmp=0.5971*exp(-0.1616*t)
      goto 105
      end
      function fmoli(x)
      external fmol
      data xsv/0./,sisv/0./
      si=sisv
      x2=xsv
      if(xsv .gt. x) go to50
   15 x1=x2
      dx=5.
      if(x1 .lt. 10.) dx=1.
      x2=x1+dx
      x2=amin1(x2,x)
      call gquadt(fmol,x1,x2,ans)
      si=si+ans
      if(x2 .ne. x) go to 15
   20 fmoli=si
      xsv=x
      sisv=si
      return
   50 x1=x2
      dx=5.
      if(x1 .lt. 20.) dx=1.
      x2=x1-dx
      x2=amax1(x2,x)
      call gquadt(fmol,x2,x1,ans)
      si=si-ans
      if(x2 .ne. x) go to 50
      go to 20
      end
      function fmolr(x)
      common/cmol/cnstu,b
      tmp=fmoli(x)
      fmolr=tmp/cnstu-1.
      return
      end
to 50
      go to 20
      end
      function fmolr(x)
      common/cmol/cnstu,b
      tmp=fmoli(x)
      fmolr=tmp/cnstu-1.
      return
      end
