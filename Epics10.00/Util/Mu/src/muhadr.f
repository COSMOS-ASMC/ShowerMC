!         muhadr: init
!         muhadv: sample fractional energy loss by muon nuc. int (v>eps)
!         muhadt: sample path for muon nuc. int. (g/cm**2)
!         muhadf: fractional energy loss by muon nuc. int.(v<eps)/(g/cm2
      subroutine muhadr(eps)
!           specify eps below which the energy loss is treated as
!           continuous one. (eps=1/100,1/30,1/20,1/10 may be given for
!           seeing energy loss and sampling
!          ( eps=1. may be given when to see b-term)
         common /$muhad/ ieps, epsi
         integer iepsa(5)/100, 30, 20, 10, 1/
         if(abs(eps-0.01)/eps .lt. 1.e-2) then
             ieps=1
         elseif(abs(eps-.0333333)/eps .lt. 1.e-2) then
             ieps=2
         elseif(abs(eps-.05)/eps .lt. 1.e-2) then
             ieps=3
         elseif(abs(eps-.1)/eps .lt. 1.e-2) then
             ieps=4
         elseif(abs(eps-1.)/eps .lt. 1.e-2) then
             ieps=5
         else
             write(*,*) ' eps=',eps, ' invalid for mupair'
             stop
         endif
         epsi=1./iepsa(ieps)
         write(*,*) ' ----had-p by muon: taken into account for',
     *              ' v > 1/',  iepsa(ieps), ' ---'
         write(*,*) '     lower energy had-p is treated as a part of ',
     *              'continuous energy loss.'
      end
!c         test muhadf .
!*include mucset
!      integer ia(5)/100, 30, 20, 10, 1/
!      character txt*70, capx*16, capy*16
!      open(13, file='c2s5001.#gd.data')
!      capx='e/tev'
!      capy='b(/10**6g/cm2)'
!      call mucset(11.,22., .5, 5.5)
!      do 100 ieps=1, 5
!          write(txt,'('' eps=1/'',i3)') ia(ieps)
!          write(13) txt
!          write(13) capx, capy
!          call muhadr(1./ia(ieps))
!          do 50 j=1, 1000
!              e=1.e-3* 10.**( (j-1)/20.)
!              if(e .gt. 1000.) goto 55
!              call muhadf(e, b)
!              write(13) e, b/1.e-6
!  50      continue
!  55      continue
!          write(13) 1.e50, 1.e50
! 100  continue
!      end
!        fractional b-term for muon hadron production.
!       dedtbe: output. fractional loss /(g/cm**2) (by nuc. int v<eps)
       subroutine muhadf(e, dedtbe)
       include  'Zmucom.f'
           call muhb10(e, b1)
           call muhb20(b2)
           call gpxs(e, xs)
           dedtbe=(b1*pointl+ (1.-pointl)*shadow*b2)
     *            * xs*gppg
       end
!c         test muhb10 .
!      integer ia(5)/100, 30, 20, 10, 1/
!      character txt*70, capx*16, capy*16
!      open(13, file='c2s5001.#gd.data')
!      capx='e/tev'
!      capy='muhb10'
!      do 100 ieps=1, 5
!          write(txt,'('' eps=1/'',i3)') ia(ieps)
!          write(13) txt
!          write(13) capx, capy
!          call muhadr(1./ia(ieps))
!          do 50 j=1, 1000
!              e=1.e-3* 10.**( (j-1)/20.)
!              if(e .gt. 1000.) goto 55
!              call muhb10(e, b)
!              write(13) e, b
!  50      continue
!  55      continue
!          write(13) 1.e50, 1.e50
! 100  continue
!      end
!        fractional b-term for muon hadron production.
!        point like part.
!        if b is  multiplied  by alfa*n/pi * sigma(gp)*
!        beta  (beta=point like fraction),
!        this gives -de/dt/e( v< eps) (/(g/cm**2))
       subroutine muhb10(e, b)
         common /$muhad/ ieps, epsi
!       integral 0-eps of   fai(inf) of had. prod. by mu
!       x=log10(e/tev)
!           eps=1/100
!       xmin=   -2.101     xmax=    3.051
      plsq1(x)=0.2288808e-01*x+0.1652847
!          eps=1/30
!       xmin=   -2.101     xmax=    3.051
      plsq2(x)=0.7546479e-01*x+0.5054637
!
!       xmin=   -2.101     xmax=    3.051
      plsq3(x)=0.1122759    *x+0.7320931
!
!       xmin=   -2.101     xmax=    3.051
      plsq4(x)=0.2191063    *x+ 1.361674
!
!       xmin=   -2.101     xmax=    3.051
      plsq5(x)= 1.535032    *x+ 7.771255
!
      x=log10(e)
      if(ieps .eq. 1) then
          b=plsq1(x)
      elseif(ieps .eq. 2) then
          b=plsq2(x)
      elseif(ieps .eq. 3) then
          b=plsq3(x)
      elseif(ieps .eq. 4) then
          b=plsq4(x)
      elseif(ieps .eq. 5) then
          b=plsq5(x)
      else
          write(*,*) ' ieps=',ieps,' undef for muhb10'
          stop
      endif
      end
!c         test muhx10 .
!      integer ia(4)/100, 30, 20, 10/
!      character txt*70, capx*16, capy*16
!      open(13, file='c2s5001.#gd.data')
!      capx='e/tev'
!      capy='muhx10'
!      do 100 ieps=1, 4
!          write(txt,'('' eps=1/'',i3)') ia(ieps)
!          write(13) txt
!          write(13) capx, capy
!          call muhadr(1./ia(ieps))
!          do 50 j=1, 1000
!              e=1.e-3* 10.**( (j-1)/20.)
!              if(e .gt. 1000.) goto 55
!              call muhx10(e, b)
!              write(13) e, b
!  50      continue
!  55      continue
!          write(13) 1.e50, 1.e50
! 100  continue
!      end
       subroutine muhx10(e, tx)
         common /$muhad/ ieps, epsi
!   integral eps-1 of fai(inf)/v of had. prod. by mu
!          x=log10(e/tev)
!          eps=1/100
!       xmin=   -2.101     xmax=    3.051
      plsq1(x)= 8.899820    *x+ 51.66591
!          eps=1/30
!       xmin=   -2.101     xmax=    3.051
      plsq2(x)= 6.180719    *x+ 33.93083
!          eps=1/20
!       xmin=   -2.101     xmax=    3.051
      plsq3(x)= 5.284678    *x+ 28.40878
!          eps=1/10
!       xmin=   -2.101     xmax=    3.051
      plsq4(x)= 3.799458    *x+ 19.62849
      x=log10(e)
      if(ieps .eq. 1) then
          tx=plsq1(x)
      elseif(ieps .eq. 2) then
          tx=plsq2(x)
      elseif(ieps .eq. 3) then
          tx=plsq3(x)
      elseif(ieps .eq. 4) then
          tx=plsq4(x)
      elseif(ieps .eq. 5) then
          tx=0.
      else
          write(*,*) ' ieps=',ieps,' undef for muhx10'
          stop
      endif
      end
!        fractional b-term for muon hadron production.
!        hadron like part.
!        if b is  multiplied  by alfa*n/pi * sigma(gp)*
!        (1-beta) a**-0.1 (beta=point like fraction,
!        a=matter mass #)
!        this gives -de/dt/e( v< eps) (/(g/cm**2))
      subroutine muhb20(b)
!         integral 0 to eps of fai(had).
         common /$muhad/ ieps, epsi
!       integral 0-eps of   fai(had) of had. prod. by mu
!
        real*4 ans(5)/.136, .3699, .5108, .8652, 2.956/
      if(ieps .ge.1 .and. ieps .le. 5) then
          b=ans(ieps)
      else
          write(*,*) ' ieps=',ieps,' undef for muhb20'
          stop
      endif
      end
      subroutine muhx20(tx)
         common /$muhad/ ieps, epsi
!       integral eps-1   of   fai(had)/v of had. prod. by mu
        real*4 ans(4)/28.53, 16.15, 12.71, 7.732/
      if(ieps .ge.1 .and. ieps .le. 4) then
         tx=ans(ieps)
      elseif(ieps .eq. 5) then
         tx=0.
      else
          write(*,*) ' ieps=',ieps,' undef for muhx20'
          stop
      endif
      end
!        total x-section of hadron production by muon(point like)
!        (if this is multiplied by alfa*n/pi *sigma(gp)-->/(g/cm**2)
      subroutine muhx1(e, tx)
       include  'Zmucom.f'
          call muhx10(e, tx)
          tx=tx*pointl
      end
!        total x-section of hadron production by muon(hdron like)
!        (if this is multiplied by alfa*n/pi *sigma(gp)-->/(g/cm**2)
!
      subroutine muhx2(tx)
       include  'Zmucom.f'
          call muhx20(tx)
          tx=tx*(1.-pointl)*shadow
      end
!c         test muhadx .
!*include mucset
!      integer ia(4)/100, 30, 20, 10/
!      character txt*70, capx*16, capy*16
!      open(13, file='c2s5001.#gd.data')
!      capx='e/tev'
!      capy='mfp(g/cm**2)'
!      call mucset(11.,22., .5, 5.5)
!      do 100 ieps=1, 4
!          write(txt,'('' had: eps=1/'',i3)') ia(ieps)
!          write(13) txt
!          write(13) capx, capy
!          call muhadr(1./ia(ieps))
!          do 50 j=1, 1000
!              e=1.e-3* 10.**( (j-1)/20.)
!              if(e .gt. 1000.) goto 55
!              call muhadx(e, tx)
!              write(13) e, 1./tx
!  50      continue
!  55      continue
!          write(13) 1.e50, 1.e50
! 100  continue
!      end
!        total x-section of hadron production by muon
!          /(g/cm**2)
      subroutine muhadx(e, tx)
       include  'Zmucom.f'
          call muhx1(e, tx1)
          call muhx2(tx2)
          call gpxs(e, sig)
          tx=gppg*sig*(tx1+tx2)
      end
!c         test gpxs
!     character*70 ttl
!     character*16 capx,capy
!       f(e)=1.647*log( 2.*.938*e*1000./88. )**2+ 114.3
!     capx='eg_lab(tev)'
!     capy='x-sec(micro_b)'
!     ttl='new version (e**0.073)'
!     open(13, file='c2s5001.#gd.data')
!     write(13) ttl
!     write(13) capx, capy
!     do 100 i=1,  600
!         e=.3e-3*10.**( (i-1)/20. )
!         if(e .gt. 5000.) goto 101
!         call gpxs(e, sig)
!         write(13) e,sig
! 100 continue
! 101 continue
!     write(13) 1.e50, 1.e50
!     ttl='old version(before 90.05.24). e**0.11'
!     write(13) ttl
!     write(13) capx, capy
!     do 200 i=1,  600
!         e=.3e-3*10.**( (i-1)/20. )
!         if(e .gt. 5000.) goto 201
!         sig=f(e)
!         write(13) e,sig
! 200 continue
! 201 continue
!     write(13) 1.e50, 1.e50
!     end
      subroutine gpxs(e, sig)
!          this gives g-p x-section in micro barn
!          e in tev
!             for the 2nd version of the sampling, the following is
!             modified to realize that the dependence at high
!             energies is the same as p-p_bar which is
!             e_lab**0.073.
        f(e)=1.647*log( 2.*.938*e*1000./88. )**2+ 114.3
!
        if(e .gt. .5) then
            sig=f(.5)*(e/.5)**0.073
        else
            sig=f(e)
        endif
      end
!           test  muhv2
!         common /$muhad/ ieps, epsi
!      open(13,file='c2s5001.#gd.data')
!      do 10 ieps=1, 4
!          do 100 j=1, 50000
!            call muhv2(v)
!            write(13) v
! 100      continue
!          write(13) 1.e50, 1.e50
!  10  continue
!      end
!        u-v/eps**u for hadron like part of muon hadron production
        subroutine muhv2(v)
          common /$muhad/ ieps, epsi
          if(ieps .eq. 1) then
             call muhv21(v)
          elseif(ieps .eq. 2) then
             call muhv22(v)
          elseif(ieps .eq. 3) then
             call muhv23(v)
          elseif(ieps .eq. 4) then
             call muhv24(v)
          elseif(ieps .eq. 5) then
             v=1.
          else
             write(*,*) ' ieps=',ieps,' undef for muhv2'
             stop
          endif
        end
        subroutine muhv21(v)
!             eps=1/100
!          x=u, plsq=v/eps**u
!       xmin=  -0.1574e-01 xmax=   0.2004
      plsq1(x)=((((( 529822.4    *x-324028.9    )*x+ 78022.75    )*x-946
     *1.043    )*x+ 624.8511    )*x-24.20937    )*x+0.9509603
!       xmin=   0.8946e-01 xmax=    1.011
      plsq2(x)=((( 1.765125    *x-3.667789    )*x+ 3.779487    )*x-1.381
     *414    )*x+0.5083277
          data eps/0.01/
!
          call rndc(u)
          if(u .gt. .15) then
             v=eps**u * plsq2(u)
          else
             v=eps**u * plsq1(u)
          endif
       end
        subroutine muhv22(v)
!             eps=1/30
!          x=u, plsq=v/eps**u
!
          data eps/0.03333333/
!       xmin=  -0.1576e-01 xmax=   0.1945
      plsq1(x)=((((-29157.52    *x+ 14869.91    )*x-2890.506    )*x+ 276
     *.3206    )*x-14.96122    )*x+0.9602987
!
!       xmin=   0.8861e-01 xmax=    1.009
      plsq2(x)=((( 1.548149    *x-3.576385    )*x+ 3.856900    )*x-1.458
     *364    )*x+0.6343563
!
          call rndc(u)
          if(u .gt. .15) then
             v=eps**u * plsq2(u)
          else
             v=eps**u * plsq1(u)
          endif
       end
        subroutine muhv23(v)
!             eps=1/20
!          x=u, plsq=v/eps**u
!
          data eps/0.05/
!
!       xmin=  -0.1577e-01 xmax=   0.1994
      plsq1(x)=((((-22026.70    *x+ 11519.39    )*x-2295.219    )*x+ 225
     *.4575    )*x-12.74661    )*x+0.9657221
!
!       xmin=   0.9861e-01 xmax=    1.009
      plsq2(x)=((( 1.350066    *x-3.227911    )*x+ 3.594836    )*x-1.389
     *973    )*x+0.6769658
          call rndc(u)
          if(u .gt. .15) then
             v=eps**u * plsq2(u)
          else
             v=eps**u * plsq1(u)
          endif
       end
        subroutine muhv24(v)
!             eps=1/10
!          x=u, plsq=v/eps**u
!
          data eps/0.10/
!
!       xmin=  -0.1564e-01 xmax=   0.2070
      plsq1(x)=((((-11924.51    *x+ 6517.480    )*x-1358.867    )*x+ 140
     *.8766    )*x-8.744789    )*x+0.9705680
!
!       xmin=   0.1130     xmax=    1.006
      plsq2(x)=((( 1.034986    *x-2.638535    )*x+ 3.092931    )*x-1.249
     *036    )*x+0.7626207
!
          call rndc(u)
          if(u .gt. .15) then
             v=eps**u * plsq2(u)
          else
             v=eps**u * plsq1(u)
          endif
       end
!c         test muhv1
!*include  mucset
!c
!      integer ia(4)/100, 30, 20, 10/
!      character txt*70, capx*16
!      open(13, file='c2s5001.#gd.data')
!      capx='v'
!      call mucset(11.,22., .5, 5.5)
!      do 100 ieps=1, 4
!          epsi=1./ia(ieps)
!          call muhadr(epsi)
!          do 50 j=1, 1000
!              e=3.16e-2* 10.**( (j-1)/1.)
!              if(e .gt.  500.) goto 55
!              if(j .eq. 1) then
!                write(txt,'('' eps=1/'',i3,'' e='',f8.4)') ia(ieps), e
!              else
!                write(txt,'('' e='',f8.4)') e
!              endif
!              write(13) txt
!              write(13) capx
!              do 30 k=1, 20000
!                  call muhv1(e, v)
!                  write(13) v
!  30          continue
!              write(13) 1.e50
!  50      continue
!  55      continue
! 100  continue
!      end
      subroutine  muhv1(ein,v)
!        sampling table of v for muon had. point like
!        v(iu, ie):  iu=1,51 (u=0 (.020) 1)
!                    ie=1, 5 (log10(e/tev)=-2 (1.0) 2  )
!        v1= v/eps**u
      parameter (iu=51, ie=5)
!
      dimension  vtbl(iu, ie, 4)
      dimension v1(iu*ie), v2(iu*ie), v3(iu*ie), v4(iu*ie)
!
      equivalence (v1(1), vtbl(1,1,1))
      equivalence (v2(1),   vtbl(1,1,2))
      equivalence (v3(1),   vtbl(1,1,3))
      equivalence (v4(1),   vtbl(1,1,4))
!
!          eps=1/100
!        point like part sampling table.  eps=1/100
      data (v1    (i),i=   1,  63)/
     1 1.00000, 0.83846, 0.74967, 0.68786, 0.64285, 0.60913, 0.58337,
     2 0.56346, 0.54801, 0.53604, 0.52688, 0.52002, 0.51509, 0.51181,
     3 0.50996, 0.50939, 0.50994, 0.51152, 0.51404, 0.51744, 0.52167,
     4 0.52667, 0.53243, 0.53891, 0.54610, 0.55397, 0.56253, 0.57177,
     5 0.58168, 0.59227, 0.60354, 0.61550, 0.62816, 0.64153, 0.65563,
     6 0.67048, 0.68608, 0.70247, 0.71967, 0.73770, 0.75658, 0.77636,
     7 0.79705, 0.81869, 0.84131, 0.86496, 0.88967, 0.91549, 0.94245,
     8 0.97060, 1.00000, 1.00000, 0.86607, 0.78412, 0.72495, 0.68088,
     9 0.64736, 0.62148, 0.60130, 0.58551, 0.57320, 0.56369, 0.55649/
      data (v1    (i),i=  64, 126)/
     1 0.55124, 0.54765, 0.54550, 0.54463, 0.54489, 0.54618, 0.54841,
     2 0.55151, 0.55542, 0.56011, 0.56552, 0.57164, 0.57844, 0.58591,
     3 0.59402, 0.60278, 0.61218, 0.62221, 0.63288, 0.64419, 0.65615,
     4 0.66876, 0.68204, 0.69599, 0.71064, 0.72599, 0.74206, 0.75888,
     5 0.77646, 0.79483, 0.81400, 0.83402, 0.85489, 0.87665, 0.89934,
     6 0.92299, 0.94762, 0.97328, 1.00000, 1.00000, 0.88206, 0.80506,
     7 0.74815, 0.70514, 0.67211, 0.64642, 0.62629, 0.61048, 0.59810,
     8 0.58850, 0.58121, 0.57585, 0.57215, 0.56988, 0.56889, 0.56903,
     9 0.57018, 0.57227, 0.57522, 0.57897, 0.58348, 0.58870, 0.59460/
      data (v1    (i),i= 127, 189)/
     1 0.60117, 0.60838, 0.61621, 0.62466, 0.63371, 0.64337, 0.65364,
     2 0.66450, 0.67597, 0.68805, 0.70076, 0.71408, 0.72805, 0.74267,
     3 0.75794, 0.77390, 0.79055, 0.80792, 0.82602, 0.84487, 0.86450,
     4 0.88493, 0.90618, 0.92829, 0.95127, 0.97516, 1.00000, 1.00000,
     5 0.89243, 0.81909, 0.76398, 0.72191, 0.68937, 0.66396, 0.64398,
     6 0.62825, 0.61590, 0.60631, 0.59901, 0.59363, 0.58990, 0.58760,
     7 0.58656, 0.58665, 0.58775, 0.58977, 0.59263, 0.59629, 0.60069,
     8 0.60579, 0.61156, 0.61798, 0.62501, 0.63265, 0.64088, 0.64970,
     9 0.65909, 0.66906, 0.67960, 0.69072, 0.70242, 0.71469, 0.72756/
      data (v1    (i),i= 190, 252)/
     1 0.74102, 0.75509, 0.76978, 0.78510, 0.80106, 0.81768, 0.83497,
     2 0.85296, 0.87166, 0.89109, 0.91127, 0.93223, 0.95399, 0.97657,
     3 1.00000, 1.00000, 0.89970, 0.82912, 0.77545, 0.73417, 0.70209,
     4 0.67695, 0.65714, 0.64151, 0.62923, 0.61969, 0.61242, 0.60706,
     5 0.60333, 0.60104, 0.59999, 0.60006, 0.60113, 0.60312, 0.60595,
     6 0.60955, 0.61389, 0.61891, 0.62459, 0.63089, 0.63781, 0.64531,
     7 0.65338, 0.66202, 0.67122, 0.68097, 0.69127, 0.70212, 0.71352,
     8 0.72547, 0.73798, 0.75106, 0.76471, 0.77894, 0.79376, 0.80919,
     9 0.82523, 0.84190, 0.85922, 0.87720, 0.89586, 0.91521, 0.93528/
      data (v1    (i),i= 253, 255)/
     1 0.95609, 0.97765, 1.00000/
!
!         eps=1/30
      data (v2    (i),i=   1,  63)/
     1 1.00000, 0.88969, 0.82585, 0.77784, 0.74051, 0.71102, 0.68750,
     2 0.66865, 0.65354, 0.64147, 0.63191, 0.62448, 0.61886, 0.61481,
     3 0.61213, 0.61068, 0.61032, 0.61096, 0.61251, 0.61490, 0.61808,
     4 0.62199, 0.62660, 0.63188, 0.63779, 0.64433, 0.65147, 0.65919,
     5 0.66750, 0.67637, 0.68581, 0.69581, 0.70637, 0.71749, 0.72917,
     6 0.74143, 0.75425, 0.76766, 0.78166, 0.79625, 0.81145, 0.82727,
     7 0.84373, 0.86083, 0.87859, 0.89703, 0.91616, 0.93600, 0.95658,
     8 0.97790, 1.00000, 1.00000, 0.91038, 0.85257, 0.80759, 0.77181,
     9 0.74306, 0.71982, 0.70099, 0.68574, 0.67344, 0.66360, 0.65585/
      data (v2    (i),i=  64, 126)/
     1 0.64989, 0.64549, 0.64245, 0.64063, 0.63990, 0.64016, 0.64133,
     2 0.64332, 0.64609, 0.64959, 0.65377, 0.65860, 0.66406, 0.67011,
     3 0.67674, 0.68394, 0.69169, 0.69998, 0.70880, 0.71815, 0.72803,
     4 0.73843, 0.74936, 0.76081, 0.77278, 0.78529, 0.79834, 0.81193,
     5 0.82606, 0.84076, 0.85603, 0.87187, 0.88830, 0.90534, 0.92298,
     6 0.94126, 0.96018, 0.97975, 1.00000, 1.00000, 0.92195, 0.86822,
     7 0.82548, 0.79098, 0.76295, 0.74012, 0.72148, 0.70630, 0.69398,
     8 0.68407, 0.67621, 0.67012, 0.66556, 0.66236, 0.66036, 0.65945,
     9 0.65951, 0.66046, 0.66224, 0.66478, 0.66803, 0.67195, 0.67652/
      data (v2    (i),i= 127, 189)/
     1 0.68169, 0.68744, 0.69376, 0.70062, 0.70801, 0.71593, 0.72435,
     2 0.73328, 0.74271, 0.75263, 0.76305, 0.77396, 0.78537, 0.79727,
     3 0.80967, 0.82258, 0.83600, 0.84994, 0.86439, 0.87938, 0.89491,
     4 0.91099, 0.92763, 0.94483, 0.96262, 0.98101, 1.00000, 1.00000,
     5 0.92930, 0.87846, 0.83739, 0.80389, 0.77647, 0.75400, 0.73559,
     6 0.72052, 0.70825, 0.69834, 0.69045, 0.68431, 0.67969, 0.67640,
     7 0.67431, 0.67329, 0.67323, 0.67406, 0.67570, 0.67810, 0.68120,
     8 0.68496, 0.68934, 0.69433, 0.69988, 0.70598, 0.71261, 0.71976,
     9 0.72741, 0.73555, 0.74418, 0.75329, 0.76288, 0.77293, 0.78346/
      data (v2    (i),i= 190, 252)/
     1 0.79446, 0.80593, 0.81787, 0.83029, 0.84319, 0.85658, 0.87045,
     2 0.88482, 0.89970, 0.91508, 0.93099, 0.94742, 0.96440, 0.98192,
     3 1.00000, 1.00000, 0.93438, 0.88567, 0.84587, 0.81316, 0.78625,
     4 0.76409, 0.74588, 0.73093, 0.71873, 0.70885, 0.70097, 0.69481,
     5 0.69016, 0.68683, 0.68468, 0.68360, 0.68347, 0.68422, 0.68577,
     6 0.68807, 0.69106, 0.69471, 0.69897, 0.70382, 0.70923, 0.71517,
     7 0.72164, 0.72860, 0.73606, 0.74400, 0.75241, 0.76128, 0.77062,
     8 0.78040, 0.79065, 0.80134, 0.81248, 0.82408, 0.83613, 0.84864,
     9 0.86161, 0.87504, 0.88894, 0.90332, 0.91819, 0.93354, 0.94939/
      data (v2    (i),i= 253, 255)/
     1 0.96574, 0.98261, 1.00000/
!
!         eps=1/20
      data (v3    (i),i=   1,  63)/
     1 1.00000, 0.90643, 0.85194, 0.81000, 0.77662, 0.74968, 0.72779,
     2 0.70995, 0.69543, 0.68365, 0.67419, 0.66669, 0.66091, 0.65660,
     3 0.65362, 0.65180, 0.65104, 0.65124, 0.65231, 0.65420, 0.65684,
     4 0.66020, 0.66422, 0.66888, 0.67415, 0.68001, 0.68644, 0.69342,
     5 0.70094, 0.70899, 0.71756, 0.72664, 0.73624, 0.74635, 0.75697,
     6 0.76810, 0.77974, 0.79190, 0.80457, 0.81777, 0.83150, 0.84577,
     7 0.86058, 0.87595, 0.89189, 0.90840, 0.92549, 0.94319, 0.96150,
     8 0.98043, 1.00000, 1.00000, 0.92476, 0.87579, 0.83683, 0.80511,
     9 0.77906, 0.75760, 0.73991, 0.72535, 0.71341, 0.70372, 0.69594/
      data (v3    (i),i=  64, 126)/
     1 0.68983, 0.68519, 0.68185, 0.67966, 0.67851, 0.67832, 0.67899,
     2 0.68047, 0.68269, 0.68560, 0.68918, 0.69338, 0.69818, 0.70354,
     3 0.70945, 0.71590, 0.72286, 0.73033, 0.73830, 0.74675, 0.75569,
     4 0.76510, 0.77500, 0.78536, 0.79620, 0.80752, 0.81932, 0.83159,
     5 0.84435, 0.85760, 0.87134, 0.88559, 0.90034, 0.91561, 0.93140,
     6 0.94772, 0.96459, 0.98201, 1.00000, 1.00000, 0.93486, 0.88956,
     7 0.85273, 0.82229, 0.79702, 0.77603, 0.75859, 0.74413, 0.73222,
     8 0.72248, 0.71461, 0.70838, 0.70359, 0.70008, 0.69771, 0.69637,
     9 0.69596, 0.69641, 0.69766, 0.69964, 0.70230, 0.70561, 0.70953/
      data (v3    (i),i= 127, 189)/
     1 0.71404, 0.71910, 0.72469, 0.73080, 0.73741, 0.74451, 0.75209,
     2 0.76013, 0.76864, 0.77760, 0.78702, 0.79689, 0.80720, 0.81797,
     3 0.82918, 0.84084, 0.85295, 0.86552, 0.87855, 0.89204, 0.90600,
     4 0.92044, 0.93536, 0.95076, 0.96667, 0.98308, 1.00000, 1.00000,
     5 0.94123, 0.89849, 0.86320, 0.83375, 0.80911, 0.78851, 0.77131,
     6 0.75700, 0.74515, 0.73543, 0.72754, 0.72126, 0.71640, 0.71280,
     7 0.71033, 0.70887, 0.70834, 0.70866, 0.70976, 0.71159, 0.71409,
     8 0.71723, 0.72097, 0.72528, 0.73013, 0.73551, 0.74139, 0.74776,
     9 0.75461, 0.76191, 0.76967, 0.77788, 0.78653, 0.79561, 0.80512/
      data (v3    (i),i= 190, 252)/
     1 0.81506, 0.82543, 0.83623, 0.84746, 0.85911, 0.87119, 0.88371,
     2 0.89666, 0.91006, 0.92390, 0.93819, 0.95294, 0.96816, 0.98384,
     3 1.00000, 1.00000, 0.94560, 0.90474, 0.87062, 0.84192, 0.81778,
     4 0.79751, 0.78053, 0.76635, 0.75458, 0.74489, 0.73701, 0.73072,
     5 0.72582, 0.72217, 0.71964, 0.71811, 0.71751, 0.71774, 0.71874,
     6 0.72046, 0.72285, 0.72587, 0.72948, 0.73365, 0.73836, 0.74358,
     7 0.74930, 0.75549, 0.76215, 0.76926, 0.77682, 0.78480, 0.79321,
     8 0.80205, 0.81130, 0.82096, 0.83104, 0.84153, 0.85242, 0.86373,
     9 0.87545, 0.88759, 0.90014, 0.91311, 0.92650, 0.94033, 0.95458/
      data (v3    (i),i= 253, 255)/
     1 0.96928, 0.98441, 1.00000/
!
!      eps=1/10
!
      data (v4    (i),i=   1,  63)/
     1 1.00000, 0.93376, 0.89563, 0.86540, 0.84045, 0.81956, 0.80198,
     2 0.78718, 0.77474, 0.76433, 0.75570, 0.74863, 0.74294, 0.73848,
     3 0.73512, 0.73277, 0.73133, 0.73073, 0.73090, 0.73179, 0.73335,
     4 0.73554, 0.73832, 0.74166, 0.74554, 0.74994, 0.75482, 0.76019,
     5 0.76601, 0.77229, 0.77900, 0.78615, 0.79372, 0.80170, 0.81010,
     6 0.81891, 0.82812, 0.83773, 0.84775, 0.85817, 0.86900, 0.88023,
     7 0.89186, 0.90391, 0.91637, 0.92924, 0.94253, 0.95625, 0.97039,
     8 0.98498, 1.00000, 1.00000, 0.94807, 0.91439, 0.88679, 0.86349,
     9 0.84364, 0.82669, 0.81223, 0.79993, 0.78953, 0.78080, 0.77355/
      data (v4    (i),i=  64, 126)/
     1 0.76762, 0.76288, 0.75922, 0.75653, 0.75473, 0.75375, 0.75352,
     2 0.75400, 0.75514, 0.75689, 0.75922, 0.76211, 0.76551, 0.76941,
     3 0.77380, 0.77864, 0.78394, 0.78966, 0.79581, 0.80237, 0.80933,
     4 0.81669, 0.82444, 0.83257, 0.84109, 0.84998, 0.85925, 0.86890,
     5 0.87892, 0.88932, 0.90009, 0.91124, 0.92276, 0.93467, 0.94695,
     6 0.95963, 0.97269, 0.98615, 1.00000, 1.00000, 0.95574, 0.92491,
     7 0.89909, 0.87697, 0.85792, 0.84149, 0.82737, 0.81527, 0.80496,
     8 0.79625, 0.78896, 0.78296, 0.77810, 0.77430, 0.77144, 0.76946,
     9 0.76828, 0.76784, 0.76809, 0.76899, 0.77049, 0.77256, 0.77517/
      data (v4    (i),i= 127, 189)/
     1 0.77829, 0.78190, 0.78598, 0.79050, 0.79547, 0.80085, 0.80664,
     2 0.81283, 0.81941, 0.82637, 0.83371, 0.84142, 0.84949, 0.85792,
     3 0.86672, 0.87587, 0.88537, 0.89523, 0.90544, 0.91600, 0.92692,
     4 0.93820, 0.94984, 0.96183, 0.97419, 0.98691, 1.00000, 1.00000,
     5 0.96049, 0.93162, 0.90706, 0.88580, 0.86734, 0.85133, 0.83748,
     6 0.82556, 0.81536, 0.80670, 0.79943, 0.79339, 0.78849, 0.78460,
     7 0.78166, 0.77957, 0.77827, 0.77770, 0.77781, 0.77855, 0.77989,
     8 0.78179, 0.78422, 0.78715, 0.79057, 0.79444, 0.79876, 0.80350,
     9 0.80865, 0.81420, 0.82014, 0.82646, 0.83315, 0.84020, 0.84761/
      data (v4    (i),i= 190, 252)/
     1 0.85538, 0.86349, 0.87195, 0.88075, 0.88989, 0.89937, 0.90919,
     2 0.91935, 0.92985, 0.94068, 0.95186, 0.96338, 0.97524, 0.98745,
     3 1.00000, 1.00000, 0.96371, 0.93626, 0.91264, 0.89203, 0.87402,
     4 0.85833, 0.84471, 0.83295, 0.82285, 0.81424, 0.80698, 0.80095,
     5 0.79602, 0.79209, 0.78909, 0.78693, 0.78555, 0.78489, 0.78491,
     6 0.78555, 0.78678, 0.78856, 0.79086, 0.79366, 0.79694, 0.80067,
     7 0.80483, 0.80941, 0.81439, 0.81977, 0.82553, 0.83165, 0.83815,
     8 0.84499, 0.85219, 0.85972, 0.86760, 0.87582, 0.88436, 0.89324,
     9 0.90244, 0.91197, 0.92183, 0.93202, 0.94253, 0.95336, 0.96453/
      data (v4    (i),i= 253, 255)/
     1 0.97602, 0.98784, 1.00000/
!
!
!
       common /$muhad/ ieps, epsi
!
       e=min(ein, 100.)
       x=log10(e)
       call rndc(u)
!
       if(ieps .ge. 1 .and. ieps .le. 4)then
          call
     *    k4ptdi(vtbl(1,1,ieps), iu, ie, iu, 0., -2., .02, 1.,
     *    u, x, ans)
          v=ans*epsi**u
       elseif(ieps .eq. 5) then
          v=1.
       else
          write(*,*) ' ieps=',ieps,' undef for muhv1'
          stop
       endif
       end
!c         test muhadt
!*include  mucset
!c
!      integer ia(4)/100, 30, 20, 10/
!      character txt*70, capx*16
!      open(13, file='c2s5001.#gd.data')
!      capx='path(g/cm**2)'
!      call mucset(11.,22., .5, 5.5)
!      do 100 ieps=1, 4
!          epsi=1./ia(ieps)
!          call muhadr(epsi)
!          do 50 j=1, 1000
!              e=3.16e-2* 10.**( (j-1)/1.)
!              if(e .gt.  500.) goto 55
!              if(j .eq. 1) then
!                write(txt,'('' eps=1/'',i3,'' e='',f8.4)') ia(ieps), e
!              else
!                write(txt,'('' e='',f8.4)') e
!              endif
!              write(13) txt
!              write(13) capx
!              do 30 k=1, 20000
!                  call muhadt(e, t)
!                  write(13) t
!  30          continue
!              write(13) 1.e50
!  50      continue
!  55      continue
! 100  continue
!      end
      subroutine muhadt(e,  fpth)
!            sample the free path until the hadron production by muon
!         e: input. tev
!      fpth: output. random path length in g/cm**2
!
          call muhadx(e, tx)
          call rndc(u)
          if(tx .gt. 0.) then
             fpth=-log(u)/tx
          else
             fpth=1.e37
          endif
       end
!c         test muhadv
!*include  mucset
!c
!      integer ia(4)/100, 30, 20, 10/
!      character txt*70, capx*16
!      open(13, file='c2s5001.#gd.data')
!      capx='v '
!      call mucset(11.,22., .5, 5.5)
!      do 100 ieps=1, 4
!          epsi=1./ia(ieps)
!          call muhadr(epsi)
!          do 50 j=1, 1000
!              e=3.16e-2* 10.**( (j-1)/1.)
!              if(e .gt.  500.) goto 55
!              if(j .eq. 1) then
!                write(txt,'('' eps=1/'',i3,'' e='',f8.4)') ia(ieps), e
!              else
!                write(txt,'('' e='',f8.4)') e
!              endif
!              write(13) txt
!              write(13) capx
!              do 30 k=1, 20000
!                  call muhadv(e, v)
!                  write(13) v
!  30          continue
!              write(13) 1.e50
!  50      continue
!  55      continue
! 100  continue
!      end
       subroutine muhadv(e, v)
!        samples fractional energy loss of muons
!   e: input. tev
!   v: output. fractional energy loss
!
           call muhx1(e, tx1)
           call muhx2(tx2)
           call rndc(u)
           if(u .lt. tx2/(tx1+tx2) ) then
               call muhv2(v)
!              write(13) 2.
           else
               call muhv1(e, v)
!              write(13) 1.
           endif
       end
