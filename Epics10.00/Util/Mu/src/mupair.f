!         mupair:  init
!         mupart:  sample path length for pair creation (v>eps)
!         muparv:  sample fractional energy loss by pair (v>eps)
!         muparf:  average fractional energy loss by pair (v<eps)
!
      subroutine mupair(eps)
         common /$mupar/ ieps, epsi
         integer iepsa(5)/100, 30, 20, 10, 1/
         if(abs(eps-0.01)/eps .lt. 1.e-2) then
             ieps=1
         elseif(abs(eps-.0333333)/eps .lt. 1.e-2) then
             ieps=2
         elseif(abs(eps-.05)/eps .lt. 1.e-2) then
             ieps=3
         elseif(abs(eps-.1)/eps .lt. 1.e-2) then
             ieps=4
         elseif(   1.-eps   .lt. 1.e-2) then
             ieps=5
         else
             write(*,*) ' eps=',eps, ' invalid for mupair'
             stop
         endif
         epsi=1./iepsa(ieps)
         write(*,*)
     *   ' ----pair cre. by muon: taken into account for',
     *   ' v > 1/', iepsa(ieps), ' ---'
         write(*,*)
     *   '     lower energy pair cre. is treated as a part of ',
     *       'continuous energy loss.'
      end
!          test muparv
!       real*4 epsa(4)/0.01, 0.0333333, .05, .1/
!       character txt*70, cap*16
!       open(13,file='c2s5001.#gd.data')
!       cap='v'
!       do 200 ie=1, 4
!          eps=epsa(ie)
!          call mupair(eps)
!          write(txt,'('' eps='',f7.4)') eps
!          write(13) txt
!          write(13) cap
!          do 100 i=1, 100000
!              call muparv(v)
!              write(13) v
! 100      continue
!          write(13) 1.e50
! 200   continue
!     end
!        *********************************************************
!        *
!        * muparv: sample fractional energy loss of mu0n by
!        *        pair creation. (/(g/cm**2)
!        *
!        *********************************************************
!
      subroutine muparv(v)
         common /$mupar/ ieps, epsi
         if(ieps .eq. 1) then
            call mupre1(v)
         elseif(ieps .eq. 2) then
            call mupre2(v)
         elseif(ieps .eq. 3) then
            call mupre3(v)
         elseif(ieps .eq. 4) then
            call mupre4(v)
         elseif(ieps .eq. 5) then
            v=1.
         else
            write(*,*) ' err ieps=',ieps
         endif
      end
!            test muparb
!            (b-term)
!*include mucset
!         open(13, file='c2s5001.#gd.data')
!         call mucset(11.,22., .5, 5.5)
!         do 100 i=1, 70
!             e=0.03*10.**((i-1)/20.)
!             call muparb(e, b)
!             write(13) e, b/1.e-6
! 100     continue
!      end
      subroutine muparb(e, b)
!             e in tev.  dedtbe in  (g/cm**2)**-1
       include  'Zmucom.f'
!
          if(e .gt. .1) then
              bpe=.451* z2pba * (log(e/mu)-5.43)/(log(e/mu)-4.34)
          else
              bpe=max( z2pba*(.060*log(e/mu) -0.156), 0.)
          endif
          b=bpe*1.e-6
      end
!            test muparx
!*include mucset
!       real*4 epsa(4)/0.01, 0.0333333, .05, .1/
!       character txt*70, cap*16, capy*16
!       open(13,file='c2s5001.#gd.data')
!       cap='e(tev)'
!       capy='mfp(v>eps)(g/cm2)'
!       call mucset(11.,22., .5, 5.5)
!       do 200 ie=1, 4
!          eps=epsa(ie)
!          call mupair(eps)
!          write(txt,'('' eps='',f7.4)') eps
!          write(13) txt
!          write(13) cap, capy
!          do 100 i=1, 1000
!              e=.01*10.** ( (i-1)/20. )
!              if(e .gt. 1000.) goto 101
!              call muparx(e,  t)
!              write(13) e, 1./t
! 100      continue
! 101      continue
!          write(13) 1.e50, 1.e50
! 200   continue
!     end
       subroutine muparx(e, tx)
!          gives total x-section/(g/cm**2) for muon pair creation with
!          v>eps
           real*8 prmui
!
           common /$mupar/ ieps, epsi
           if(ieps .eq. 5)then
               tx=0.
           else
               call muparb(e, b)
               tx=b* ( prmui(1.d0)-prmui(dble(epsi)))
           endif
       end
       function prmui(v)
!c
!           indefinite integral of pair cre function of muon
         implicit real*8 (a-h,o-z)
          data dp/5.1d-3/
          prmui=  (1.d0+dp)*(1.d0/(v+dp) - log( (v+dp)/v )/dp)
      end
!            test mupart
!*include mucset
!       real*4 epsa(4)/0.01, 0.0333333, .05, .1/
!       character txt*70, cap*16
!       open(13,file='c2s5001.#gd.data')
!       cap='sampled path'
!       call mucset(11.,22., .5, 5.5)
!       do 200 ie=1, 4
!          eps=epsa(ie)
!          call mupair(eps)
!          do 100 i=1, 1000
!              e=.01*10.** ( (i-1) )
!              if(e .gt. 1000.) goto 101
!              write(txt,'('' eps='',f7.4, '' e='',f7.2)') eps, e
!              write(13) txt
!              write(13) cap
!              do 33 k=1, 20000
!                  call mupart(e,  t)
!                  write(13) t
!  33          continue
!              write(13) 1.e50
! 100      continue
! 101      continue
! 200   continue
!     end
      subroutine mupart(e, t)
!        samples free path for pair creation by muon.  energy loss
!        considered for fractional energy > epsi(given by mupair)
         common /$mupar/ ieps, epsi
!           t is given in unit of g/cm**2
!             standard rock.
!           eps       e     mfp(g/cm**2)
!           1/30      0.01    3821
!           1/30      .1      1754
!           1/30      1       1258
!           1/30      100     1122
!
!           1/100     1       166.5
!
!           1/20      1       2717
!
!           1/10      1       10300
!         total x-section
         call muparx(e, tx)
         if(tx .eq. 0.) then
            t=1.e37
         else
            call rndc(u)
            t=-log(u)/tx
!           write(13) u, t
         endif
       end
!     open(13, file='c2s5001.#gd.data')
!     do 100 i=1, 10000
!         call mupre1(v)
! 100 continue
!     end
      subroutine mupre1(v)
!        muon pair eps=1/100
!
!       xmin=   -5.189     xmax=   -3.483
      plsq1(x)=((-.1818992e-01*x-.3865440    )*x-2.546626    )*x
     *-4.370152
!
!       xmin=   -3.826     xmax=   -2.086
      plsq2(x)=((0.3488257e-01*x+0.4117552    )*x+ 1.151911    )*x+ 1.07
     *5700
!
!       xmin=   -2.553     xmax=  -0.9087
      plsq3(x)=((-.2146920e-01*x-.4787301e-01)*x-.8291966e-01)*x-.183397
     *9e-01
!
!       xmin=   -1.134     xmax=   0.3508e-01
      plsq4(x)=((-.6359987e-02*x+0.6492972e-02)*x-.1562950e-01)*x+0.9963
     *192e-02
!
!
!                            f(1.)-f(eps)      eps*(eps+dp)**2
          data  eps/.01/,cnst1/2858.454/, cnst2/2.2801e-6/
          data  dp/5.1e-3/
!
         call rndc(u)
         if(u .gt. .99) then
             alfa= cnst1*(1.-u)*cnst2
             v=eps+alfa
         elseif(u .lt. 2.e-5) then
             alfa=cnst1* u* (1.+dp)**2
             v=1.- alfa
         else
             x=log10(u)
             if(x .gt. -1.) then
                 v=plsq4(x)
             elseif(x .gt. -2.5) then
                 v=plsq3(x)
             elseif(x .gt. -3.5) then
                 v=plsq2(x)
             else
                 v=plsq1(x)
             endif
         endif
!        write(13) u, v
      end
!     open(13, file='c2s5001.#gd.data')
!     do 100 i=1, 100000
!         call mupre2(v)
! 100 continue
!     end
      subroutine mupre2(v)
!
!       muon pair eps=1/30
!
!       xmin=   -4.453     xmax=   -2.921
      plsq1(x)=((-.4670496e-01*x-.6458234    )*x-3.002569    )*x-3.70142
     *1
!
!       xmin=   -3.225     xmax=   -1.533
      plsq2(x)=((0.8478802e-01*x+0.6513133    )*x+ 1.230665    )*x+0.874
     *9610
!
!       xmin=   -1.896     xmax=   0.3221e-01
      plsq3(x)=(((0.3310572e-02*x-.1010634e-01)*x+0.2237669e-01)*x-.4283
     *034e-01)*x+0.3327289e-01
!
!                            f(1.)-f(eps)      eps*(eps+dp)**2
          data  eps/.033333333/,cnst1/371.28077/, cnst2/4.923737e-5/
          data  dp/5.1e-3/
!
         call rndc(u)
         if(u .gt. .98) then
             alfa= cnst1*(1.-u)*cnst2
             v=eps+alfa
         elseif(u .lt. 2.e-5) then
             alfa=cnst1* u* (1.+dp)**2
             v=1.- alfa
         else
             x=log10(u)
             if(x .gt. -1.8) then
                 v=plsq3(x)
             elseif(x .gt. -3.1) then
                 v=plsq2(x)
             else
                 v=plsq1(x)
             endif
         endif
!        write(13) u, v
      end
!     open(13, file='c2s5001.#gd.data')
!     do 100 i=1, 100000
!         call mupre3(v)
! 100 continue
!     end
      subroutine mupre3(v)
!
!           eps=1/20
!       xmin=   -4.973     xmax=   -2.530
      plsq1(x)=(((-.1935341e-02*x-.6902999e-01)*x-.6927310    )*x-2.7571
     *48    )*x-2.884690
!
!       xmin=   -2.944     xmax=   -1.439
      plsq2(x)=((0.1031246    *x+0.6872867    )*x+ 1.085623    )*x+0.715
     *8033
!
!       xmin=   -1.715     xmax=   0.4106e-01
      plsq3(x)=((-.2469236e-01*x+0.2412995e-01)*x-.6468838e-01)*x+0.4974
     *359e-01
!                            f(1.)-f(eps)      eps*(eps+dp)**2
          data  eps/.05/,cnst1/175.11765/, cnst2/1.518005e-4/
          data  dp/5.1e-3/
!
         call rndc(u)
         if(u .gt. .96) then
             alfa= cnst1*(1.-u)*cnst2
             v=eps+alfa
         elseif(u .lt. 1.0e-5) then
             alfa=cnst1* u* (1.+dp)**2
             v=1.- alfa
         else
             x=log10(u)
             if(x .gt. -1.50) then
                 v=plsq3(x)
             elseif(x .gt. -2.7) then
                 v=plsq2(x)
             else
                 v=plsq1(x)
             endif
         endif
!        write(13) u, v
      end
!     open(13, file='c2s5001.#gd.data')
!     do 100 i=1, 100000
!         call mupre4(v)
! 100 continue
!     end
      subroutine mupre4(v)
!
!
!           eps=1/10
!
!       xmin=   -3.705     xmax=   -1.873
      plsq1(x)=((-.4124709e-01*x-.4717042    )*x-1.813875    )*x-1.35316
     *9
!
!       xmin=   -2.116     xmax=  -0.6956
      plsq2(x)=((0.7702798e-01*x+0.3945251    )*x+0.2330459    )*x+0.223
     *8593
!
!       xmin=  -0.9640     xmax=  -0.8318e-02
      plsq3(x)=(0.1034741    *x-.1019843    )*x+0.1015701
!                            f(1.)-f(eps)      eps*(eps+dp)**2
          data  eps/.1/,cnst1/46.288379/, cnst2/1.104601e-3/
          data  dp/5.1e-3/
!
         call rndc(u)
         if(u .gt. .95) then
             alfa= cnst1*(1.-u)*cnst2
             v=eps+alfa
         elseif(u .lt. 2.0e-4) then
             alfa=cnst1* u* (1.+dp)**2
             v=1.- alfa
         else
             x=log10(u)
             if(x .gt. -0.90) then
                 v=plsq3(x)
             elseif(x .gt. -2.0) then
                 v=plsq2(x)
             else
                 v=plsq1(x)
             endif
         endif
!        write(13) u, v
      end
!          test muparf
!*include mucset
!       real*4 epsa(4)/0.01, 0.0333333, .05, .1/
!       character txt*70, cap*16, capy*16
!       open(13,file='c2s5001.#gd.data')
!       cap='e(tev)'
!       capy='b(v<eps)'
!       call mucset(11.,22., .5, 5.5)
!       do 200 ie=1, 4
!          eps=epsa(ie)
!          call mupair(eps)
!          write(txt,'('' eps='',f7.4)') eps
!          write(13) txt
!          write(13) cap, capy
!          do 100 i=1, 1000
!              e=.01*10.** ( (i-1)/20. )
!              if(e .gt. 1000.) goto 101
!              call muparf(e, be)
!              write(13) e, be/1.e-6
! 100      continue
! 101      continue
!          write(13) 1.e50, 1.e50
! 200   continue
!     end
      subroutine muparf(e, dedtbe)
!         fractional energy loss by v<eps (/(g/cm**2))
          common /$mupar/ ieps, epsi
          data  dp/5.1e-3/
          call muparb(e, b)
          dedtbe=b* (1.+dp)/(epsi+dp)*epsi
      end
