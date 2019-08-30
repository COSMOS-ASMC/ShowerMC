!//job job c2s5001,class=b
!//  exec fortclg,
!//   parm.fort='nos,num'
!//fort.sysinc dd dsn=c2g5100.cosmos.gem,disp=shr
!//fort.sysin dd *
!c **********************************************************
!c               test mufolw
!c    use mujcl0 to make library in #load.load
!c    at tss mode, use mylib lib1(#load.load)
!c                 alloc f(sysinc) da('c... cosmos.gem')
!c   for do loop by emu; use batch job
!c            (with #load.load and c2g5600.lib.load)
!c    ***********************************************************
!c    see $$$$   part where the depth of particle depth is obtained
!c    ***********************************************************
!c
!c----------------------------------------
!      parameter (n=100)
!      dimension sdep(n), eo(n), xo(n), yo(n), txo(n), tyo(n), tzo(n)
!c     data sdep/1000.e2, 2000.e2, 3000.e2, 4000.e2,
!c   *      5000.e2, 6000.e2, 7000.e2, 8000.e2/
!       character ttl*70, capx*16/'r(10**3 hg/cm**2)'/
!       open(13,file='c2s5001.#gd.data')
!        mode=3
!        epsp=0.1
!        epsb=0.05
!        epsh=0.05
!c          standard rock
!        z=11.
!        a=22.
!        zba=.5
!        z2ba=5.5
!        rho=2.8
!c         frejus rock
!        z=11.
!        a=22.
!        zba=.5
!        z2ba=5.87
!        rho=2.73
!c          kgf rock
!        z=12.8
!        a=25.8
!        zba=.493
!        z2ba=6.30
!        rho=3.02
!c
!        emin=.5e-3
!c
!c     do 10 i=1,n
!c        sdep(i)=1000.e2 + (i-1)*100.e2
!c 10  continue
!      sdep(1)=20000.e2
!c
!      call mucset(z, a, zba, z2ba)
!      call mupair(epsp)
!      call mubrem(epsb)
!      call muhadr(epsh)
!c
!      call muradl(z, zba, z2ba,   x0ing)
!      write(*,*) ' x0ing=',x0ing
!      call mufol0(rho, x0ing, emin,  mode)
!c     call mufol1(sdep, n)
!      call mufol1(sdep, 1)
!      do 200 ie=1, 12
!          emu= 100./10.**( (ie-1)/3.)
!          write(ttl,'('' emu='',f7.3,''tev'')') emu
!          write(13) ttl
!          write(13) capx
!          do 100 j=1, 2000
!             call mufolw(emu, 0., 0., 0.d0, 0.d0, 1.d0,
!    *        eo, xo, yo, txo, tyo, tzo, m)
!c            write(*,*) ' m=',m
!c            write(*,*) ' e=',(eo(i), i=1, m)
!c            write(*,*) ' x=',(xo(i), i=1, m)
!c            write(*,*) ' y=',(yo(i), i=1, m)
!c            write(*,*) 'tx=',(txo(i), i=1, m)
!c            write(*,*) 'ty=',(tyo(i), i=1, m)
!c
!c            write(*,*) ' -------------------------- '
!c            write(13) (1000.e2 + (m-1)*100.e2)/100.
! 100     continue
!         write(13) 1.e50
! 200  continue
!      end
!     ************************************************************
!     *
!     * mufolw: follow a muon in rock.
!     *     to follow 2000 muons of 10 tev until they die,
!     *     this  need 16 sec. (mode=3)
!     *                10 sec. (mode=1)
!     *
!     * 100 tev: mode=3 --> 22 sec
!     *
!     *   processes considered: energy loss by
!     *     ionization, knockon, brems, paircre, nuclear interactions.
!     *     only the muon is traced.  particles produced by
!     *     muon interactions are not traced.
!     *
!     ************************************************************
!       before calling subroutines here, the user must
!      call mucset(z, a, zba, z2ba)         ;to set consts
!      call mupair(epsp)                    ;to set min. e for pair
!      call mubrem(epsb)                    ;to set min. e for brem
!      call muhadr(epsh)                    ;to set min. e for n.i
!      call muradl(z, zba, z2ba,   x0ing)   ;to compute rad.l of rock
! see the above test prog. for further details.
!
!  /usage/ call mufol0(rho, x0ingi, emin, mode)    ; init
!          call mufol1(sdep, n)                    ; init
!          call mufolw(ein, xin, yin, txin, tyin, tzin,
!         *            eo, xo, yo, txo, tyo,
!         *            m)                          ; for each mu
!
!   rho:   input.  average rock density in g/cm**3
!   x0ingi:   input.  usual radiation length in g/cm**2
!   emin:  input.  minimum energy of particle to be followed (tev)
!   mode:   input.  1---> 1 dimensional
!                   3---> three dimensional
!                   4---> 4 dimensional
!
!   sdep(n):  input. slant depth along 1ry where detector is located
!                   (g/cm**2).
!      n:  input.  # of observation depth.
!
!   xin,yin:  input. muon position in 1ry system.
!   txin,tyin, tzin: input  real*8
!              muon direction cos in 1ry system.
!        ein:  input. muon energy in tev
!   xo(m),yo(m):  output. xo(i),yo(i) muon position at i-th  sdep
!             in 1ry system. (z-axis is on the detector)
!
!  txo(m),tyo(m): ouput.
!             txo(i), tyo(i) is the muon direction cos (1st and 2nd)
!             in 1ry system at i-th sdep.
!    eo(m):    output. eo(i) is the muon energy at i-th sdep.
!       m:     ouput. # of sdep where muon is observed.
       subroutine mufol0(rhoin, x0ingi, eminin, modein)
       include  'Zmufol.f'
            rho=rhoin
            x0ing=x0ingi
            emin=eminin
            x0=x0ing/rho
            mu=105.6e-6
            mode=modein
            if(mode .eq. 1 .or. mode .eq. 3 .or. mode .eq. 4) then
            else
                write(*,*) ' mode=',mode,' invalide'
                stop
            endif
       end
       subroutine mufol1(sdep, n)
       include  'Zmufol.f'
           dimension sdep(n)
           ndep=min(n, ndepmx)
           if(ndep .lt. n) then
               write(*,*) ' too large # of sdep '
               write(*,*) ' only first ',ndep,' are used'
           endif
!              convert sdepth into cm
            do   i=1, ndep
              sdepth(i)=sdep(i)/rho
!$$$$$$$$$$$$
!              write(*,*) ' sdepth=',sdepth(i)/100, ' m'
!$$$$$$$$$$$$
            enddo
       end
!
       subroutine
     * mufolw(ein, xin, yin, txin, tyin, tzin, eo,  xo, yo, txo, tyo,
     *   tzo,  m)
       include  'Zmufol.f'
        dimension xo(*), yo(*), txo(*), tyo(*), tzo(*), eo(*)
        real*8 txin, tyin, tzin
!
         x=xin
         y=yin
         z=0.d0
         e=ein
         wx=txin
         wy=tyin
         wz=tzin
         tm=0.d0
         m=0
!           depth counter
         idep=0
!
!         *** until loop*** 
         do while (.true.)
!                 fix process and free path dt(g/cm**2)
!                 and dl=dt*x0(cm)
!                 and compute tentative new pos. xp, yp, zp
              call munewp
!                 check depth if pass the detector depth
!                 if so, truncate the pass there.
              call mucdep(icon)
              if(icon .eq.  0) then
!                   not reach observation depth
!                   treate scattering, energy loss, interaction
                  call muint0(jcon)
              elseif(icon .ge. 1) then
!                      accross observation depth
!                      path is trancated at the observation point
                  dl=sqrt((xp-x)**2+ (yp-y)**2 + (zp-z)**2 )
                  dt=dl*rho
                  trunc=.true.
!                    treat scattering, energy loss
                  call muint0(jcon)
                  if(jcon .eq. 0) then
!                      observation
                      call muobsv(eo, xo, yo, txo, tyo, tzo)
                      m=idep
                  endif
              endif
         if         (icon .eq. 2 .or. jcon .ne. 0)
     *                      goto 100
         enddo
  100    continue
       end
       subroutine muobsv(eo, xo, yo, txo, tyo, tzo)
       include  'Zmufol.f'
         dimension eo(*), xo(*), yo(*), txo(*), tyo(*), tzo(*)
         eo(idep)=e
         xo(idep)=x
         yo(idep)=y
         txo(idep)=wx
         tyo(idep)=wy
         tzo(idep)=wz
         if(mode .eq. 4) then
             tim(idep)=tm
         endif
       end
       subroutine mutime(to)
       include  'Zmufol.f'
            parameter (cv=2.99e10/1.e-9)
!             cm lenght / cv  ---> nsec
         real*8 to(*)
!          to(i), i=1, m
         if(mode .eq. 4) then
              do   i=1, idep
                to(i)=tim(i)/cv
              enddo
         endif
       end
       subroutine munewp
!          fix the process and sample the path and compute
!          new tentative pos. xp,yp,zp
       include  'Zmufol.f'
!
!           save current position
          xbmv=x
          ybmv=y
          zbmv=z
          call muproc
!              trancate if dt is too long (set trunc=t/f)
          call mutrnc
          dl=dt/rho
          xp=x + dl*wx
          yp=y + dl*wy
          zp=z + dl*wz
       end
       subroutine mutrce
!            heading must be given by the user be for each
!            shower.
       include  'Zmufol.f'
!              trace informaion
           write(iotrc)  sngl(e), sngl(xbmv), sngl(ybmv),
     *                   sngl(zbmv), sngl(x),
     *                   sngl(y), sngl(z)
           return
       entry mutrcz
!              1 shwoer end
           write(iotrc)    -sngl(e), sngl(xbmv), sngl(ybmv),
     *                  sngl(zbmv), sngl(x),
     *                  sngl(y), sngl(z)
       end
       subroutine muproc
!
!      4    processes, i.e, pair creation, brems, nucrea int.
!            muon decay  are considered
!
!
       include  'Zmufol.f'
!
        parameter ( c=2.99e10,
     1   t0mu=2.197e-6,
     3   t0muc=t0mu*c)
           if(e .lt. enoint) then
               proc='decy'
               call mudcyl(sngl(e), mu, t0muc, dt)
               dt=dt*rho
           else
               call mudcyl(sngl(e), mu, t0muc, td)
               td=td*rho
               call mubrmt(sngl(e), tb)
               call mupart(sngl(e), tp)
               call muhadt(sngl(e), th)
               dt=min(td,tb,tp,th)
               if(dt .eq. tp) then
                   proc='pair'
               elseif(dt .eq. tb) then
                   proc='brem'
               elseif(dt .eq. th) then
                   proc='nuci'
               else
                   proc='decy'
               endif
           endif
       end
       subroutine mutrnc
       include  'Zmufol.f'
!          truncate path if it is too long
!           max length movable
          if(e .lt. enoint) then
              if(mode .ge. 3)then
                  tmax=max( tcoef*e,  tmin)
              else
                  tmax=1.e37
              endif
          elseif(e .gt. .50) then
              tmax=2.e4
          else
              tmax=4.e4*e
          endif
          if(dt .lt. tmax) then
               trunc=.false.
          else
               dt=tmax
               trunc=.true.
          endif
       end
       subroutine mudcyl(e, am, t0c,  t)
!         sample decay length in cm
         t0cm=t0c*e/am
         call rndc(u)
         t =-log(u)*t0cm
       end
       subroutine muint0(icon)
!            interaction routine
!             icon=0:  active ptcl exists
!                  1: no active ptcl exists
       include  'Zmufol.f'
!
!                scattering and energy loss
              call muscel(icon)
!                   update position
              call muudtp
              if(icon .eq. 0) then
                  if(trunc) then
                     icon=0
                  else
                     call muint1(icon)
                  endif
              endif
!$$$$$$$$$$$$$$$$
!             if(icon .ne. 0) then
!                write(13) sngl(z*rho/1.e5)
!             endif
!$$$$$$$$$$$$$$$$$
       end
       subroutine muscel(icon)
       include  'Zmufol.f'
!          scattering and energy loss
          ee1=e
!              compute energy loss rate
          call mudedx(sngl(e), dedt)
          if(e .ge. enoint) then
             call muparf(sngl(e), dedtp)
             call mubrmf(sngl(e), dedtb)
             call muhadf(sngl(e), dedth)
             dedt= (dedt + (dedtp+dedtb+dedth)*e)
          endif
          e=ee1-dedt*dt
          if(e .lt. emin) then
!                  death
              icon=1
!              get stop pos.
              dt=(ee1-emin)/dedt
              e=emin
              dl=dt/rho
              xp=x+dl*wx
              yp=y+dl*wy
              zp=z+dl*wz
              if(mode .ge. 3) then
                  call muscat(dt/x0ing)
              endif
          else
              icon=0
              if(mode .ge. 3) then
                  call muscat(dt/x0ing)
              endif
          endif
      end
      subroutine muudtp
       include  'Zmufol.f'
!               update position when new pos. is in the same
!               media
          real*8 beta
!
          if(mode .eq. 4) then
              if(e .gt. mu *35.)then
                  beta=1.d0
              elseif(e .gt. mu*10.) then
                  beta=1. - (mu/e)**2/2
              else
                  beta=sqrt(1. - (mu/e)**2)
              endif
              tm=tm+dl/beta
          endif
          x=xp
          y=yp
          z=zp
!                now x,y,z are updated
          if(trace) then
              call mutrce
          endif
       end
       subroutine muint1(icon)
       include  'Zmufol.f'
          if(proc .eq. 'pair') then
             call muparv( v)
          elseif(proc .eq. 'brem') then
             call mubrmv(sngl(e), v)
          elseif(proc .eq. 'nuci') then
             call muhadv(sngl(e), v)
          elseif(proc .eq. 'decy') then
             v=1.
          else
             write(*,*) ' proc=',proc, ' undef in muint1'
             stop
          endif
          e=(1.-v)*e
          if(e .lt. emin) then
             icon=1
          else
             icon=0
          endif
      end
!     ****************************************************************
!     *                                                              *
!     * muscat: cause charged ptcl scattering and compute new coord. *
!     *                                                              *
!     ****************************************************************
!
!  /usage/
!         call muscat(t)
!
!      t:  path of charged ptcl in r.l.
!
!
      subroutine muscat(t)
!
       include  'Zmufol.f'
           parameter (pi=3.141592)
           real*4 t
           real*8 w1, w2, w3
           real*8 wa, wb, wc
!             sample theta
           if(t .gt. 0.) then
               if(molier) then
                  call muang1(t,sngl(e*ee1), theta)
               else
                  call muang2(t, sngl(e*ee1), theta)
               endif
               if(theta .lt. 0.01) then
!                     cos
                   wc=1.-theta**2/2
                   sint=theta
               else
                   theta=min(theta, pi)
                   wc=cos(theta)
                   sint=sin(theta)
               endif
!
               call kcossn(cs, sn)
!
!                 wa,wb,wc: direction cos of scattering angle
!                             sample displacement correlated to theta
!
               wa=sint*cs
               wb=sint*sn
               tmp=t/2
               avx=tmp*wa
               avy=tmp*wb
!                   dispersion
               disp=sqrt(t/(6.*ee1*e) )*es*t / 2
!                     sample 2 independent gaussian variables
!                 with mean 0 and var 1
               call kgauss(0.,1.,g1,g2)
               dx=g1*disp+avx
               dy=g2*disp+avy
!
!                    displacement
!
               r=sqrt(dx*dx+dy*dy)
!                    direction cos of vector r in original sys.
               if(r .gt. 0.) then
                   w1=dx/r
                   w2=dy/r
                   w3=0.
!                       transform w1,w2,w3 to original sys.
                   call mutrns(2, wx, wy, wz, w1, w2, w3)
!                     r is in cm.
                   r=r*x0
!                    xp etc is pos.
!                    without scattering; add scattering effect.
!                    r*w1 etc is displacement by scattering
                   xp=r*w1+xp
                   yp=r*w2+yp
               endif
!
!                  convert scattering angle at end of path to
!                  original system
!                  wx,  etc are new direction cosine
                call mutrns(1, wx, wy, wz, wa, wb, wc)
            endif
      end
!     ****************************************************************
!     *                                                              *
!     * muang1: samples scattering angle by molier theory
!     *                                                              *
!     ****************************************************************
!
! /usage/
!        call muang1(t, e1*e2, theta)
!
!     t:  ptcl     path length in r.l.
!    e1*e2: path top and end energy (tev**2)
! theta:  sampled angle in radian
!
!
      subroutine muang1(t, e2in, teta)
!
         real*4 t, e2in, teta
!
!            used in moliere thory
          data esd/3.3e-5/
!          2 dimensional as ubscat(7,7),ubsca2(11,7)
          real ubscat( 49 ),ubsca2( 77 )
!
!      ubscat:  containes x for u=0.6 to 0.9 step 0.05 and for b=4.5 to
!               16.5 step 2.
!      ubsca2:  containes log(x) from u=0.9 to 1 step 0.01 in sqrt(1-
!        u)
!
!
      data (ubscat(i),i=   1,  49)/
     1 0.9370, 1.1119, 1.3331, 1.6267, 2.0406, 2.6703, 3.7352, 0.9313,
     2 1.0890, 1.2809, 1.5239, 1.8492, 2.3236, 3.1242, 0.9285, 1.0790,
     3 1.2593, 1.4832, 1.7753, 2.1877, 2.8587, 0.9268, 1.0734, 1.2476,
     4 1.4614, 1.7365, 2.1173, 2.7196, 0.9257, 1.0698, 1.2402, 1.4479,
     5 1.7126, 2.0745, 2.6355, 0.9249, 1.0674, 1.2351, 1.4387, 1.6965,
     6 2.0457, 2.5796, 0.9243, 1.0656, 1.2315, 1.4320, 1.6849, 2.0251,
     7 2.5398/
!
!
      data (ubsca2(i),i=   1,  72)/
     1 1.3178, 1.4751, 1.6416, 1.8340, 2.0663, 2.3518, 2.6934, 3.1142,
     2 3.6467, 4.2571, 4.6052, 1.1392, 1.2850, 1.4435, 1.6184, 1.8246,
     3 2.0805, 2.4108, 2.8389, 3.4102, 4.1338, 4.6052, 1.0504, 1.1821,
     4 1.3262, 1.4879, 1.6753, 1.9080, 2.2150, 2.6362, 3.2238, 4.0232,
     5 4.6052, 1.0005, 1.1225, 1.2553, 1.4042, 1.5766, 1.7895, 2.0734,
     6 2.4788, 3.0718, 3.9256, 4.6052, 0.9691, 1.0845, 1.2094, 1.3483,
     7 1.5087, 1.7049, 1.9674, 2.3529, 2.9440, 3.8376, 4.6052, 0.9476,
     8 1.0585, 1.1777, 1.3092, 1.4601, 1.6429, 1.8863, 2.2503, 2.8342,
     9 3.7572, 4.6052, 0.9321, 1.0397, 1.1546, 1.2806, 1.4241, 1.5963/
!
          data (ubsca2(i),i= 73, 77)/
     1     1.8234, 2.1658, 2.7384, 3.6835, 4.6052/
!
!                 sqrt(0.1)/10
!
           data du2/3.162278e-2/
!               to gev**2
          e2=e2in*1.e6
          call rndc(u)
!               by moliere theory
          sb=log(t) +10.33
          if(sb .lt. 3.) then
!                 t < 1.e-3
!                 gaussian approx.
              call muang2(t, e2in, teta)
          else
              b=(-0.01451*sb+1.32)*sb+0.7
              if(u .lt. 0.6) then
                   xtmp=((1.563*u-0.0625)*u+1.025)*u
              elseif(u .lt. .9) then
!                     2 dim. interpolation
                  call
     *            k4ptdi(ubscat, 7, 7, 7, 0.6, 4.5, 0.05, 2., u, b,
     *            xtmp)
              else
                  usq=0.3162278-sqrt(1.-u)
                  call
     *            k4ptdi(ubsca2,11, 7,11, 0.,  4.5, du2,  2., usq,b,
     *            xtmp)
                  xtmp=exp(xtmp)
              endif
              teta=sqrt(xtmp*esd/e2      *t*b)
          endif
       end
       subroutine muang2(t, e2, teta)
!                gaussian  approx. e2 tev**2
       include  'Zmufol.f'
           real*4 t, e2, teta
           call rndc(u)
           teta=es*sqrt(-log(u)*t/e2)
       end
      subroutine  mutrns(j, wx, wy, wz, w1, w2, w3)
!
!      angle transformation:
!         let suppose direction cosines wx,wy,wz are given in
!        a certain system (a-system=1ry system). w1,w2,w3 are direction
!        cosines in a system (b-system) whose z axis is given by
!        wx,wy,wz. (x- and y- axies of the b-system
!        are rather arbitrary).  this routine
!        transforms w1,w2, w3 to those in the a-system.
!
!    j: input.  1---> transformed values are returned in wx,
!                     wy, wz.
!               2---> transformed values are returned in
!                     w1, w2, w3.
!      wx,... w3 are all real*8
!
!
          implicit real*8 (a-h,o-z)
!
          data epsx/1.d-4/
!
          el2=wx**2
          em2=wy**2
          eps=el2+em2
          d=1.+wz
          if(abs(d) .gt. epsx) then
              a=el2/d - 1.
              b=wx*wy/d
              c=em2/d - 1.
              tmpa=a*w1 + b*w2 + wx*w3
              tmpb=b*w1 + c*w2 + wy*w3
          else
              tmpa=w2
              tmpb=w1
          endif
          tmpc=wx*w1 + wy*w2 + wz*w3
!              check result
          eps=tmpa**2 + tmpb**2 + tmpc**2 - 1.d0
          if(abs(eps) .gt. epsx ) then
!                 renormalize
               anrm=sqrt(1.d0+eps)
               tmpa=tmpa/anrm
               tmpb=tmpb/anrm
               tmpc=tmpc/anrm
          endif
          if(j .eq. 1) then
              wx=tmpa
              wy=tmpb
              wz=tmpc
          elseif(j .eq. 2) then
              w1=tmpa
              w2=tmpb
              w3=tmpc
          else
              write(*,*) ' error of j=',j,' to mutrns'
              stop
          endif
       end
       subroutine mucdep( icon)
       include  'Zmufol.f'
          if(zp .ge. sdepth(idep+1)) then
              idep=idep+1
              tang=(sdepth(idep)-zbmv)/(zp-zbmv)
              xp= tang * (xp-xbmv) + xbmv
              yp= tang * (yp-ybmv) + ybmv
              zp= sdepth(idep)
              if(idep .eq. ndep)then
                  icon=2
              else
                  icon=1
              endif
          else
              icon=0
          endif
      end
      block data
       include  'Zmufol.f'
         data es/20.e-6/, enoint/30.e-3/, molier/.true./ ,
     *   tmin/30./,  tcoef/1.4e5/,  mu/105.6e-6/,
     *   iotrc/11/, trace/.false./
      end
!//*
!//lked.syslib dd dsn=c2s5001.#load.load,disp=shr
!//  dd
!//  dd
!//  dd
!//  dd dsn=c2g5600.lib.load,disp=shr
!//go.ft13f001 dd dsn=c2s5001.#gd.data,disp=shr
!//
       subroutine dummy
       end
