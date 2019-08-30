!c       test murang
!       character*70 titl
!       character*16 capx, capy
!       open(13,file='c2s5001.#gd.data')
!       emu=100.
!       capx='rang(1000 hg)'
!       capy='prob'
!       do  10 i=1, 20
!          write(titl,'(''e='',f6.2)') emu
!          write(13)titl
!          write(13) capx, capy
!          do 100 r=0., 20., .1
!             call murang(log10(emu), r, f)
!             write(13) r, f
! 100      continue
!          emu=emu/10.**.1
!           write(13) 1.e50, 1.e50
!  10   continue
!       end
!      frg rock: range distribtution fucntion.
      subroutine murang(emul, r, f)
!      emul : input. log10 of energy of muon in tev
!      r    : input. range of mu in 1000 hg
!      f    : output. probability of having range r  (integral fdr from
!         0 to  inf)=1.
!    -- ----------------- formula valid at e>4. tev---------
!         x=log10(e/tev)
!       a*r**n * exp( - b*r**m)dr where r is range express in
!       1000 hg/cm**2)
!            log10(a)
!
!
!       xmin=   0.2828     xmax=    2.717
      plsqa(x)=((-.6522435e-01*x+0.3872343    )*x-1.673499    )*x-1.0411
     *14
!
!
!                n
!       xmin=   0.2315     xmax=    2.767
       plsqn(x)=0.4242518    *x+ 2.001287
!
!          log10(b)
!
!       xmin=   0.2828     xmax=    1.717
      plsqb1(x)=((-.8638042    *x+ 1.847070    )*x+ 1.503570    )*x
     *-8.435040
!
!       xmin=    1.283     xmax=    2.717
      plsqb2(x)=((0.2062451    *x-2.806771    )*x+ 8.216535    )*x
     *-11.64697
!
!              m
!
!       xmin=   0.2828     xmax=    2.717
      plsqm(x)=((-.4801480    *x+ 4.493530    )*x-12.30122    )*x+ 15.12
     *132
!
        if(r .le. 0.)then
            f=0.
        else
            if(emul .gt. 0.176) then
                x=emul
                if(emul .lt. .6020 .and. r .gt. rngmax(x))  then
                    f=0.
                else
                    a=10.**(plsqa(x))
                    en=plsqn(x)
                    if(x .lt. 1.333) then
                       b=10.**(plsqb1(x))
                    else
                       b=10.**(plsqb2(x))
                    endif
                    em=plsqm(x)
!                   f=a*r**en * exp(-b*r**em)
                    f=a* exp( -b*r**em + en*log(r))
                endif
            else
!                low energy. range can be approximated by average
                write(*,*) ' use delta function '
                f=0.
            endif
        endif
        end
        function rngmax(x)
!          max range  for given log10(emu/tev)
           rngmax=10.** ( (-7.149e-2*x + 0.7407)*x + 0.438 )
        end
        function rmxtoe( x)
!          inverse of  x=log10(rngmax).  rmxtoe in tev
          rmxtoe=10.**(  (.1221*x+1.24)*x -.5584 )
        end
        function avrang(emulog)
!        emulog: input. log10(emu/tev).
!        avrang: output. average range in 1000 h.g/cm**2
!
!          average range vs energy frjus.
!           elog---><>range
!       xmin=   -1.743     xmax=    2.742
      plsq1(x)=((0.1360355e-01*x-.1447803    )*x+0.5489631    )*x+0.3441
     *143
!
        avrang=10.**(plsq1(emulog))
      end
      function avrngp(e)
      plsq (x)= (0.1360355e-01*3*x-.1447803*2 )*x+0.5489631
!           derivative of avrang by e
        x=log10(e)
        avrngp=1./e * avrang(x) * plsq(x)
      end
!             <> range--->emu (tev)
      function avrtoe(avrl)
!        avrl: input. log10 of average range (/1000 h g/cm**2)
! avrtoe: output. corresponding energy (tev)
!       xmin=   -1.059     xmax=   0.4084
      plsq1(x)=(0.2084617    *x+ 1.337060    )*x-.5196140
!
!       xmin=   0.9690e-01 xmax=    1.100
      plsq2(x)=(( 3.910639    *x-3.671625    )*x+ 2.708297    )*x-.66777
     *75
!
        if(avrl .lt..3) then
            avrtoe=10.**(plsq1(avrl))
        else
            avrtoe=10.**(plsq2(avrl))
        endif
      end
