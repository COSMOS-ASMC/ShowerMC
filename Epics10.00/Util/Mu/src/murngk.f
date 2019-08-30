!        test murang
! !!!!!!!! note:  avrang etc here is not for kgf but for frej. they can
!                 be used for mudepi. (depth - intensity)
!       character*70 titl
!       character*16 capx, capy
!       open(13,file='c2s5001.#gd.data')
!       emu=464.1577
!       capx='rang(1000 hg)'
!       capy='prob'
!       do  10 i=1, 10
!          write(titl,'(''e='',f6.2)') emu
!          write(13)titl
!          write(13) capx, capy
!          do 100 r=0., 20., .1
!             call murang(log10(emu), r, f)
!             write(13) r, f
! 100      continue
!          emu=emu/10.**.3333333
!           write(13) 1.e50, 1.e50
!  10   continue
!       end
      subroutine murang(emu, r, f)
!      emu  : input. log10(energy of muon in tev)
!      r    : input. range of mu in 1000 hg
!      f    : output. probability of having range r  (integral fdr from
!         0 to  inf)=1.
!      kgf rock: range distribtution fucntion coeff.
!       a*r**n * exp( - b*r**m)dr where r is range express in
!       1000 hg/cm**2)
!    -------------------- formula valid at e>4. tev---------
!         x=log10(e/tev)
!            log10(a)
!
!
!       xmin=   0.2828     xmax=    2.717
      plsqa(x)=((-.1188245    *x+0.6235009    )*x-1.941816    )*x-.92198
     *18
!                n
!       xmin=   0.2315     xmax=    2.767
       plsqn(x)=0.4242518    *x+ 2.001287
!
!          log10(b)
!
!       xmin=   0.2828     xmax=    1.717
      plsqb1(x)=((-1.264644    *x+ 2.772780    )*x+ 1.547551    )*x
     *-8.980374
!
!       xmin=    1.616     xmax=    2.717
      plsqb2(x)=(-1.278649    *x+ 4.513782    )*x-8.513817
!              m
!       xmin=   0.2828     xmax=    2.717
      plsqm(x)=((-.8051989    *x+ 6.399025    )*x-15.89226    )*x+ 17.19
     *075
!
!
        if(r .le. 0.)then
            f=0.
        else
            x=emu
            if(x .gt. 0.176) then
                a=10.**(plsqa(x))
                en=plsqn(x)
                if(x .lt. 1.666) then
                   b=10.**(plsqb1(x))
                else
                   b=10.**(plsqb2(x))
                endif
                em=plsqm(x)
!               f=a*r**en * exp(-b*r**em)
                f=a* exp( -b*r**em + en*log(r))
            else
!                low energy. range can be approximated by average
                write(*,*) ' use delta function '
                f=0.
            endif
        endif
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
