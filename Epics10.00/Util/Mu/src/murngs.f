!        test murngs
!c      character*70 titl
!c      character*16 capx, capy
!c      open(13,file='c2s5001.#gd.data')
!c      emu=464.1577
!c      capx='rang(1000 hg)'
!c      capy='prob'
!c      do  10 i=1, 14
!c         write(titl,'(''e='',f6.2)') emu
!c         write(13)titl
!c         write(13) capx, capy
!c         do 100 r=0., 20., .1
!c            call murngs(emu, r, f)
!c            write(13) r, f
!c100      continue
!c         emu=emu/10.**.3333333
!c          write(13) 1.e50, 1.e50
!c 10   continue
!c      end
!      standard rock: range distribtution fucntion coeff.
!       a*r**n * exp( - b*r**m)dr where r is range express in
!       1000 hg/cm**2)
      subroutine murngs(emu, r, f)
!      emu  : input. energy of muon in tev
!      r    : input. range of mu in 1000 hg
!      f    : output. probability of having range r  (integral fdr from
!         0 to  inf)=1.
!    -------------------- formula valid at e>4. tev---------
!         x=log10(e/tev)
!            log10(a)
!       xmin=   0.2315     xmax=    2.767
        plsqa(x)=(0.9661061e-01*x-1.322639    )*x-1.175355
!                n
!       xmin=   0.2315     xmax=    2.767
        plsqn(x)=0.4242518    *x+ 2.001287
!          log10(b)
!       xmin=   0.2315     xmax=    1.433
        plsqb1(x)= 3.685186    *x-9.420974
!          log10(b)
!       xmin=    1.232     xmax=    2.767
        plsqb2(x)=(0.1709792    *x-1.417317    )*x-2.953930
!              m
!       xmin=   0.2315     xmax=    1.433
        plsqm1(x)=-6.390082    *x+ 13.15197
!             m
!       xmin=    1.232     xmax=    2.767
        plsqm2(x)=0.3898764e-01*x+ 4.794041
!
!
        if(r .le. 0.)then
            f=0.
        else
            x=log10(emu)
            a=10.**(plsqa(x))
            en=plsqn(x)
            if(x .lt. 1.3333) then
               b=10.**(plsqb1(x))
               em=plsqm1(x)
            else
               b=10.**(plsqb2(x))
               em=plsqm2(x)
            endif
!           f=a*r**en * exp(-b*r**em)
            f=a*  exp( -b*r**em + en*log(r))
        endif
        end
