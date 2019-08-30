!             test brmmd
!*include mucset
!        implicit real*8 (a-h, o-z)
!        character txt*70, capx*16, capy*16
!        open(13,file='c2s5001.#gd.data')
!        call mucset(11., 22.,  .5, 5.5)
!        capx='v'
!        capy='v*brem(v)/d'
!        do 200 i=1, 4
!           e=.03162*10.**( (i-1)/1.)
!           write(txt,'('' e='',f8.2)') e
!           write(13) txt
!           write(13) capx, capy
!           do 100 j=1, 1000
!               v=1.d-5*10.**( (j-1)/40.)
!               if(v .ge. 1.) goto 101
!               write(13) sngl(v), sngl( v*brmmd(e, v) )
! 100       continue
! 101       continue
!           write(13) 1.e50, 1.e50
! 200    continue
!        end
!c---------------------------------------------------------------
!c          test <v> for exact value.
!*include  mucset
!         implicit real*8 (a-h, o-z)
!         real*4 z, a, zba, z2ba, rho
!         eps=.01d0
!         e=.5
!         z=11.
!         a=22.
!         zba=.5
!         z2ba=5.87
!         rho=2.73
!         call mucset(z, a, zba, z2ba)
!            xs=  brmmi(e, eps, 1.d0)
!            av=(vbrmmi(e, 1.d0)-vbrmmi(e, eps))/xs
!            write(*,*) ' <v>=', av
!         end
         function brmmd(e,v)
!              differential muon brem x-section.  / d
!              (rosental with correction)
!           (d=3.10 z'**2/ax10**-8 /(g/cm**2) )
            implicit real*8 (a-h, o-z)
       include  'Zmucom.f'
!
            data ak/191./, emass/.511e-6/
!                2/3 not 3/2 in rosental
            akm=2*ak*mu/emass/3.
            akm2=ak *sqrt(2.7182)/2 * mu/emass
            brmmd= (4.*(1.-v)/3.d0+v**2)/v *
     *            log(  akm/z3/z3 / (akm2*(mu/e)*v/(1.-v)/z3+1.))
         end
!             test brmmi
!*include  mucset
!        implicit real*8 (a-h, o-z)
!        character txt*70, capx*16, capy*16
!        open(13,file='c2s5001.#gd.data')
!        call mucset(11., 22., .5, 5.5)
!        capx='vmin'
!        capy='integral brem'
!        do 200 i=1, 7
!           e=.01*10.**( (i-1)/2.)
!           write(txt,'('' e='',f8.2)') e
!           write(13) txt
!           write(13) capx, capy
!           do 100 j=1, 1000
!               v=1.d-5*10.**( (j-1)/5.)
!               if(v .ge. 1.) goto 101
!               write(13) sngl(v), sngl( brmmi(e, v, 1.d0-1.d-7))
! 100       continue
! 101       continue
!           write(13) 1.e50, 1.e50
! 200    continue
!        end
      function brmmi(e, v1, v2)
!            integral of brmmd
         implicit real*8 (a-h, o-z)
         external brmk
         common /$/ ec
         data  epsa/1.d-8/
         ec=e

         call kdexpIntF(brmk, v1, v2, epsa, s, err, icon)
         if(icon .ne. 0) then
             write(*,*) ' icon=',icon, ' err=',err
         endif
         brmmi=s
      end
      function brmk(x)
         implicit real*8 (a-h, o-z)
         dimension x(2)
         common /$/ ec
         brmk=brmmd(ec, x(1))
      end
!             test vbrmmi
!*include  mucset
!        implicit real*8 (a-h, o-z)
!        character txt*70, capx*16, capy*16
!        real*8 va(5)/1.d-2, 0.033333333333d0, 0.05d0, 0.1d0, .99999d0/
!
!        open(13,file='c2s5001.#gd.data')
!        call mucset(11.,22., .5, 5.5)
!        capx='emu(tev)'
!        capy='f.e.l/d'
!        do 100 j=1, 5
!           v=va(j)
!           write(txt,'('' v='',f8.2)') v
!           write(13) txt
!           write(13) capx, capy
!           do 200 i=1, 1000
!               e=.001*10.**( (i-1)/40.)
!               if(e .gt. 1000.) goto 201
!               write(13) sngl(e), sngl(vbrmmi(e, v))
! 200       continue
! 201       continue
!           write(13) 1.e50, 1.e50
! 100    continue
!        end
!             vbremi for total x-section
!*include  mucset
!           implicit real*8 (a-h, o-z)
!           character txt*70, capx*16, capy*16
!          real*8 va(4)/1.d-2, 0.033333333333d0, 0.05d0, 0.1d0/
!*include $mucom
!c
!        open(13,file='c2s5001.#gd.data')
!        call mucset(11., 22., .5, 5.5)
!        capx='emu(tev)'
!        capy='xsec(v>eps)/d'
!        do 100 j=1, 4
!           v=va(j)
!           write(txt,'('' eps='',f8.5)') v
!           write(13) txt
!           write(13) capx, capy
!           do 200 i=1, 1000
!               e=.01*10.**( (i-1)/40.)
!               if(e .gt. 1000.) goto 201
!               write(13) sngl(e), sngl(brmmi(e, v, 1.d0-1.d-5))
! 200       continue
! 201       continue
!           write(13) 1.e50, 1.e50
! 100    continue
!        end
!c         to make sampling talbe for (u,v) for mu brem.
!*include mucset
!        implicit real*8 (a-h,o-z)
!        character txt*70, capx*16, capy*16
!        parameter (ie=15,iu=41, iv=1000)
!        dimension  vtbl(iu, ie), vu(iv), vw(4*iv),  ua(iv)
!        eps=1./100.d0
!        call mucset(11., 22., .5, 5.5)
!        open(13,file='c2s5001.#gd.data')
!        open(14,file='c2s5001.#gd2.data')
!        open( 7,file='c2s5001.#h.fort(v1)')
!        write(*,*)' eps=',eps
!        capx='u'
!        capy='v/eps**u'
!     do 300 j=1, 50
!        e=.01*10.**( (j-1)/3. )
!        write(*,*) ' e=',e
!        if(e .gt.  50.) goto 301
!        write(txt,'('' eps='',f8.5,'' e='',f7.2)') eps, e
!        write(13) txt
!        write(13) capx, capy
!        write(14) txt
!        write(14) capx, capy
!        v=1.
!        nv=0
!        do 100 i=1,  1000
!          v=v/10.d0**0.005
!          if(v .le. eps ) goto 101
!          u= brmmi(e, v, 1.d0-1.d-5)/ brmmi(e, eps,1.d0-1.d-5)
!          write(*,*) ' u=',u, ' v=',v
!          vn=eps**u
!          write(14) sngl(u), sngl(v/vn)
!          if(u .gt. 0. .and. u .lt. 1.)then
!              nv=nv+1
!              vu(nv)=v/vn
!              ua(nv)=u
!          endif
! 100   continue
! 101   continue
!c        make u->v table
!       write(13) 0.e0, 1.e0
!       iuu=2
!       do 120 u=.025d0, .975001d0, .025d0
!          m=5
!          epsx=1.d-4
!          call daklag(ua, vu, nv, u, m, epsx, f, vw, icon)
!          if(icon .ne. 0)then
!              write(*,*) ' icon=',icon
!          endif
!          vtbl(iuu, j) =f
!          iuu=iuu+1
!          write(13) sngl(u), sngl(f)
! 120   continue
!       write(13) 1.e0, 1.e0
!       write(13) 1.e50, 1.e50
!       vtbl(1, j)=1.
!       vtbl(41, j)=1.
!       write(14) 1.e50, 1.e50
! 300   continue
! 301   continue
!       call mkdt('v1  ',vtbl, 2, 41*(j-1),'f8.5,  ', 8,  1)
!     end
      function vbrmmi(e, v)
         implicit real*8 (a-h, o-z)
         external vbrmk
         common /$/ ec
         data  epsa/1.d-8/
         ec=e

         call kdexpIntF(vbrmk, 0.d0, v, epsa, s, err, icon)
         if(icon .ne. 0) then
             write(*,*) ' icon=',icon, ' err=',err
         endif
         vbrmmi=s
      end
      function vbrmk(x)
         implicit real*8 (a-h, o-z)
         dimension x(2)
         common /$/ ec
         vbrmk=brmmd(ec, x(1))*x(1)
      end
