!             test  vmhd1   or vmhd2
!*include mucset
!        implicit real*8 (a-h, o-z)
!       character*70 txt
!       character*16 capx,capy
!       open(13, file='c2s5001.#gd.data')
!       capx=' v'
!       capy='v*fai1'
!       call mucset(11.,22., .5, 5.5)
!       do 100 j=1, 10
!          e=.01*10.**( (j-1)/2.)
!          if(j .eq. 1) then
!              write(txt,'('' v*fai(inf) of muon nuclear interaction'',
!    *         '' e='',f8.2)') e
!          else
!              write(txt,'('' e='',f8.2)') e
!          endif
!          write(13) txt
!          write(13) capx, capy
!          do 80 k=1, 10000
!               v=1.e-4*10.**( (k-1)/40.)
!               if(v .gt. 0.999) goto 81
!               tmp=vmhd1(e,v)
!c              tmp=vmhd2(v)
!               write(13) sngl(v), sngl(tmp)
!  80      continue
!  81      continue
!          write(13) 1.e50, 1.e50
! 100   continue
!       end
        function vmhd1(e, v)
!             v*fai1(e,v) for nuclear interaction  of mu
         implicit real*8 (a-h, o-z)
       include  'Zmucom.f'
         vmhd1= ((1.d0+(1.d0-v)**2)*log(2*(e/mu)*8.88 * (1.-v)/v)
     *       + 2*(v-1.)
!    *       + ( (mu/e)**2*v**2/(1.-v) - 2*8.88*mu/e*v)/2
!    *       + v*mu/e/8.88
     *                )/2
        end
        function vmhd2(v)
!             v*fai2(e,v) for nuclear interaction  of mu
         implicit real*8 (a-h, o-z)
       include  'Zmucom.f'
         parameter (alm2=0.365e-6)
         vmhd2=   (v-1.)
!    *    + v*mu/e/8.88/2
!    *           +alm2/e**2/4*log(  ( 2*v+alm2/mu /8.88/e)/
!    *          (v**2/(1.-v)*mu/e/8.88+ alm2/mu/8.88/e) )
     *       +   (1.-v+ v**2/2*(1.+2* mu**2/alm2) ) *
     *           log( (1.+ (1.-v)/v**2 * alm2/mu**2) )
!    *  log( ( v+ (1.-v)/v*alm2/mu**2 )/(v + alm2/2/.938/e) )
        end
!*include mucset
!        implicit real*8 (a-h, o-z)
!       character*70 txt
!       character*16 capx,capy
!       real*4  epsa(5)/0.01, .03333333, .05, .1, 0.999999/
!       open(13, file='c2s5001.#gd.data')
!       capx=' e'
!       capy='integ(fai1)'
!       call mucset(11.,22., .5, 5.5)
!       do 80 k=1, 5
!           eps=epsa(k)
!           if(k .eq. 1) then
!              write(txt,
!    *          '(''integral 0-eps of v*fai(inf) of had. prod. by mu'',
!    *           '' eps='',f7.4)') eps
!           else
!              write(txt,'('' eps='',f7.4)') eps
!           endif
!           write(13) txt
!           write(13) capx, capy
!           do 100 j=1, 1000
!               e=.01*10.**( (j-1)/20.)
!               if(e .gt. 1000.) goto 101
!                 write(13) sngl(e), sngl(vmhi1(e, 1.e-5, eps))
! 100       continue
! 101       continue
!           write(13) 1.e50, 1.e50
!  80   continue
!       end
        function vmhi1(e, v1, v2)
!           integral of vmhd1 from v1 to v2
         implicit real*8 (a-h, o-z)
         external vmhd1k
         common /$/ ec
         data epsa/1.d-8/
         ec=e
         call kdexpIntF(vmhd1k, v1, v2, epsa, s, err, icon)
         if(icon .ne. 0) then
             write(*,*) ' icon=',icon, ' err=',err
         endif
         vmhi1=s
      end
      function vmhd1k(x)
         implicit real*8 (a-h, o-z)
         common /$/ ec
         vmhd1k=vmhd1(ec, x)
      end
!            test vmhi2
!*include mucset
!        implicit real*8 (a-h, o-z)
!       real*4  epsa(5)/0.01, .03333333, .05, .1, 0.999999/
!       call mucset(11.,22., .5, 5.5)
!       do 80 k=1, 5
!           eps=epsa(k)
!           write(*,*) ' eps=',eps, ' int(0-eps)=', vmhi2(1.e-5, eps)
!  80   continue
!       end
        function vmhi2(v1, v2)
!           integral of vmhd2 from v1 to v2
         implicit real*8 (a-h, o-z)
         external vmhd2
         data  epsa/1.d-8/
         call kdexpIntF(vmhd2, v1, v2, epsa, s, err, icon)
         if(icon .ne. 0) then
             write(*,*) ' icon=',icon, ' err=',err
         endif
         vmhi2=s
      end
!            test vmhtx2
!*include mucset
!        implicit real*8 (a-h, o-z)
!       real*4  epsa(4)/0.01, .03333333, .05, .1/
!       call mucset(11.,22., .5, 5.5)
!       do 80 k=1, 4
!           eps=epsa(k)
!           write(*,*) ' eps=',eps, ' int(eps-1)=',
!    *      vmhtx2(eps)
!  80   continue
!       end
        function vmhtx2(v1)
!           integral of vmhd2/v from v1 to 1
         implicit real*8 (a-h, o-z)
         external vmhd2v
         data  epsa/1.d-8/
         v2=1.d0-1.d-5

         call kdexpIntF(vmhd2v, v1, v2, epsa, s, err, icon)
         if(icon .ne. 0) then
             write(*,*) ' icon=',icon, ' err=',err
         endif
         vmhtx2=s
      end
      function vmhd2v(x)
          vmhd2v=vmhd2(x)/x
      end
!           test vmhtx1
!*include mucset
!        implicit real*8 (a-h, o-z)
!       character*70 txt
!       character*16 capx,capy
!       real*4  epsa(4)/0.01, .03333333, .05, .1/
!       open(13, file='c2s5001.#gd.data')
!       capx=' e'
!       capy='integ(fai1)/v(>eps)'
!       call mucset(11.,22., .5, 5.5)
!       do 80 k=1, 4
!           eps=epsa(k)
!           if(k .eq. 1) then
!              write(txt,
!    *          '(''integral eps-1 of fai(inf)/v of had. prod. by mu'',
!    *           '' eps='',f7.4)') eps
!           else
!              write(txt,'('' eps='',f7.4)') eps
!           endif
!           write(13) txt
!           write(13) capx, capy
!           do 100 j=1, 1000
!               e=.01*10.**( (j-1)/20.)
!               if(e .gt. 1000.) goto 101
!                 write(13) sngl(e), sngl(vmhtx1(e, eps))
! 100       continue
! 101       continue
!           write(13) 1.e50, 1.e50
!  80   continue
!       end
      function vmhtx1(e, v1)
!        integral v1 to 1 of  fai1/v
         implicit real*8 (a-h, o-z)
         external vmhd1a
         common /$/ ec
         data  epsa/1.d-8/
         data v2/0.9999999999999d0/
         ec=e

         call kdexpIntF(vmhd1a, v1, v2, epsa, s, err, icon)
         if(icon .ne. 0) then
             write(*,*) ' icon=',icon, ' err=',err
         endif
         vmhtx1=s
      end
      function vmhd1a(x)
         implicit real*8 (a-h, o-z)
         common /$/ ec
         vmhd1a=vmhd1(ec, x)/x
      end
!        make sampling table for hadron like part
!*include mucset
!       implicit real*8 (a-h,o-z)
!       eps=1./ 10.d0
!       open(13,file='c2s5001.#gd.data')
!       call mucset(11.,22., .5, 5.5)
!       v=eps
!       write(*,*)' eps=',eps
!       do 100 i=1,  1000
!          v=v*10.d0**.005d0
!          if(v .gt. .999) goto 101
!          u= vmhtx2(v)/ vmhtx2(eps)
!          epsu=eps**u
!          write(*,*) ' u=',u, ' v=',v
!          write(13) sngl(u), sngl(v/epsu)
! 100   continue
! 101   continue
!     end
!c         to make sampling talbe for (u,v) for mu hadron prod. point
!          like
       include  'mucset.f'
         implicit real*8 (a-h,o-z)
         character txt*70, capx*16, capy*16
         parameter (ie=10,iu=51, iv=1000)
         dimension  vtbl(iu, ie), vu(iv), vw(4*iv),  ua(iv)
         eps=1./ 10.d0
         call mucset(11.,22., .5, 5.5)
         open(13,file='c2s5001.#gd.data')
         open(14,file='c2s5001.#gd2.data')
         open( 7,file='c2s5001.#h.fort(v4)')
         write(*,*)' eps=',eps
         capx='u'
         capy='v/eps**u'
       do   j=1, 50
         e=.01*10.**( (j-1)   )
         if(e .gt. 300.) goto 301
         write(*,*) ' e=',e
         write(txt,'('' eps='',f8.5,'' e='',f7.2)') eps, e
         write(13) txt
         write(13) capx, capy
         write(14) txt
         write(14) capx, capy
         v=1.
         nv=0
          do   i=1, 1000
           v=v/10.d0**0.005
           if(v .le. eps ) goto 101
           u= vmhtx1(e, v)/ vmhtx1(e, eps)
!c         write(*,*) ' u=',u, ' v=',v
           vn=eps**u
           write(14) sngl(u), sngl(v/vn)
           if(u .gt. 0. .and. u .lt. 1.)then
               nv=nv+1
               vu(nv)=v/vn
               ua(nv)=u
           endif
          enddo
  101   continue
!c        make u->v table
        write(13) 0.e0, 1.e0
        iuu=2
         do   u=.020d0, .980001d0, .020d0
           m=5
           epsx=1.d-4
           call daklag(ua, vu, nv, u, m, epsx, f, vw, icon)
           if(icon .ne. 0)then
               write(*,*) ' icon=',icon
           endif
           vtbl(iuu, j) =f
           iuu=iuu+1
           write(13) sngl(u), sngl(f)
         enddo
        write(13) 1.e0, 1.e0
        write(13) 1.e50, 1.e50
        vtbl(1, j)=1.
        vtbl(iu, j)=1.
        write(14) 1.e50, 1.e50
       enddo
  301   continue
        call mkdt('v4  ',vtbl, 2, iu*(j-1),'f8.5,  ', 8,  1)
      end
