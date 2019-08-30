       include  'muhadr.f'
!     compute average muon range
!      **************************************************************
!      *  usage:  if you want to run this program at tss terminal,
!      *          0)    use mujcl0 to make load module library in
!      *                #load.load. (once in a day).
!      *          1)    mylib lib1(#load.load)
!      *                 1) is needed for each session.
!      *          2)    issue the fort command.
!
       external dxdemu
        character ttl*70, capx*16, capy*16
        data nmin/20/, nmax/641/, epsa/1.e-3/, epsr/1.e-3/
        open(13,file='c2s5001.#gd%data')
         epsp=1.0
         epsb=1.0
         epsh=1.0
!c          kgf rock
         z=12.8
         a=25.8
         zba=.493
         z2ba=6.30
         rho=3.02
!c          standard rock
         z=11.
         a=22.
         zba=.5
         z2ba=5.5
         rho=2.8
!c         frejus rock
         z=11.
         a=22.
         zba=.5
         z2ba=5.87
         rho=2.73
!c
         emin=.5e-3
       ttl='muon average range for frejus rock:  these are old g-p xsec'
       capx='e_mu(tev)'
       capy='<r>(1000hg/cm**2)'
       write(13) ttl
       write(13) capx, capy
!c
       call mucset(z, a, zba, z2ba)
       call mupair(epsp)
       call mubrem(epsb)
       call muhadr(epsh)
       e=.1
!       *** until loop*** 
       do while (.true.)
          call aqe(emin, e, dxdemu,  epsa, epsr, nmin, nmax,
     *             s, err, n, icon)
          if(icon .ne. 0) then
              write(*,*) ' e=',e,' icon=',icon, ' err=',err
          endif
          write(13) e, s/1.e5
          e=e*10.**.2
       if         (e .gt. 500.)
     *                    goto 100
       enddo
  100  continue
       write(13) 1.e50, 1.e50
       end
       function dxdemu(e)
!
         call muparf(e, temp)
         dedtp=temp*e
         call mubrmf(e, temp)
         dedtb=temp*e
         call muhadf(e, temp)
         dedth=temp*e
         call mudedx(e, temp)
         dedti=temp
         dxdemu=1./ (dedti+dedtp+dedtb+dedth)
      end
