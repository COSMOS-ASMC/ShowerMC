!        *********************************************************
!        *
!        * mupair0:
!        *
!c         to make sampling talbe for (u,v) for mu pair.
!       implicit real*8 (a-h,o-z)
!       eps=1./100.d0
!       open(13,file='c2s5001.#gd.data')
!       v=eps
!       dp=5.1d-3
!       write(*,*)' eps=',eps, ' total xsection/b(pair)=',
!    *  (prmui(1.)-prmui(eps)),   eps*(eps+dp)**2
!       do 100 i=1,  1000
!          v=v*10.d0**.005d0
!          if(v .gt. .999) goto 101
!          u=(prmui(1.)- prmui(v))/ (prmui(1.)-prmui(eps))
!          write(*,*) ' u=',u, ' v=',v
!          write(13) sngl(u), sngl(v)
!          if(u .gt. .95) then
!              alfa=(prmui(1.)-prmui(eps))*(1.-u)*eps*(eps+dp)**2
!              vv=eps+alfa
!              write(*,*) '             v=',vv
!          elseif(u .lt. 2.e-5) then
!              alfa=(prmui(1.)-prmui(eps))*u*(1.+dp)**2
!              vv=1.-alfa
!              write(*,*) '             v=',vv
!          endif
! 100   continue
! 101   continue
!     end
!        test muon pair function *v**2
!     open(13,file='c2s5001.#gd.data')
!     dp=5.1d-3
!     do 100 i=1, 1000
!          v=0.001*10.**( (i-1)/40. )
!          if(v .gt. .9999) goto 101
!          write(13) v, dp*(1.+dp)/(v+dp)**2/v  * (v**2)
! 100 continue
! 101 continue
!      end
       function prmui(v)
!c
!           indefinite integral of pair cre function of muon
         implicit real*8 (a-h,o-z)
          data dp/5.1d-3/
          prmui=  (1.d0+dp)*(1.d0/(v+dp) - log( (v+dp)/v )/dp)
      end
