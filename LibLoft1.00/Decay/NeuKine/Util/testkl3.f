!           kl3               temporary level=2    date=88.07.27
!          k--->pi+mu+neu or pi+e+neu decay, sampling talbe making
!           and distribution function.
!
          real*4  kl3, kl3min, kl3max
          external kl3
!@@@@@@@@@@@@@@@@@@@
          parameter (nnn=1001)
          dimension fa(nnn), ua(nnn), fax(101)
          character ttl*70, capx*16,capy*16
          open(13, file='c2s5001.#gd.data',status='shr',
     *    action='write')
!@@@@@@@@@@@@@@@@@@@
          open(07, file='c2s5001.#h.fort(fbmu)',status='shr',
     *    action='write')
          fmin=kl3min(f)
          fm.p(1)=kl3max(f)
          epsa=1.e-5
          epsr=1.e-5
          nmin=20
          nmax=641
          call aqe(fmin, fm.p(1), kl3, epsa, epsr, nmin, nmax,
     *    snorm, err, nn, icon)
          if(icon .ne. 0) then
             write(*,*) ' icon=',icon
          endif
          write(*,*) ' norm=',snorm
!
!@@@@@@@@@@@@@@@@@@@
          ttl='energy distribution of mu for kc--->pi0,mu,neu;gzai=-.35'
          write(13) ttl
          capx='f=e/mk'
          capy='prob'
          write(13) capx, capy
           do   f=fmin, fm.p(1), .01
              y=kl3(f)
              write(13) f, y/snorm
           enddo
          write(13) 1.e50, 1.e50
!
!@@@@@@@@@@@@@@@@@@@
          ttl='samplign function for e of mu  for kc--->pi0,mu,neu'
          write(13) ttl
          capx='u'
          capy='f<'
          write(13) capx, capy
          write(13) 0., fmin
          i=1
          fa(1)=fmin
          ua(1)=0.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
           do   f=fmin+.001,fm.p(1), .001
             call aqe(fmin, f,    kl3, epsa, epsr, nmin, nmax,
     *       s, err, nn, icon)
             if(icon .ne. 0) then
               write(*,*) ' icon=',icon
             endif
             write(13)    s/snorm, f
!            write(*,*) ' u=',s/snorm, ' f=',f
             i=i+1
             fa(i)=f
             ua(i)=s/snorm
           enddo
          write(13) 1., fm.p(1)
          fa(nnn)=fm.p(1)
          ua(nnn)=1.
          write(13) 1.e50, 1.e50
          ttl='by interpolation'
          write(13) ttl
          write(13) capx, capy
          write(13) 0.,fmin
          fax(1)=fmin
          u=0.01
           do   i=2, 100
             call kfrge(ua, 1, nnn, u, l,icon)
             f=(fa(l)-fa(l-1))/(ua(l)-ua(l-1)) * (u-ua(l-1))
     *          + fa(l-1)
             fax(i)=f
             write(13) u, f
             u=u+.01
           enddo
          fax(101)=fm.p(1)
          write(13) 1., fm.p(1)
          write(13) 1.e50, 1.e50
          call mkdt('fb    ', fax, 1, 101, 'f7.4,  ', 7, 1)
        end
        real function kl3(f)
!           see n.p 22(1961)553-578
! ****************************************************************
!       put ml=object lepton mass (neu or mu)
!           mlp=other lepton  mass
! ****************************************************************
          real*4  mk,  ml, mlp, mmu, me, mpimk
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          parameter ( mk=493.667, mmu=105.66, me=.511,
     *               mpi=134.96,  ml=mmu, alfa=ml/mk,
     *               mpimk=mpi/mk, mlp=0, gz=-.35, a2=alfa**2)
!
         t= mpimk**2/ ( 1.-2*f+ a2)
         kl3= sqrt(f**2-a2)*(1.-t)**2 *(     4* f *(1.-2*f)
     *     + 5*a2*f-a2**2  +gz*a2*(4.-6*f+2*a2)+
     *     gz**2 * a2 * (f-a2)    )
          return
        entry kl3min(f)
            kl3=alfa
            return
        entry kl3max(f)
            kl3=(1.+a2-mpimk**2)/2
        end
