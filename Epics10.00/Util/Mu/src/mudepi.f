!       muon depth intesity relation calculation.
!
!          change the next include to one of murngf, murngk, murngs
!          for frejus, kgf, standard rock.
!    note:  for frjus, every module is ready.  for kgf and s.r,
!          avrang etc are for frejus. and usable only at low energies
!          as this program uses.
       include  'murngk.f'
       include  'mumin.f'
!
        external func
        character*70  ttl
        character*16  capx, capy
        real eba( 3)/10., 31.6, 100./
        dimension ra(200), di(200)
        common /$/ b1, b2, eb, gb, r
        logical gb
        open(13,file='c2s5001.#gd.data')
        gb=.true.
        capx='depth(1000hg/cm2)'
        capy='intesity'
         do   b1=2.50, 2.701,.1
             do   b2=2.90, 3.2010, .1
                do   ib=1, 3
                   eb=eba(ib)
                   write(ttl,'(''kgf    rock: b1='',f4.2,'' b2='',
     *             f4.2, '' eb='',f5.0)') b1, b2, eb
                   write(13) ttl
                   write(13) capx, capy
                   ir=0
                    do   r= .25, 17., .1
                      call mumin(r*100000., emin)
!                     write(*,*)' r=',r, ' emin=',emin
                      call muinte(func, emin, ans)
                      ir=ir+1
                      ra(ir)=r
                      di(ir)=ans
!                     write(13) r, ans
                    enddo
!                  sum=0.
                    do   ic=ir-5, 1, -1
                       call trap(ra(ic), di(ic), ir-ic+1, s, icon)
                       if(icon .ne. 0) then
                            write(*,*) ' error'
                       endif
                       write(13) ra(ic), s
                    enddo
                   write(13) 1.e50, 1.e50
                   write(*,*)  ' b1=',b1,' b2=',b2, ' s=',s
                enddo
             enddo
         enddo
        end
        subroutine muinte(f, emin, ans)
        common /$/ b1, b2, eb, gb, r
        logical gb
           data nmin/20/, nmax/641/, epsa/1.e-4/, epsr/1.e-4/
!
                sum=0.
            if(emin .gt. 1.5) then
                e2=emin
!                *** until loop*** 
                do while (.true.)
                    e1=e2
                    e2=e1*4.
                    call aqe(e1, e2, f, epsa, epsr, nmin, nmax,
     1                       s, err, n, icon)
                    if(icon .gt.25000) then
                        write(*,*) ' icon=',icon, ' n=',n, ' err=',err
                    endif
                    sum=sum+s
                    if(sum .eq. 0.) then
                       eps=1.
                    else
                       eps=s/sum
                    endif
                if         (eps .lt. 1.e-4)
     *                             goto 100
                enddo
  100           continue
            else
                ex=avrtoe(log10(r))
                gp=avrngp(ex)
                sum=fd1ry(b1, b2, eb, gb, ex)/gp
            endif
            ans=sum
       end
!
                    function fd1ry( b1, b2, eb, gbend, e )
!
      logical gbend
!
      if( b1 .eq. b2 ) then
          fd1ry=e**(-b1-1.)
      elseif(gbend) then
          fd1ry=e**(-b1-1.) * (1. + e/eb)**(b1 - b2)
      elseif( e .lt. eb) then
          fd1ry=e**(-b1-1.)
      else
          fd1ry=eb**(b2-b1) * e**(-b2-1.)
      endif
      end
      function func(e)
        common /$/ b1, b2, eb, gb, r
        logical gb
        x=log10(e)
        call murang(x, r, f)
        func=fd1ry(b1, b2, eb, gb, e)*f
      end
