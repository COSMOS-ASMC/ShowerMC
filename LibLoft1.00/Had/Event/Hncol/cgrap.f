!      ****************************************************** 
       subroutine cgrap(w,  ptav, ntp,  a, icon)
!          generation of rapidity for mass missing. mass
!      ****************************************************** 
       implicit none
!----       include  '../../Zptcl.h'
#include  "Zptcl.h"
       integer ntp, icon
       type(ptcl):: a(ntp)
       real*8  ptav, w   ! w is missing mass
!
       integer i, maxi, mini
       real*8 y, z
!                 sample proto-rapidity
         call cprap(ptav,  a,  ntp)
!                 normalize proto-rapidity to 0 to 1.
         call cnprap(a, ntp, maxi, mini)
!                 compute transverse mass.(note; save in pt pos.)
         do   i=1, ntp
!              a(i).fm.tm= sqrt(a(i).mass**2 + a(i).fm.p(3)**2)
              a(i)%fm%p(3)= sqrt(a(i)%mass**2 + a(i)%fm%p(3)**2)
         enddo
!              get coef. for y and z to modify rapidity
!              to conserve 4 momenta.
         call ccmrap(w, a, ntp, maxi, mini, y, z, icon)
         if(icon .eq. 0) then
!               convert to true rapidity satisfing 4 mom.
!               conservation.
             call cctrap(a, ntp, y, z)
!            ____________________________________________________
!                         check conservation
!            call cccrap(a,  ntp, sume, sump)
!            write(*,*) ' sume=',sume,' m  =',pj.mass,  ' sump=',sump
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         endif
       end
!      *******************************************************
       subroutine cctrap(g, n, y, z)
!           convert rapidity into true rapidity satisfing e-p
!           conservation.
!      *******************************************************
       implicit none
!----       include '../../Zptcl.h'
#include  "Zptcl.h"
       integer n
       type(ptcl):: g(n)
       real*8 y, z
!       
       integer i
!
          do   i=1, n
!               g(i).fm.rap = g(i).fm.rap*y + z
               g(i)%fm%p(4) = g(i)%fm%p(4)*y + z
          enddo
       end
!      *******************************************************
       subroutine cccrap(g,  n, sume, sump)
!                 check conservation
!      *******************************************************
       implicit none
!----       include '../../Zptcl.h'
#include  "Zptcl.h"
       integer n
       type(ptcl):: g(n)
       real*8 sume, sump
!
       integer i
       real*8 yr
!
        sume=0.d0
        sump=0.d0
        do   i=1, n
!               yr=g(i).fm.rap
                yr=g(i)%fm%p(4)
!               sume=sume+ cosh(yr)*g(i).fm.tm
!               sump=sump+ sinh(yr)*g(i).fm.tm
               sume=sume+ cosh(yr)*g(i)%fm%p(3)
               sump=sump+ sinh(yr)*g(i)%fm%p(3)
        enddo
       end
!      *********************************************
       subroutine cnprap(g,  n, maxi, mini)
!            normalize proto-rapidity in 0 to 1       
!      *********************************************
!       g(n): /ptcl/ Input.  g(i).fm.rap has y and
!             is normalized to have values in (0, 1.0).
!             the max value becomes 1.0 and the minimum one 0.0
!         n : integer Input. # of ptcls in g
!       maxi: integer. Output. g(maxi).fm.rap has  max y(=1)
!       mini: integer. Output. g(mini).fm.rap has  min y(=0)
!
       implicit none
!----       include '../../Zptcl.h'
#include  "Zptcl.h"
       integer n, maxi, mini
       type(ptcl):: g(n)
!
       integer i
       real*8 gmx, gmn
       maxi=1
       mini=1
       do   i=1, n
!          if(g(i).fm.rap .gt. g(maxi).fm.rap ) maxi=i
          if(g(i)%fm%p(4) .gt. g(maxi)%fm%p(4) ) maxi=i
!          if(g(i).fm.rap .lt. g(mini).fm.rap ) mini=i
          if(g(i)%fm%p(4) .lt. g(mini)%fm%p(4) ) mini=i
       enddo
!       gmx=g(maxi).fm.rap
!       gmn=g(mini).fm.rap
       gmx=g(maxi)%fm%p(4)
       gmn=g(mini)%fm%p(4)
       do   i=1, n
!          g(i).fm.rap = (g(i).fm.rap - gmn )/(gmx-gmn)
          g(i)%fm%p(4) = (g(i)%fm%p(4) - gmn )/(gmx-gmn)
       enddo
      end
       subroutine ccmrap(w,  g,  n, maxi, mini, y, z, icon)
!        get coefficients y and z to modify rapidities so that
!        the total energy and pz should be conserved.
!
!
!       w: real*8. input.  available energy
!    g(n): /ptcl/  input.  g(i).fm.tm and g(i).fm.rap are used.
!       n: integer. input.  number of particles
!    maxi: integer. input. g(maxi).fm.rap is the max y
!    mini: integer. input. g(mini).fm.rap is the min y
!       y: real*8   output. coefficient in rap <- z + y* rap'
!                   where rap is the true rapidity which
!                   satisfy the energy-momentum conservation.
!       z: real*8   output. coefficient in rap <- z + y* rap'
!    icon: integer  output. 0--> o.k
!                           1--> n.g. retry
!  see. cpc. vol9. (1975). 297 by Jadach.
!
!          function to be solved for y is  symmetric around y=0
!          and has a form like  f(y)= c - y**2 ( c> 0)
!          (y > 0 is obtained as a solution)
!          in some stragne input, c becomes <0 with no solution
        implicit none
!----        include '../../Zptcl.h'
#include  "Zptcl.h"
        integer n, icon
        type(ptcl):: g(n)
!
        real*8 w, y, z
        integer maxi, mini
!
        real*8 eps1/0.0010d0/, w2, alw2, y1, epsx
        real*8 sump, summ, sumgp, sumgm, tmp
        real*8 expgyp, expgym, fy1, fy1p, dy, eps
        integer lp, i
!
         w2=w*w
         alw2=log(w2)
!         y1=log(w2/g(maxi).fm.tm/g(mini).fm.tm)
         y1=log(w2/g(maxi)%fm%p(3)/g(mini)%fm%p(3))
         epsx=eps1/sqrt(dble(n))
!
         lp=0
!         *** until loop*** 
         do while (.true.)
              lp=lp+1
              y = y1
              sump=0.
              summ=0.
              sumgp=0.
              sumgm=0.
              do   i=1, n
!                  tmp=g(i).fm.rap*y
                  tmp=g(i)%fm%p(4)*y
                  if(tmp .gt. 100.d0 .or. 
     *                  tmp .lt. -100.d0) then
!                        no solution case: c < 0
                      lp=100
                      goto 100
                  else

                      expgyp=exp(tmp)
                      expgym=1.d0/expgyp
!                      sump = sump+ g(i).fm.tm*expgyp
!                      summ = summ+ g(i).fm.tm*expgym
                      sump = sump+ g(i)%fm%p(3)*expgyp
                      summ = summ+ g(i)%fm%p(3)*expgym
                      sumgp = sumgp +
!     *                     g(i).fm.tm*g(i).fm.rap*expgyp
     *                     g(i)%fm%p(3)*g(i)%fm%p(4)*expgyp
                      sumgm=sumgm +
!     *                     g(i).fm.tm*g(i).fm.rap*expgym
     *                     g(i)%fm%p(3)*g(i)%fm%p(4)*expgym
                  endif
               enddo

              fy1=alw2 -log(sump*summ)
              fy1p= - sumgp/sump + sumgm/summ
              dy =fy1/fy1p
              y1=y- dy
              eps=dy/y1
              if(abs(eps) .lt. epsx .or. lp .gt. 10) goto 100
         enddo
  100    continue
         if(lp .gt. 10) then
             icon=1
         else
             icon=0
             y=y1
             z=log(w/sump)
         endif
       end
!      *************************************
       subroutine cprap(ptav, pc, n)
!      *************************************
!            generate n proto-rapidities
!   ptav: real*8. input. avarage pt of this event
!  pc(n): /ptcl/  input/output. pc(i).fm.rap is created.
!      n: input. # of particles
!
!
!
       implicit none
!----       include  '../../Zptcl.h'
#include  "Zptcl.h"
       integer n
       type(ptcl):: pc(n)
       real*8 ptav
!
       integer i
       real*8 y, b, u
!
       do   i=1, n
!             b=(pc(i).fm.p(3)/ptav)**(-0.1) * 1.5   ! goot at 900 GeV

!             b=(pc(i).fm.p(3)/ptav)**(-0.2) * 1.5  !  goot at 900 GeV
!             b=(pc(i).fm.p(3)/0.3)**(-0.2) * 0.5   ! good at 200 GeV !
!              this  simple one is best.!!!
!             b= 1.  !  good at almost every where
             b = 1.0

!             call cprap0(b, y)

             call rndc(u)
             y =min(u*(ptav/pc(i)%fm%p(3))**0.40d0,  2.5d0)
             call rndc(u)
             if(u .lt.  0.5) then
                y = -y
             endif
                
!             pc(i).fm.rap = y
             pc(i)%fm%p(4) = y
       enddo
       end
!     ***************************
       subroutine cprap0(a, y)
!         proto rapidity sampling
!     ***************************
       implicit none
       real*8 a, y
!
!             !
!           1 !*********!
!             !         !  *
!             !         !    *
!             !         !      *
!             !         !        *
!             !         !          *
!             0~~~~~~~~~1~~~~~~~~~~a+1
!
       real*8 s, u
          s=1.+a/2
          call rndc(u)
          if(u .lt. 1.d0/s) then
             call rndc(u)
             y=2*u-1.d0
          else
             call rndc(u)
             y=a*(1.d0 -sqrt(u)) + 1.d0
             call rndc(u)
             if(u .lt. 0.5d0) then
                y=-y
             endif
          endif
       end
