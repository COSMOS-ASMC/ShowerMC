!      include '../../KKlib/rnd.f'
!      include '../../KKlib/kcossn.f'
!      include 'cpxyzp.f'
!      include 'cmkptc.f'
!c  ----------------------------------
!c       to test cnbdcy
!      include '../Zptcl.h'
!      include '../Zcode.h'
!      implicit none
!      integer  n, i, j, icon
!      parameter(n = 10)
!      type(ptcl):: p(n)
!      real *8  ecm/5./, w, sumx, sumy, sumz, sume
!      j = 1
!      do i=1, 3
!         call cmkptc(kpion, 0, 1,  p(j))
!         call cmkptc(kkaon, k0s, 0, p(j+1))
!         call cmkptc(kpion, 0, 0, p(j+2))
!         j = j+3
!      enddo
!      call cmkptc(komega, 0, 0,  p(n))
!      do j=1, 100
!            call cnbdcy(n, ecm,  p, 0,  w, icon)
!            if(icon .ne. 0) stop 111
!            sumx=0.
!            sumy=0.
!            sumz=0.
!            sume=0.
!            write(*, *) ' ----------w=', w
!            do i=1, n
!c             ---------------------- to draw momentum balance graph
!c                 write(*,*) 0., 0.
!c                 write(*,*) sngl(p(i).fm.p(1)), sngl(p(i).fm.p(3))
!c                 write(*,*)
!c            --------------------------
!c                /////////// to see momentum conservation
!               sumx=sumx + p(i).fm.p(1)
!               sumy=sumy + p(i).fm.p(2)
!               sumz=sumz + p(i).fm.p(3)
!               sume=sume + p(i).fm.p(4)
!c              ///////////////////////
!           enddo
!           write(*,*) sumx, sumy, sumz, sume
!       enddo
!      end
!      ***********************************************************
       subroutine cnbdcy(n, ecm, p, jw,  w, icon)
       implicit none
!      ***********************************************************
!
!        ref:  CPC.  40(1986)p359.  Kleiss, Stirling and Ellis
!
!       n: input.  number of ptcls >=2 (see however, for n=2,
!                  c2bdcy and for n=3, c3bdcy)
!     ecm: input.  cms energy.
!       p: input.  /ptcl/ p(i).mass gives the mass of the i-th ptcl
!                  in the same unit of ecm.
!         output.  /ptcl/ p(i).fm.p(1)
!                                 py
!                                 pz
!                                 e   of the i-th ptcl
!      jw: input.  0--->unweighted event (w=1) obtained. the event
!                  generated need not be discarded.
!                  1--> weighted event( w changes event to event )
!                  the event must be discarded according to the
!                  acceptance probability of w=weight/wax weight).
!      w: output.  see jw
!   icon: output.  0-->event generated successfully
!                  1-->ecm < sum of mass
!                  2-->could not generate (weight problem)
!        With n=10  ( 6 pions, 3 kaons, 1 omega ) Ecm=5 GeV,
!        to get icon=0 for unweighted events,
!        average number of trials is 7;  the distribution is
!        very well approximated  by exp(-ntry/7)dntry
!        With Ecm=10 GeV, <ntry> becomes 2.1
!             Ecm=4 GeV, <ntry> = 13
!      ==============
!        With n=4, (2 pions, 1 kaon, 1 omega) Ecm=4 GeV,
!         <ntry> is almost 1
!----       include '../Zptcl.h'
#include  "Zptcl.h"
       integer n, jw, icon
       type(ptcl):: p(n)
       real*8 ecm, w
!      ------------------
       real*8 mu(1, 1)/0.d0/, wx, w0, wmax, gzai, u
!
       logical ok
       integer nc
!
	nc = 0    ! counter to break inf. loop
!       *** until loop*** 
       do while (.true.)
!            generate massless ptcls isotropically without conservation
           call cnbdc1(n, p)
!             conformal transformation to conserve 4-momentum
           call cnbdc2(n, ecm, p)
!$$$$$$$$$$$
!          call cnbdct(n, p)
!$$$$$$$$$$
!             get gzai to transform massive case
           call cnbdc3(n, ecm, p, mu, 0,  gzai, icon)
!          **********************
           if(icon .ne. 0) return
!          **********************
!             tranform to massive case
           call cnbdc4(n, p,  mu, 0, gzai)
!$$$$$$$$$$$
!          call cnbdct(n, p)
!$$$$$$$$$$
!             compute weight for massive  case
           call cnbdc5(n, ecm,p, wx)
!             compute weight for massless case
           call cnbdc6(n, ecm, w0)
!$$$$$$$$$$$$$$
!          write(*,*) ' wx=',wx,' w0=',w0
!$$$$$$$$$$$
           w=wx*w0
!                compute max possible weight
           call cnbdc7(n, ecm, p,  wmax)
           wmax=wmax*w0
           if(jw .eq. 0) then
!$$$$$$$$$$$$$$
!          write(*,*) ' wmax=',wmax
!$$$$$$$$$$$
!                judge if the event is to be accepted
              call rndc(u)
              if(wmax .eq. 0.d0) then
                 ok=.true.
              else
                 ok = u .lt. w/wmax
              endif
              w=1.
           else
              if(wmax .eq. 0.d0) then
                 w=1.d0
              else
                 w=w/wmax
              endif
              ok=.true.
           endif
           if(ok) goto 100
	   nc = nc +1
	   if(nc .gt. 20) then
	      icon = 2
	      goto 100
           endif
       enddo
  100  continue
       end
       subroutine cnbdct(n, p)
       implicit none
!----       include '../Zptcl.h'
#include  "Zptcl.h"
       integer n
       type(ptcl):: p(n)
!
       real*8 sumx, sumy, sumz, sume
       integer i
!
           sumx=0.d0
           sumy=0.d0
           sumz=0.d0
           sume=0.d0
           do   i=1, n
              sumx=sumx+p(i)%fm%p(1)
              sumy=sumy+p(i)%fm%p(2)
              sumz=sumz+p(i)%fm%p(3)
              sume=sume+p(i)%fm%p(4)
           enddo
           write(*,*) ' sumx,y,z=',sumx, sumy, sumz, ' sume=',sume
       end
       subroutine cnbdc1(n, p)
       implicit none
!            generate massless ptcls isotropically without conservation
!----       include '../Zptcl.h'
#include  "Zptcl.h"
       integer n
       type(ptcl):: p(n)
!
       integer i
       real*8 u1, u2, u, cs, sn, cst, snt
       do   i=1, n
!             *** until loop*** 
             do while (.true.)
                 call rndc(u1)
                 call rndc(u2)
                 u=u1*u2
                if(u .gt. 0.) goto 10
             enddo
   10        continue
             p(i)%fm%p(4) = -log(u)
             call kcossn(cs, sn)
             call rndc(u)
             cst=2*u-1.d0
             snt=sqrt(1. - cst**2)
             p(i)%fm%p(1) = p(i)%fm%p(4)*snt*cs
             p(i)%fm%p(2) = p(i)%fm%p(4)*snt*sn
             p(i)%fm%p(3) = p(i)%fm%p(4)*cst
           enddo
       end
       subroutine cnbdc2(n, ecm, p)
!             conformal transformation to conserve 4-momentum
       implicit none
!----       include '../Zptcl.h'
#include  "Zptcl.h"
       integer n
       type(ptcl):: p(n)
       real*8 ecm
!
       real*8 sumx, sumy, sumz, sume, em, g
       real*8 a, x, bx, by, bz, bq, pe, tmp, px, py, pz
       integer i
!
       sumx=0.d0
       sumy=0.d0
       sumz=0.d0
       sume=0.d0
       do   i=1, n
          sumx=sumx+p(i)%fm%p(1)
          sumy=sumy+p(i)%fm%p(2)
          sumz=sumz+p(i)%fm%p(3)
          sume=sume+p(i)%fm%p(4)
       enddo
       em=sqrt( sume**2 - (sumx**2+sumy**2+sumz**2) )
       g=sume/em

       a=1.d0/(1.d0+g)
       x=ecm/em
       bx=-sumx/em
       by=-sumy/em
       bz=-sumz/em
!
       do   i=1, n
          bq=bx*p(i)%fm%p(1) + by*p(i)%fm%p(2) + bz*p(i)%fm%p(3)
          pe=x*(g*p(i)%fm%p(4) +bq)
          tmp=p(i)%fm%p(4)+a*bq
          px=x*(p(i)%fm%p(1) +   tmp*bx)
          py=x*(p(i)%fm%p(2) +   tmp*by)
          pz=x*(p(i)%fm%p(3) +   tmp*bz)
          p(i)%fm%p(1)=px
          p(i)%fm%p(2)=py
          p(i)%fm%p(3)=pz
          p(i)%fm%p(4)=pe
        enddo
       end
!      ***********************************************
       subroutine  cnbdc3(n, ecm, p, mu, inm, gzai, icon)
!             get gzai to transform massive case
!             put inm=0 if all mu are the same.
!      ***********************************************
       implicit none
!----       include '../Zptcl.h'
#include  "Zptcl.h"
       integer n, inm, icon
       type(ptcl):: p(n)
       real*8 ecm, mu(inm, n), gzai
!
       real*8 eps/1.d-3/, f, fp, fow
       integer nr
!           initial guess of gzai
       gzai=.85d0
       nr=0
!          *** until loop*** 
       do while (.true.)
               call cnbdcf(n, ecm, p,  mu, inm,  gzai, f, fp)
               gzai=  gzai - f/fp
               fow=f
               nr=nr+1
!$$$$$$$$$$$$
!              write(*,*) ' fow=',fow
!$$$$$$$$$$$$
              if(abs(fow) .lt. eps .or. nr .gt. 15) goto 100
       enddo
  100  continue
       if(nr .gt. 15) then
           icon=1
       else
           icon=0
       endif
       end
!      *********************************************
       subroutine cnbdcf(n, ecm, p, 
     *            mu, inm, gzai, f, fp)
!      *********************************************
       implicit none
!----       include '../Zptcl.h'
#include  "Zptcl.h"
       integer n, inm
       type(ptcl):: p(n)
       real*8 gzai, f, fp, mu(inm, n)
!
       real*8 mux, fx, tmp, ecm
       integer i
!
       fx=0.d0
       fp=0.d0
       do   i=1, n
!                  if compiler is good, we can use mu(1,i)
!                  even for inm=0; next is for safty.
              if(inm .eq. 0) then
!                 mux=mu(1,1)
                  mux=0.d0
              else
                  mux=mu(1,i)
              endif

              tmp=  sqrt(p(i)%mass**2+
     *        gzai**2 *( p(i)%fm%p(4)**2 -mux**2 ) )
              fx=fx + tmp
              fp=fp + ( p(i)%fm%p(4)**2- mux**2)/ tmp
           enddo
          f=log(fx/ecm)
          fp=fp*gzai/fx
       end
!      *********************************************
       subroutine  cnbdc4(n, p, mu, inm,  gzai)
!             tranform to massive case
!      *********************************************
       implicit none
!----       include '../Zptcl.h'
#include  "Zptcl.h"
       integer n, inm
       type(ptcl):: p(n)
       real*8 mu(inm, n), gzai

       real*8 mux
       integer i
       do   i=1,n
             p(i)%fm%p(1) = gzai*p(i)%fm%p(1)
             p(i)%fm%p(2) = gzai*p(i)%fm%p(2)
             p(i)%fm%p(3) = gzai*p(i)%fm%p(3)
!                next treatment is for safty
             if(inm .eq. 0) then
!                mux=mu(1,1)
                 mux=0.
             else
                 mux=mu(1,i)
             endif
             p(i)%fm%p(4) = sqrt(p(i)%mass**2 +
     *       gzai**2*( p(i)%fm%p(4)**2-mux**2 ) )
           enddo
       end
       subroutine  cnbdc5(n, ecm, p, wx)
!            compute weig.p(4) for massive case
       implicit none
!----       include '../Zptcl.h'
#include  "Zptcl.h"
       integer n
       type(ptcl):: p(n)
       real*8 ecm, wx
!
       real*8 sum1, pro2, sum3, pab
       integer i
!
          sum1=0.
          pro2=1.
          sum3=0.
          do   i=1, n
             call cpxyzp(p(i)%fm,  pab)
             sum1=sum1+pab
             pro2=pro2* pab/p(i)%fm%p(4)
             sum3=sum3+pab**2/p(i)%fm%p(4)
          enddo
          wx=(sum1/ecm)**(2*n-3)*pro2 /sum3
       end
       subroutine cnbdc6(n, ecm, w0)
!             compute weig.p(4) for massless case
       implicit none
!----       include '../Zptcl.h'
#include  "Zptcl.h"
       integer n
       real*8 ecm, w0
!
       real*8 pi, hpi, gn1
       integer i
       parameter (pi=3.14159265d0, hpi=pi/2)
!
          gn1=1.
          do   i=1, n-2
             gn1=gn1*i
          enddo
          w0=  hpi**(n-1) * ecm**(2*n-4)/(n-1)/gn1/gn1
       end
       subroutine cnbdc7(n, ecm, p,  wmax)
!                compute max possible weig.p(4)
       implicit none
!----       include '../Zptcl.h'
#include  "Zptcl.h"
       integer n
       type(ptcl):: p(n)
       real*8 ecm, wmax

       integer idx(2), nm, i
       real*8 summ, en, beta
!          count massive ptcls
         nm=0
         do   i=1,n
            if(p(i)%mass .gt. 0.d0) then
               nm=nm+1
               if(nm .le. 2) then
                  idx(nm)=i
               endif
            endif
         enddo
         if(nm .eq. 1) then
             wmax=(1. - p(idx(1))%mass/ecm)**(2*n-3)
         elseif(nm .eq. 2)  then
             wmax=
     *        (1. + (p(idx(1))%mass/ecm)**2 -
     *        (p(idx(2))%mass/ecm)**2 )**2
     *       -4*(p(idx(1))%mass/ecm)**2
             if(wmax .le. 0.d0)then
                wmax=1.d-30
             else
                wmax= wmax**(n-1.5d0)
             endif
         else
             summ=0.d0
             do   i=1, n
                 if(p(i)%mass .gt. 0.d0) then
                    en=p(i)%mass/ecm
                    summ=summ+en
                 endif
             enddo
             beta=1. - summ**2
             if(beta .le. 0.d0) then
                wmax=1.d-30
             else
                wmax=sqrt(beta)**(2*n+nm-5)
             endif
         endif
       end
