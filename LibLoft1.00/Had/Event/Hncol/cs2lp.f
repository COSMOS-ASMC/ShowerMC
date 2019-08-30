!      *******************************************************
!            generate 2 leading ptcls
!      *******************************************************
      subroutine cs2lp(proj, trgt, icon)
!           proj: /ptcl/ Input. projectile in lab
!
!           trgt: /ptcl/ Input. target in lab
!           icon :   integer. Output. if 0, ok.
!                             if non 0, sampling failed after
!                             20 trials. or energy is too low
!      *** Note ***
!    After this call, leading particle infomation is set in
!    variables in ../Zevhnv.h.
!  Projectile 
       implicit none

#include  "Zptcl.h"
#include  "Zmass.h"
#include  "Zevhnv.h"

       logical first/.true./
       save first



!      ***************
       external cblkEvhnp           ! block common name
!      ****************
       type(ptcl):: proj, trgt
       integer icon
!
       type(fmom):: gc
       real*8 xpmin, xpmax, xtmin, xtmax, c2, dest1,
     *        dest2, den

       real*8 maslimit2   ! missing mass is too small or not    
       integer maxtry, count
       parameter (maxtry = 20, maslimit2 =(maspic*1.1)**2 )
       character*70 msg
       logical noferm    ! if target is at rest, no Fermi momentum.

       if(first) then
          call cinippx  ! make sampling table for pp->p+x
          call cinippxn  !  pp->n+x
          call cinipipx
          first = .false.
       endif
!
       count = 0
       icon = 0
       Pjlab = proj
       Tglab = trgt
       noferm = trgt%fm%p(4) .eq. trgt%mass


!          get cms equivlent mass and 4 momentum  
       call cgeqm(Pjlab, Tglab, Cmsp, icon)
       if(icon .ne. 0) then
          write(msg, *) ' cms cannot be formed in cs2lp; proj and ',
     *                'target are '
          call cerrorMsg(msg, 1)
          call cprptc(Pjlab, 1)
          call cprptc(Tglab, 1)
          stop 9999
       endif   
!          get Lorentz factor of cms
       call cgetlf(Cmsp,  gc)
!          boost pj into cms.
       call cbst0(1, gc, Pjlab, Pjcms)     
!            boost target into cms
       call cbst0(2, gc, Tglab, Tgcms)
!          boost proj into target rest system
       if(noferm) then
           Pjtatr = Pjlab
       else    
           call cbst1(1, Tglab, Pjlab, Pjtatr)       
       endif
!            boost target into projectile rest system
       call cbst1(1, Pjlab, Tglab, Tgpatr )
!            get possible max and min x of leading particles
       call cgextx(xpmin, xpmax, xtmin, xtmax)
       if(xpmin .ge. xpmax .or. xtmin .ge. xtmax) then
          icon = 1
          return
       endif

!      ----------------------------------------------
!      *** until loop*** until virtual ptcl that balances
!                              two outgoing leadings become
!                              timelike. 
      do while (.true.)
!          *** until loop***  generation of projectile and target leading
!             ptcls.
          do while (.true.)
!               sample 1 leading ptcl and set it in 
!                  Rpjtatr(target at rest).  Note Rpjtatr is
!                  should be rotated later   
              call cslp(Pjtatr, xpmin, xpmax, Rpjtatr)
!               sample recoil target  and set it in Rtgpatr
!               the same note as above.
              call cslp(Tgpatr, xtmin, xtmax, Rtgpatr)
!                       next should not be used.
!                    some dirty trick to make the strange evnet less
!                    (may not be needed)
!                  make target pt colinear with projectile Pt
!                  while keeping the magnitude as it is
!              c2=Rtgpatr.fm.p(1)**2 + Rtgpatr.fm.p(2)**2
!              den=sqrt(( Rpjtatr.fm.p(1)**2 +  Rpjtatr.fm.p(2)**2)/c2 )
!              Rtgpatr.fm.p(1) = Rpjtatr.fm.p(1)/den
!              Rtgpatr.fm.p(2) = Rpjtatr.fm.p(2)/den
!     If the following two call's were omitted and equivalent ones
!     were placed inside 'cslp', 
!         Absoft compile fials to compile it; it shows
!/home01/kasahara/f/cosmos/sun/Particle/Event/Hncol/cs2lp.f -o cs2lp.o
! error on line 256 of /tmp/temp10493.f: synch error in intermediate code
!         This is not an  ordinary error.  Resolution is give by
!         putting a dummy line relating to rcord /ptcl/Rtgpart
!                rotate   recoils so that they are seen in
!                in a frame where Pjatr or Tgpatr is seen.
             call crot3vec(Pjtatr%fm, Rpjtatr%fm, Rpjtatr%fm)
             call crot3vec(Tgpatr%fm, Rtgpatr%fm, Rtgpatr%fm)
!               next is a dummy substitution to avoid stupid Absoft
!               compiler error.
! **          Rpjtatr = Rpjtatr
!                boost it to lab
             if(noferm) then
                Rpjlab = Rpjtatr
             else
                call cibst1(1, Tglab, Rpjtatr, Rpjlab)
             endif   
!                   boost to cms             
             call cbst1(1, Cmsp, Rpjlab, Rpjcms)
!                  energy libarated by projectile in cms
             dest1= Pjcms%fm%p(4) - Rpjcms%fm%p(4)
!
!                boost to lab 
             call cibst1(1, Pjlab, Rtgpatr, Rtglab)
!                   boost to cms
             call cbst1(1, Cmsp, Rtglab, Rtgcms)
!                  energy libarated by target in cms
             dest2=Tgcms%fm%p(4) - Rtgcms%fm%p(4)
             if(dest1 .gt. maspic .or. dest2 .gt. maspic)  goto 5
             count = count + 1
             if( count .gt. maxtry) then
                icon =1
                goto 5
             endif   
          enddo
    5     continue
!           form a missing mass particle
          Missingp%fm%p(1) = - (Rpjcms%fm%p(1) + Rtgcms%fm%p(1))
          Missingp%fm%p(2) = - (Rpjcms%fm%p(2) + Rtgcms%fm%p(2))
          Missingp%fm%p(3) = - (Rpjcms%fm%p(3) + Rtgcms%fm%p(3))
          Missingp%fm%p(4) = Cmsp%mass - Rpjcms%fm%p(4) - Rtgcms%fm%p(4)
          Missingp%mass = Missingp%fm%p(4)**2
     *     -(Missingp%fm%p(1)**2 + Missingp%fm%p(2)**2 +
     *       Missingp%fm%p(3)**2)
          if(Missingp%mass .lt. maslimit2 ) then
               count = count + 1
               if(count .gt. maxtry)then
                  icon = 1
                  goto 10
               endif   
          else
               Missingp%mass = sqrt (Missingp%mass)
               goto 10
          endif
      enddo
   10 continue
      end
!      ****************************************
       subroutine cgextx(xpmin, xpmax, xtmin, xtmax)
!        get extream of recoil x, defined as the ratio 
!        of incoming and outgoing leading partilce,
!        where the counter  particle is at rest.
!  xpmin:  real*8. Output. minimum x of projectile.
!  xpmax:   //             maximum //
!  xtmin:  real*8. Output. minimum x of target
!  xtmax:   //             maximum //
!
!        min is when projectile after coll. is at rest in cms.
!        max is when projectile after coll. loses mass of 1 pion
!
       implicit none
!

#include  "Zptcl.h"
#include  "Zmass.h"
#include  "Zevhnv.h"
!
       real *8  xpmin, xpmax, xtmin, xtmax
!
       type(ptcl):: rest   ! resting particle

       type(ptcl):: temp, temp2, temp3
!            
       rest%fm%p(1)=0.
       rest%fm%p(2)=0.
       rest%fm%p(3)=0.
       rest%mass=Pjlab%mass
       rest%fm%p(4) = rest%mass
!        min of projectile. 
!           boost stopped proj in cms into lab.
       call cibst1(1, Cmsp, rest, temp)
       temp%mass = rest%mass
!         boost it to target rest system
       call cbst1(1, Tglab, temp, temp2)
       xpmin= temp2%fm%p(4)/Pjtatr%fm%p(4)
!        max
!         get proj. cms energy - mass of pion
       temp=Pjcms
       temp%fm%p(4) =max(temp%fm%p(4) - maspic, Pjlab%mass)
       call cadjm(temp, temp)  ! adjust momentum along with E
!         boost it into lab
       call cibst1(1, Cmsp, temp, temp2)
       temp2%mass=Pjlab%mass
!         boost to target rest system
       call cbst1(1, Tglab, temp2, temp3)
       xpmax= temp3%fm%p(4)/Pjtatr%fm%p(4)
!                       ............
!        max and min x  of target
!            min
!          boost stopped  target in cms into lab system
       rest%mass = Tglab%mass
       rest%fm%p(4) = rest%mass
       call cibst1(1, Cmsp, rest, temp)
       temp%mass = Tglab%mass
!          boost to projectile rest system
       call cbst1(1, Pjlab, temp, temp2)
       xtmin = temp2%fm%p(4)/Tgpatr%fm%p(4)
!            max
!          get cms energy - mass of pion
       temp = Tgcms
       temp%fm%p(4) = max(temp%fm%p(4) - maspic, Tglab%mass)
!          boost it in cms  into lab
       call cibst1(1, Cmsp, temp,  temp2)
       temp2%mass = Tglab%mass
!          boost it to projectile rest system.
       call cbst1(1, Pjlab, temp2, temp3)
       xtmax =  temp3%fm%p(4)/Tgpatr%fm%p(4)
      end
!         

!     *****************************************************************
!     *                                                               *
!     * cslp:    leading ptcl sampling
!     *                                                               *
!     *****************************************************************
!
!
       subroutine cslp(p, akmin, akmax,  a)
!           p:  type ptcl. Input.  Particle
!               given at the rest system of the counter ptcl.
!        akmin:  real*8. Input. min of x of the leading ptcl
!        akmax:  real*8. Input. max of x o//
!            a:  type ptcl. Output.  sampled leading ptcl.
!          Note that the momentum of "a" is defined in
!          a system whose z-axis is the direction of p.fm
!          so that you have to rotate it after calling this,
!          if p has non-zero x, y component of momentum.
!
        implicit none


#include  "Zptcl.h"
        type(ptcl):: p, a
        integer nc, icon
        real*8 xp, avpt, ptn, tmsq, u, akmin, akmax
        logical notfirst  
!


        a = p
!
        nc=0
!       *** until loop*** 
        do while (.true.)
!            sample leading ptcl pt: avpt. output <pt>
!                                    ptn.  output sampled pt

           call cslppt(p, avpt,  ptn)


!            sample leading particle xp with  pt
           tmsq=ptn**2 + p%mass**2
           call cslpx(p, tmsq, akmin, akmax,  xp, notfirst, icon)


           if(icon .eq. 0 .and.  xp-akmin .lt. .2 ) then
!                  xp is small; reject some large pt 
               if(ptn .gt. avpt) then
                   call rndc(u)
                   if(u .gt. avpt/ptn) then
                       icon=1
                   endif
               endif
           endif
           nc=nc+1
           if(icon .eq. 0 .or. nc .gt. 20) goto 5
        enddo
    5   continue
        if(nc .gt. 20) then
           call cerrorMsg(' nc>20 in cslp', 0)
        endif
        a%fm%p(4)=p%fm%p(4)*xp
!           set pt tentatively in pt
        a%fm%p(3) = ptn
!           convert it to ptx, pty
        call csptxy(a,  1)
!               set pz
        a%fm%p(3) = sqrt(a%fm%p(4)**2 - a%mass**2 - ptn**2)
!         fix chacge after collision
        if(notfirst) then
!           keep the same charge if the 2nd,3rd coll. inside A.
!           a.charge = p.charge  ! not needed since a = p
        else
           call cfclp(p, xp, a)
        endif
!          this may be needed if crot3vec is not called
!          after cslp; 
!        call crot3vec(p.fm, a.fm, a.fm)
      end
      subroutine cxtuln(x, ux)
!           get normalized integral (from 0 to x) for given x
!           of leading ptcl (pp-->p)
!           u for x=0 to 1 step .01
         implicit none
         integer i
         real*8 x, ux

#include "Zcinippxc.h"

!
         i=x*nx+1.

         if(i .eq. n) then
             ux=1.
         else
!            ux=(intendndx(i+1)-intendndx(i))*nx * (x - (i-1)*dx)
!             + intedndx(i)
             ux=(intendndx(i+1)-intendndx(i))* (x*nx - i+1) +
     *           intendndx(i)
         endif
      end
      subroutine cxtulnpi(x, ux)
!           get normalized integral (from 0 to x) for given x
!           of leading pi (pp-->p)
!           u for x=0 to 1 step .01
         implicit none
         integer i
         real*8 x, ux

#include "Zcinippxc.h"

!
         i=x*nx+1.

         if(i .eq. n) then
             ux=1.
         else
             ux=(intendndx2(i+1)-intendndx2(i))* (x*nx - i+1) +
     *           intendndx2(i)
         endif
      end
!     *****************************************************************
!     *                                                               *
!     * cfclp:   fix charge of a leading particle
!     *                                                               *
!     *****************************************************************
!                            =   =   =   =
!
        subroutine cfclp(pj, xp, p)
!
        implicit none
!
!----        include  '../../Zptcl.h'
#include  "Zptcl.h"
!----        include  '../../Zcode.h'
#include  "Zcode.h"
!----        include  '../Zevhnp.h'
#include  "Zevhnp.h"
!
      logical pnchgex
      common /Zchgex/ pnchgex
        type(ptcl):: pj, p
        real*8 xp
!
        real*8 rf, u
        integer k0
!c        character*70  msg
!        
        k0=pj%code
        call rndc(u)
!            branch by ptcl kind
        if(k0 .eq. kpion) then
!                      pion; more inelastic one is 
!                           more likely chargeexchanged
           rf=sqrt(1.-xp)
!            if(u .gt. Cepic0*rf) then
!          if(u .gt. 0.3* sqrt(rf) ) then
          if(u .gt. 0.35* rf ) then
!                no charge exc.
             p%charge = pj%charge
          else
             if(pj%charge .eq. 0) then
!                   0--> + or -
                call rndc(u)
                if(u .lt. .5) then
                   p%charge = 1
                else
                   p%charge = -1
                endif
             else
!                 charge-->0 or opposite charge
                call rndc(u)
                if(u .lt. rf*0.30) then
                   p%charge = -pj%charge
                else
                   p%charge = 0
                endif
             endif
          endif
      elseif(k0 .eq. kkaon) then
!                                         kaon
         rf=sqrt(1.-xp)
         if(u .gt. 0.35*rf) then
              p%charge = pj%charge
         else
              p%charge = 0
              call rndc(u) 
              p%subcode = pj%subcode
         endif
      elseif(k0 .eq. knuc) then
!                                        nucleon
         if( .not. pnchgex )  then
!              same charge
              p%charge = pj%charge
         else
              if(pj%charge .eq. 0) then
                  if(pj%subcode .eq. regptcl) then
                      p%charge = 1
                  else
                      p%charge = -1
                  endif
              else
                  p%charge = 0
                  p%subcode = pj%subcode
              endif
         endif
      elseif(k0 .eq. krho) then
         p%charge = 0
      elseif(k0 .eq. komega)then
         p%charge = 0
      elseif(k0 .eq. kphi) then
         p%charge = 0
      elseif(k0 .eq. keta) then
         p%charge = 0
      else
!         write(msg,*) ' code=',k0,' undef. in cfclp'
!         call cerrorMsg(msg, 1)
!           same charge as input
      endif
      end
      subroutine cslpx(pj,  tmsq, akmin, akmax,  x, notfirst, icon)
!           sampling of x 
!            pj: type ptcl. Input.
!          tmsq: input.incident transverse mass square after collision.
!         akmin: input. min x allowed
!         akmax: input. max x allowedn
!             x: output.  sampled x
!       notfirst: output.   becomes t if this is 2nd, 3rd coll. in A
!          icon:  0  x sampled
!                 1  x not sampled. kinematically impossible.
!    **** note ***  If the target is a nucleus and the collision is
!    2nd, 3rd , ...  times inside the nucleus, the x distribution is
!    changed to  x**SucPw dx to have smaller inelasticity.  
!    (SucPw=1.5 is default;
!     this corressponds to alfa=2.5 to Date et al's paper.
!    (PRD1985,vol.32. 619)  This should be 
!    managed by calling cslpx2
!
      implicit none
#include "Zcode.h"
#include  "Zptcl.h"
#include  "Zmass.h"
#include "Zcinippxc.h"
#include "Zevhnp.h"

       logical pnchgex
       common /Zchgex/ pnchgex

       type(ptcl):: pj
       real*8 tmsq, x,  akmin, akmax
       real*8 umin, umax, temp1, temp2
       integer i, icon
       real*8 u, uc
       logical lessInela/.false./, makeless, notfirst
       save lessInela
!

         if(pj%fm%p(4)**2 .le. tmsq) then
             icon=1
         elseif(.not. lessInela) then
!             cxtuln(x0, ans) ; ans= integral of dn/dx from 0, x0
            if(pj%code .ne. knuc) then
               call cxtulnpi(akmin, umin)
               call cxtulnpi(akmax, umax)
            else
               call cxtuln(akmin, umin)
               call cxtuln(akmax, umax)
            endif
!                  uniform random number should be between
!                  umin and umax

            call rndc(u)
            u=(umax-umin)*u + umin
            i=u*nx +1
            if(pj%code .eq. knuc) then
!              
               call rndc(uc)
               if(uc .lt. Ceneuc) then
                  pnchgex= .true.
                  x=(ppsxn(i+1) - ppsxn(i))*nx*(u- (i-1)*dx)
     *              + ppsxn(i)
               else   
                  x=(ppsx(i+1) - ppsx(i))*nx*(u- (i-1)*dx)
     *              + ppsx(i)
                  pnchgex= .false.
               endif
            else
               x=(pipsx(i+1) - pipsx(i))*nx*(u- (i-1)*dx)
     *              + pipsx(i)
            endif
         else
            call rndc(u)
            if(pj%code .ne. knuc) then
!                for mesons, make more inelastic
               temp1 = SucPw + 0.5
            else
               temp1 = SucPw + 1.
            endif
            temp2 = akmin**temp1
            x = ( (akmax**temp1 - temp2 )*u + temp2 )**(1./temp1)
         endif
         if((pj%fm%p(4)*x)**2 .le. tmsq) then
            icon=1
         else
            icon=0
         endif
         notfirst = lessInela
         return
!     ************ call this before 2nd, 3rd coll. inside nucleus
!                  with .true.  and  after that, call with .false.
      entry  cslpx2(makeless)
!     *************
         lessInela = makeless
      end
!     *****************************************************************
!     *                                                               *
!     * cslppt:  samples leading ptcl pt                              *
!     *                                                               *
!     *****************************************************************
!
!
!
        subroutine cslppt(pj, avpt,  ptn)
        implicit none
!
!       pj: strucutre /ptcl/. Input. Projectile particle at
!                                    the rest system of target.
!     avpt: real*8.  Output.  average pt at this energy.
!      ptn: real*8.  Output.  sampled pt in GeV.
!
!----       include  '../../Zptcl.h'
#include  "Zptcl.h"
!
       type(ptcl):: pj
       real*8 avpt, ptn, pw
!
      avpt=226.d-3* pj%fm%p(4)**0.1d0         ! energy is GeV

      pw=2.59d0/pj%fm%p(4)**0.1d0 
!          pt**pw * epx(-pt)dpt  type
      call ksgmrm(pw, avpt,  ptn)
      end
