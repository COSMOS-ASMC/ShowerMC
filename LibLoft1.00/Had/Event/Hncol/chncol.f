!     this is a trial to suppress very large x part by making an
!     intermediate routine chncol2 which is the original chncol.
!     And on top of it, chncol is made.
!     *****************************************************************
!     *                                                               *
!     * chncol: hadron nucleon collision
!     *                                                               *
!     *****************************************************************

       subroutine chncol(pj, tg, a, ntp, icon)
!             pj:  type ptcl. Input. to give the incident
!     projectile hadron in the Lab. system.
!             tg:  type ptcl. Input. to give the target
!                  nucleon in the Lab. system.
!              a:  array of type ptcl. Output. to get
!                  produced particles.
!            ntp:  integer. Output. to get the total number of produced
!                  particles in a.
!           icon:  Integer. Output. if 0, particle generation is ok
!                   if non 0, you may give up generation.
!
!                 Note: target recoil is put in a(ntp-1),
!                       projectile recoil is put in a(ntp)
!
       implicit none

#include  "Zptcl.h"
#include  "Zevhnv.h"

!
       type(ptcl):: pj, tg,  a(*)
       integer ntp, icon
       integer i
       real*8 u, xmax

!       real*8 BigXRejCnst/1./, BigXRejPw/2./
!         default
        real*8 BigXRejCnst/.4/, BigXRejPw/2.2/
!         enhance
!        real*8 BigXRejCnst/.8/, BigXRejPw/2.2/
!         denhance
!       real*8 BigXRejCnst/.01/, BigXRejPw/1./
!       real*8 BigXRejCnst/0.0003/, BigXRejPw/3.5/

       do while (.true.)
          call chncol2(pj, tg, a, ntp, icon)
          if(icon .ne. 0) goto 100
!             find max secondary energy.
          xmax = 0.
          do i = 1, ntp-2
             if(xmax .lt. a(i)%fm%p(4)) xmax = a(i)%fm%p(4)
          enddo
          xmax = xmax/pj%fm%p(4)
!              if very large x appears, discard it
!              with some probabilty to adjust
!              the x-distribution.
!
!c          if(xmax .gt. 0.4) then
!c             call rndc(u)
!c             if(u .lt.  (1. - (xmax-0.40)/0.6 )**2 ) then
!c                goto 100
!c             endif
!c          else
!c             goto 100
!c          endif
          call rndc(u)
          if(u  .lt. 
     *         (  BigXRejCnst/(BigXRejCnst +
     *         (xmax/(1.0-xmax))**BigXRejPw) ) )  goto 100

       enddo
 100   continue
       end

!     *****************************************************************
!     *                                                               *
!     * chncol2: hadron nucleon collision
!     *                                                               *
!     *****************************************************************
!
!

       subroutine chncol2(pj, tg, a, ntp, icon)
!             pj:  type ptcl. Input. to give the incident
!                  projectile hadron in the Lab. system.
!             tg:  type ptcl. Input. to give the target
!                  nucleon in the Lab. system.
!              a:  array of type ptcl. Output. to get
!                  produced particles.
!            ntp:  integer. Output. to get the total number of produced
!                  particles in a.
!           icon:  Integer. Output. if 0, particle generation is ok
!                   if non 0, you may give up generation.
!
!                 Note: target recoil is put in a(ntp-1),
!                       projectile recoil is put in a(ntp)
!
       implicit none

#include  "Zptcl.h"
#include  "Zevhnv.h"
!
!
       type(ptcl):: pj, tg,  a(*)
       integer ntp, icon, jcon, nfin
       integer i, outc
       type(ptcl):: pjin
!         use pt=0 incident
      pjin = pj
      pjin%fm%p(1) = 0.
      pjin%fm%p(2) = 0.
      pjin%fm%p(3) = sqrt(pjin%fm%p(4)**2 - pjin%mass**2)
!

!        *** until loop**until succeed
      do while (.true.)
!            generate 2 leading ptcls  Zevhnv become ready for use.
          call cs2lp(pjin, tg, icon)
          if(icon .ne. 0) then
!              neglect this event; because of very low energy
!              no particle can be produced
               ntp = 0
              jcon = 0
          else
!               generation of  particles other than the leadings.
              call cgnlp(a, ntp, jcon)
          endif

          if(jcon .eq. 0 ) goto 10
       enddo   
   10 continue
!           ptcls in 'a' should have cms energy, here.
!           give final ptcl charge/ subcode for 
!           pi+-,  k0,k0~,k+,k-
      call cpikcd(a, ntp)
      outc =Rtglab%charge+Rpjlab%charge - pj%charge- tg%charge
      call cconsvChg(outc, a, ntp, icon)
      if(icon .eq. 0) then
!              boost  to lab system
         do   i=1, ntp
            call cibst1(i, Cmsp, a(i), a(i))
         enddo
!          decay of "composit ptcls" (nn~, dd~)
         call cdcycp(a, ntp, nfin)
         ntp=nfin
!            store recoils in a
         a(ntp+1) = Rtglab
         a(ntp+2) = Rpjlab
         ntp = ntp +2
!     -------------- rotate 
         call crot3mom(pj, a, ntp)
      endif

      end
!      ****************************************************
!           decay of nn~  or  DD~ in the projectile system
       subroutine cdcycp(a, nin, n)
!      ****************************************************
       implicit none
!----       include  '../../Zptcl.h'
#include  "Zptcl.h"
!----       include '../../Zcode.h'
#include  "Zcode.h"
       integer nin, n
       type(ptcl):: a(*)
       integer i, k, nx
       type(ptcl):: p(2)
          n=nin
          do i=1, nin
             k=a(i)%code
             if(k .eq. knnb .or. k .eq. kddb) then
                 if(k .eq. knnb) then
                     call cnnbdc(a(i), p, nx)
                 elseif(k .eq. kddb) then
                     call cddbdc(a(i), p, nx)
                 endif
!                 put n or D in the i-th pos. and append c antiptcl 
!                 at the bottome and increase n
                 a(i) = p(1)
                 a(n +1 ) = p(2)
                 n=n+1
             endif    
          enddo
      end
! ******************************************************* 
!       give final ptcl code for pi+,-, k0,k0~,k+,k-
      subroutine cpikcd(a, ntp) 
!
!*******************************************************
      implicit none
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
!----      include '../../Zcode.h'
#include  "Zcode.h"
!----      include '../Zevhnv.h'
#include  "Zevhnv.h"
      integer ntp
      type(ptcl)::  a(ntp) 
!
      integer i, k
      real*8 x
        if(ntp .eq. 1) then
!              nothing to do. charge assignment should have
!              been done in c1pion
        else
             do i=1, ntp
                x = a(i)%fm%p(3)/Pjcms%fm%p(3)
                k = a(i)%code
                if(k .eq. kpion .and. a(i)%charge .ne. 0) then
                     if(x .gt. 0.) then
                         call cspipm(Pjcms, x, a(i))
                     else
                         call cspipm(Tgcms, -x, a(i))
                     endif
                elseif(k .eq. kkaon) then 
!                            set kaon charge
                     if(x .gt. 0.) then
                        call cskchg(Pjcms, x, a(i))
                     else
                        call cskchg(Tgcms, -x, a(i))
                     endif
                endif
             enddo
        endif
      end 
!   ******************************************************* 
!         nn~ decay.  
!        a: projectile.  b: decay product 
!   *******************************************************
      subroutine cnnbdc(a, b, n)
      implicit none
!----       include '../../Zptcl.h'
#include  "Zptcl.h"
!----       include '../../Zcode.h'
#include  "Zcode.h"
       integer n
       type(ptcl):: a, b(2) 
!
       integer ic
       real*8 u 
!
         call rndc(u)
         ic=u*2         ! charge
         if(ic .eq. 0) then
             call cmkptc(knuc,  kneutron, ic,   b(1))
             call cmkptc(knuc,  kneutronb, -ic,  b(2))
         else    
             call cmkptc(knuc, 0, ic,  b(1))
             call cmkptc(knuc, 0, -ic, b(2))
         endif    
         call c2bdcy(a, b(1), b(2))
         n=2
      end 
! ******************************************************* 
!       decay of (dd~) 
!       a: parent.  b: decay product
! *******************************************************
      subroutine cddbdc(a, b, n)
      implicit none
!----       include '../../Zptcl.h'
#include  "Zptcl.h"
!----       include '../../Zcode.h'
#include  "Zcode.h"
       integer n
       type(ptcl):: a, b(2) 
!
       integer ic
       real*8 u 
!
         call rndc(u)
         ic=int(u*3)-1 ! charge
         if(ic .eq. 0) then
             call cmkptc(kdmes,  kd0, ic,  b(1))
             call cmkptc(kdmes,  kd0b, -ic, b(2))
         else    
             call cmkptc(kdmes,  0, ic,  b(1))
             call cmkptc(kdmes,  0, -ic,  b(2))
         endif
         call c2bdcy(a, b(1), b(2))
         n=2
      end 
!***************************************************************** 
!  cspipm: set pi+/pi- 
!
!***************************************************************** 
! the ratio npi+/npi- for p incident case is
!       1+b*x with b=3.5 (x<.6) and 3.1exp(4.4(x-.6)) (x>.6)
!       if  the ratio,   pi+/pi- = f(x),
!          prob of pi+ = (f/(1+f))

       subroutine cspipm(pj, x, a)

       implicit none

#include  "Zptcl.h"
#include  "Zcode.h"

       type(ptcl):: pj, a
       real*8 x
!
       real*8 f, up, u
!
         if(x .lt. .6d0) then
             f = 1.d0+3.33d0*x
         else
             f = 3.0d0* exp( 4.4d0*(x-.6d0))
         endif
         if(pj%charge .eq. 1) then
             up=f/(1.d0+f)
         elseif(pj%charge .eq. -1) then
             up=1.d0/(1.d0+f)
         else
             if(pj%subcode .eq. 0) then
                 up=0.5
             elseif(pj%subcode .eq. regptcl) then
                 up=1.d0/(1.d0+f)
             else
                 up=f/(1.d0+f)
             endif
         endif
         call rndc(u)
         if(u .lt. up) then
              a%charge = 1
              a%subcode = regptcl
         else
              a%charge = -1
              a%subcode = antip
         endif
      end
!     ************************************************
!           set kaon charge/subcode
      subroutine cskchg(pj, x, a)
!     ************************************************
      implicit none

#include  "Zptcl.h"
#include  "Zcode.h"

      type(ptcl):: pj, a
      real*8 x
!
      real*8 f, u, up
       if(a%charge .ne. 0) then
!             k+/k- for p incident
           if(x .lt. .3d0) then
              f=exp(4.36d0*x)
           elseif(x .lt. .6d0) then
              f=3.7d0*exp(6.5d0*(x-.3d0))
           else
              f=27.d0*exp(11.3d0*(x-.6d0))
           endif
           if(pj%charge .eq. 1) then
              up=f/(1.d0+f)
           elseif(pj%charge .eq. -1) then
              up=1./(1.d0+f)
           else
              if(pj%subcode .eq. 0) then
                 up=0.5d0
              elseif(pj%subcode .eq. regptcl) then
                 up=1.d0/(1.d0+f)
              else
                 up=f/(1.d0+f)
              endif
           endif
           call rndc(u)
           if(u .lt. up) then
               call cmkptc(kkaon, 0, 1,  a)
           else
               call cmkptc(kkaon, 0,  -1,  a)
           endif
       else
!            k0
           call rndc(u)
           if(u .lt. .50d0) then
              call rndc(u)
              if(u .lt. 0.5) then
!                 k0 short
                 call cmkptc(kkaon, k0s, 0,  a)
              else
                 call cmkptc(kkaon, -k0s, 0,  a)
              endif
           else
              call rndc(u)
              if(u .lt. 0.5) then
                 call cmkptc(kkaon, k0l, 0, a)
              else
                 call cmkptc(kkaon, -k0l, 0, a)
              endif
           endif
       endif
      end
