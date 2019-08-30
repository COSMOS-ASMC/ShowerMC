!         
!         
       subroutine cgnlp(a, ntp, icon)
!     ********************************************************
!            a: /ptcl/    Output. to get produced particles in cms
!          ntp:  output.  Integer. The number of particle in a.
!         icon:  Integer. Output.0--> o.k
!                                1--> n.g
       implicit none
#include  "Zptcl.h"
#include  "Zmass.h"
#include  "Zevhnv.h"
!
       type(ptcl):: a(*)
       integer ntp, icon
!      -------------------------------- 
       real*8   missgm, roots
       integer jcon
!


!       integer maxiso/10/  ! max # of ptcl from isotropic p.s   
       integer maxiso/0/  ! max # of ptcl from isotropic p%s   
!       logical ok
!
!              sample # of charged ptcls
!
       missgm = Missingp%mass
       roots =  Cmsp%mass
!        *** until loop*** 

       do while (.true.)
!                get average # and sampled # of charged ptcl
!                Avncharged & Nch
          call csnchp(jcon)
          if(jcon .ne. 0) then
             icon=1
             goto 900           ! return
          endif
!              fix # of pi+-, pi0, k+-,k0, nn~, dd~,
!              and give ptcl mass and code in a;
!              # of pi+- etc is obtained by calling
!              cqnptc(code, charge, nout)
          call cfnptc(a,  ntp)
           if(ntp .ge. 1 ) goto 10
!           if(ntp .eq. 1 .and. missgm .lt. 15.*maspic) goto 10
!           if(ntp .eq. 2 .and. missgm .lt. 30.*maspic) goto 10
       enddo
 10    continue
       if(ntp .le. maxiso .and. ntp .ge. 2) then
!               use isotripc p.s
            call ciso(ntp, a, icon)
            if(icon .eq. 1) then
!                   <Pt> is too large so use cylindrical p.s
               call ccylps(ntp, a, icon)
            elseif(icon .eq. 2) then
                icon=1
            endif
        elseif(ntp .ge. 2) then
           call ccylps(ntp,  a, icon)
        else
!              assume 1 (or 0) pion production in cms
            call c1pion(a, ntp, icon)
        endif
  900  continue
       end
!     *****************************
      subroutine ccpmul(roots, avn)
!     *****************************
!      average charged particle multiplicity at root(s)
!     /**** UA5 parameterization ***/
!      root(s)  GeV
       implicit none
#include "Zevhnp.h"

       real*8 roots, avn
!       real*8 lambda/0.3/, no/0.6135/, qcdErg/1000./  ! up to 1000 GeV.
!                              | this is wrong (total  multiplicity)
        real*8 lambda/0.3/, no/0.34/, qcdErg/1000./  
!                           
!               e UA5 data.
!                   original
!       avn= 7.2* roots**(2*0.127) -7.   including  leading.
!
         avn = (no* exp(sqrt(23./18. * log( (roots/lambda))*2))
     *          -3.5)*0.8 
         if(avn .lt. 0.1) then
            avn = 0.1
         endif
!
!
      end
!      **************ccylps**************************************
!           generation of ptcls by cylindrical p.s
!         icon:  0---> o.k
!                1---> n.g
!      ****************************************************
       subroutine ccylps(ntp, a, icon)
       implicit none
!----       include  '../../Zptcl.h'
#include  "Zptcl.h"
!----       include  '../Zevhnv.h'
#include  "Zevhnv.h"
       integer ntp, icon
       type(ptcl)::  a(ntp)
!
       real*8 ptav, w
       integer i
!
!          assign pt 
!              init for const.
         call caspti
!                            loc pt = loc pz
         call caspt(a, ntp)
!          pt---> ptx,pty
         call csptxy(a, ntp)
!            forced conservation of pt.
         call cptcns(a, ntp, ptav)
!           generation of rapidity for missing mass
         w = Missingp%mass
!                             loc tm = loc pz
!                             loc rap = loc e
         call cgrap(w, ptav, ntp,  a, icon)
         if(icon .eq. 0) then
!             convert y into to missing mass system energy
            call cytoe(a, ntp)
!           ________________________________
!           call cchk(' after cytoe', a, ntp)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                boosting ntp ptcls to cms
             do   i=1, ntp
                call cibst1(i, Missingp, a(i), a(i))
             enddo
         endif
       end
!      ********************************************
!            generate 1 pion
!            leadings are reset here.
!          icon:  0--> o.k
!                 1--> n.g
       subroutine c1pion(a, ntp, icon)
!         Method:
!          Now, the recoil leadings and missing mass and its
!          4 momentum are given.  We change this missing mass
!          to be a pion and make a transformaton so that
!          3 particles satify the 4 momentum conservation.            
!       a(1): /ptcl/. Outupt. produced pion
!        ntp: integer. Output.  1  if icon =0 (produced)
!                               0  if icon =1 (n.g)
!       icon: integer. Output.  0--> ok. 1--> n.g
!      ********************************************
       implicit none
!----       include '../../Zptcl.h'
#include  "Zptcl.h"
!----       include '../../Zmass.h'
#include  "Zmass.h"
!----       include '../Zevhnv.h'
#include  "Zevhnv.h"
#include  "Zcode.h"

!
       integer ntp, icon
       type(ptcl)::  a(*)

       type(ptcl):: p3(3)   ! to store 3 particles
       real*8 amu(3), roots, gzai
       integer  csum0, csum1, charge
       character*80  msg 
!
!       modify E,P following cpc p.364
!       get gzai for q~=gzai*p~
!       q0=sqrt(m**2 + gzai**2(p0**2-mu**2) )
!     
        p3(1) = Rpjcms   ! reooil proj.  in cms
        p3(2) = Rtgcms   ! recoil trgt.  in cms
        p3(3) = Missingp ! missing mass in cms
        roots = Cmsp%mass  

!           present mass
         amu(1)=Pjlab%mass
         amu(2)=Tglab%mass
         amu(3)=Missingp%mass  ! missing mass
!            

         csum0 = Pjlab%charge + Tglab%charge
         csum1 = Rpjcms%charge +Rtgcms%charge
         if( abs(csum1-csum0) .gt. 1) then
!                retry once more
            icon =1 
            goto 900
         else
            charge = csum0 - csum1
!             true mass    ! we want missing --> pion mass.
!             modify missing mass to be the pion mass
            call cmkptc(kpion, 0, charge, p3(3))
         endif   
!             get convesion factor
!          note: in cnbdcy, amu is not used, because 
!                it is always 0. (start from 0 mass ptcl)
!                and transform to massive case (true mass
!                is in /ptcl/ data.)
!                However, here, the trial mass is missing mass 
!                and must be given in amu. True mass is
!                in /ptcl/ p3
         call cnbdc3(3, roots, p3, amu, 1,
     *   gzai, icon)
         if(icon .eq. 0) then
             call cnbdc4(3, p3,  amu, 1, gzai)
!            ________________________________
!            write(*,*) ' p3 sumpx',
!    *       (p3(1).fm.p(1)+p3(2).fm.p(1)+p3(3).fm.p(1))
!            write(*,*) ' p3 sumpy',
!    *       (p3(1).fm.p(2)+p3(2).fm.p(2)+p3(3).fm.p(2))
!            write(*,*) ' p3 sumpz',
!    *       (p3(1).fm.p(3)+p3(2).fm.p(3)+p3(3).fm.p(3))
!            write(*,*) ' p3 sum e (gev)',
!    *       ( p3(1).fm.p(4)+p3(2).fm.p(4)+p3(3).fm.p(4))
!            write(*,*) ' roots=',roots
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               reset leadings
             Rpjcms = p3(1)
             Rtgcms = p3(2)
!                boost  proj. recoil into lab (cms -> lab)
!                note:  tg(7).fm = cms 4 mom. in lab
             call cibst1(1, Cmsp, Rpjcms, Rpjlab)
!                boost target recoil into lab
             call cibst1(2, Cmsp, Rtgcms, Rtglab)
!               boost proj to target rest system. (may not be used)
             call cbst1(1, Tglab, Rpjlab, Rpjtatr)
!               same for projectile rest system. (may not be used)
             call cbst1(1, Pjlab, Rtglab, Rtgpatr)
!                move pion
             a(1) = p3(3)
             ntp = 1
         else
!            ____________________________________________
             write(0,*) ' failed to adjust  missing mass=',
     *       Missingp%mass, ' into pion mass. '
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             ntp = 0
             icon = 1
         endif
 900     continue
       end
!     *********************************************************
      subroutine caspti
!     *********************************************************
!          preparation for assigning  pt
      implicit none
!----        include '../../Zptcl.h'
#include  "Zptcl.h"
!----        include '../../Zmass.h'
#include  "Zmass.h"
!----        include '../Zevhnv.h'
#include  "Zevhnv.h"
!
      real*8  efe0, xx
!      real*8  ptbase/330.d-3/
      real*8  ptbase/180.d-3/
!      real*8  ptbase/170.d-3/
!
!              some sort of effective E0, which I cannot
!           remember the reasoning now.
!        efe0 = (Efrs**2 + 2*(masn + Pjlab.mass)*Efrs)/2
!     *          /masn
        efe0 =max(1.d0,  (Efrs**2 - masn*2- Pjlab%mass**2)/2
     *          /masn)

!
!              Powerexp in   pt**pw2 * exp(-pt/..) dpt case.
        Powerexp=1.0d0 + (efe0)**(-0.05d0)
        if(efe0 .lt. 40.) then
           xx = (efe0 - 40.)/5.5 + 6.
           Ptnorm = ptbase *exp(xx)/(1+exp(xx))
        elseif(efe0 .lt. 40000.) then
           Ptnorm= ptbase
        else
           Ptnorm= ptbase
!           Ptnorm = ptbase *( efe0/40000.)**0.025
        endif
        Probpower=min(0.04d0*Nch, 0.33d0)  
        Powerp=3.d0 +
     *        1.d0/
     *     (0.01 + 0.01*3.29* Missingp%mass**0.3)
        if(Powerp .gt. 100.d0) then
             Probpower=0.d0
        endif
       end
!
!      *****************
       subroutine caspt(a, ntp)  
!           generation of pt
!      *****************
!----        include '../../Zptcl.h'
#include  "Zptcl.h"
!----        include '../../Zcode.h'
#include  "Zcode.h"
!----        include '../Zevhnv.h'
#include  "Zevhnv.h"
        integer ntp
        type(ptcl):: a(ntp)

        real*8  pttmp
        integer ntpc
!           <>pt for NN~,          DD~
!        real*8 ptnnb/500.d-3/, ptddb/750.d-3/, ptavpi/330.d-3/,
!     *         ptavk/479.d-3/, ptaveta/480.d-3/
        real*8 ptnnb/250.d-3/, ptddb/370.d-3/, ptavpi/180.d-3/,
     *         ptavk/240.d-3/, ptaveta/240.d-3/
!        real*8 ptnnb/240.d-3/, ptddb/360.d-3/, ptavpi/170.d-3/,
!     *         ptavk/230.d-3/, ptaveta/230.d-3/
!
        ntpc = 0
!                 (n,n~)
        pttmp=ptnnb* Ptnorm/ptavpi
!                     sample pt
        call cspt(pttmp, Nnnb, a, ntpc)

!                (d,d~)
        pttmp = ptddb* Ptnorm/ptavpi
        call cspt(pttmp, Nddb, a, ntpc)
!                 pi+/-
        pttmp=Ptnorm
        call cspt(pttmp, Npic,  a, ntpc)
!                 pi0
        call cspt(pttmp, Npi0,  a, ntpc)
!                eta
        pttmp = ptaveta*Ptnorm/ptavpi
        call cspt(pttmp, Neta, a, ntpc)
!                 kaon +/ -
        pttmp = ptavk*Ptnorm/ptavpi
        call cspt(pttmp, Nkch, a, ntpc)
!                 k0
        call cspt(pttmp, Nk0,  a, ntpc)
       end
!      ****************************************************
       subroutine ciso( ntp,  a, icon)
!           isotropic phase space
!          icon: 0---> o.k
!                1---> <pt> is big,   better to use cylindrical p.s
!                2---> n.g
!      ****************************************************
       implicit none
!----       include  '../../Zptcl.h'
#include  "Zptcl.h"
!----       include  '../Zevhnv.h'
#include  "Zevhnv.h"
       integer ntp, icon
       type(ptcl)::   a(ntp)
!
       integer nc, jcon, i
       real*8 ptc1/500.d-3/, ptc2/700.d-3/
!
       real*8 ptlm1, ptlm2, missgm, wg, sumpt, avpt
!
       character*80  msg

       missgm = Missingp%mass
       nc=0
       if(Pjlab%fm%p(4) .lt. 1000.d0) then
             ptlm1=ptc1*(Pjlab%fm%p(4)/1000.d0)**0.04
             ptlm2=ptc2*(Pjlab%fm%p(4)/1000.d0)**0.04
       else
             ptlm1=ptc1
             ptlm2=ptc2
       endif
       if(3.1415/4.0* missgm/ntp .gt. ptlm1) then
             jcon=1
       else
!             *** until loop*** 
             do while (.true.)
!                  isotropic p.s decay; take almost all weights
!                  wg=.95 does not make any diff.(slower speed)
!                call cnbdcy(ntp, missgm, a, 1, wg, jcon)
                call cnbdcy(ntp, missgm, a, 0, wg, jcon)
                nc=nc+1
                if(nc .gt. 15) then
                    jcon=3
                    write(msg,*) ' cnbdcy fail but try further'
                    call cerrorMsg(msg, 1)
                endif
                if(jcon .ne. 0 .or. wg .gt. .05d0) goto 99
             enddo
   99        continue
             if(jcon .eq. 0) then
!                 boost 4 momentum into cms
                 do i=1, ntp
                    call cibst1(i, Missingp, a(i), a(i))
                 enddo
!                   get <pt> and see if it is too large
!                        get sum of  pt
                 sumpt=0.d0
                 do   i=1, ntp
                     sumpt=sumpt +
     *                  sqrt(a(i)%fm%p(1)**2+a(i)%fm%p(2)**2)
                 enddo
                 avpt=sumpt/ntp
                 if(avpt .gt. ptlm2) then
                    jcon = 2
                 elseif(avpt .gt. ptlm1) then
                    jcon = 1
                 else
                    jcon = 0
                 endif
             else
                 jcon = 2
             endif
         endif
         icon=jcon
      end
!     ***********************************************
!         convert rapidity y into to cms energy
!
      subroutine cytoe(a, n)
      implicit none
!----       include  '../../Zptcl.h'
#include  "Zptcl.h"
       integer n
       type(ptcl):: a(n)
!
       integer i
       real*8 etemp
          do   i=1, n
!              note that rap and e have the same pos.
!             etemp=a(i).fm.tm *  cosh( a(i).fm.rap  )
             etemp=a(i)%fm%p(3)*  cosh( a(i)%fm%p(4)  )
!             a(i).fm.p(3) = a(i).fm.tm * sinh( a(i).fm.rap )
             a(i)%fm%p(3) = a(i)%fm%p(3)* sinh( a(i)%fm%p(4))
             a(i)%fm%p(4) = etemp
          enddo
       end
! *****************************************************
      subroutine cspt(avpt,  nptcl, a, ntpc) 
! to sample a pt 
! avpt: real*8. Input. average pt in exponential part.  
! nptcl: integer. Input. # of ptcls to be assigned
!     a: /ptcl/.  output. a.fm.p(3) is given a pt value
!  ntpc: integer. input/oututp. a(ntpc+1) is the first
!                 ptcl pos.  ntpc is incremented by nptcl
!                 on return.

        implicit none
!----        include '../../Zptcl.h'
#include  "Zptcl.h"
!----        include '../Zevhnv.h'
#include  "Zevhnv.h"
        integer nptcl, ntpc
        type(ptcl):: a(*)
        real*8 avpt
!
!        real*8 u, bpt/1.7d0/, pt
!        real*8 u, bpt/1.5d0/, pt
        real*8 u, bpt/1.0d0/, pt


        integer nc
        
        do nc = 1, nptcl
           call rndc(u)
           if(u .lt. Probpower) then 
!                  power pt
                 call cspwpt(bpt, Powerp, pt)
           else 
!                   pt**pw2 * exp(-pt/..)dpt

                call ksgmrm(Powerexp, avpt, pt)
           endif
           ntpc = ntpc +1
           a(ntpc)%fm%p(3) = pt
        enddo
        end
!       ***********************************************************
! 
!         dptcns: do forced conservation of pt * c
!       ***********************************************************
! 
!/usage/ call cptcns(a, nt, ptav) 
! 
! 1) compute sum of ptx and pty and distribute them to each ptcl 
!    propotionally to pt to have zero-sum ptx and pty.  
! 2) adjust pt, ptx, pty so that the sum of pt becomes the same 
!    values as that of original one 
! 
! a(nt) : /ptcl/. Input/Output.  
! ptav: real*8. Output.  <pt> 
! 
! 
! before calling this routine, pz should not be 
! set as z-component. it is the location of
! pt.
       subroutine cptcns(a, nt, ptav)
       implicit none
!----       include '../../Zptcl.h'
#include  "Zptcl.h"
       integer nt
       type(ptcl):: a(nt)
       real*8 ptav 
!
       real*8 sumpx, sumpy, sumpt, cfx, cfy
       real*8 sumpt2, cf
       integer i 
!
      if(nt .ge. 2) then 
!           get sum of pt, ptx, pty
         sumpx = 0.d0
         sumpy = 0.d0
         sumpt = 0.d0
         do i=1, nt
            sumpt = sumpt + a(i)%fm%p(3)
            sumpx = sumpx + a(i)%fm%p(1)
            sumpy = sumpy + a(i)%fm%p(2)
         enddo
         if(sumpt .gt. 0.d0) then 
!              correction factor
            cfx=sumpx/sumpt
            cfy=sumpy/sumpt 
!
            sumpt2=0.d0
            do i=1, nt
               a(i)%fm%p(1) = a(i)%fm%p(1) - a(i)%fm%p(3) * cfx
               a(i)%fm%p(2) = a(i)%fm%p(2) - a(i)%fm%p(3) * cfy
               a(i)%fm%p(3) = sqrt(a(i)%fm%p(1)**2 + a(i)%fm%p(2)**2 )
               sumpt2=a(i)%fm%p(3) + sumpt2
            enddo 
!                adjust: sum (pt) = original value
            cf = sumpt/sumpt2 
!                  multipliy cf to pt and ptx,pty
            do i=1, nt
                a(i)%fm%p(3) = a(i)%fm%p(3)*cf
                a(i)%fm%p(1) = a(i)%fm%p(1)*cf
                a(i)%fm%p(2) = a(i)%fm%p(2)*cf
            enddo
            ptav=sumpt/nt
        else
            ptav = 1.d-1
        endif
      elseif(nt .gt. 0) then
         ptav=a(1)%fm%p(3)
      else
         ptav = 1.d-1
      endif
      end 
