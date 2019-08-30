!     Now this will not be used since sofia is default
!     (HowPhotoP= 1)
!     make next 1 if vector meson after collsion is to be decayed
! else put 0 then, vector meson is replaced by photon
#define VECMESDECAY   1
!  uncomment until   cgpHad and use make -f test.mk
!
! #include "BlockData/cblkGene.h"
!      program main
!#include "ZcosmosExt.h"
!      call testprog
!      end
      

!            test cgpHad
!      subroutine testprog
!        implicit none
!#include  "Zptcl.h"
!#include  "Zmass.h"
!#include  "Zcode.h"

!     integer  massN
!      integer atomicN
!      integer icin
!      integer ntp
!      type(ptcl):: pj
!      integer  nmax
!      parameter (nmax=5000)
!      type(ptcl):: a(nmax)
!      real*8  sumP(4), Eg
!      integer i, j, k

!      massN=14
!      atomicN=7
!      icin = 2
!      call cmkptc(kphoton, 0, 0, pj)
!      write(0,*) 'Enter Eg'
!      read(*,*) Eg
!      pj.fm.p(4)=Eg
!
!      pj.fm.p(1)= 0.
!      pj.fm.p(2)= 0.
!      pj.fm.p(3)=Eg
!      do i = 1, 10000
!         call cgpHad(massN, atomicN, pj, icin, a, ntp)
!         do j= 1, 4
!            sumP(j) = 0.
!         enddo
!         do j = 1, ntp
!            do k = 1, 4
!               sumP(k)  = sumP(k) + a(j).fm.p(k) 
!            enddo
!            write(*,'(2i3, 4g12.3)') a(j).code, a(j).charge,
!     *                      (a(j).fm.p(k),k=1,4)
!         enddo
!         write(*,'(4g12.3)') (sumP(k), k=1,4)
!         write(*,*)
!         write(*,*) 'n= ', ntp-1
!      enddo
      
!      end
!         gamma-n(p or A)-->hadrons
!       This cgpHad is called when whichcode = "current"
!       in cphotop. 
!       "current" means if Eg < 2.5, experimetnal data is used
!        Eg>2.5 GeV, current Active interaction model is 
!        basically used with the following  projectile which is
!        made from incident photon.
!       rho, omega or phi.  if the model can accept one of these
!       pi0                 if not, if the model accept  this
!       pi+/-               if not, use this and leading pi+/- is
!                           replaced by pi0 in the collision prod.

!          a: /ptcl/ output. container of produced ptcls
!        ntp: integer. output. # of produced ptcls.
        subroutine cgpHad(pj, massN, atomicN,  a, ntp)
        implicit none
#include  "Zptcl.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zmass.h"
#include  "Zcode.h"

        integer ntp, icin
        type(ptcl):: pj,  a(*)
        integer  massN, atomicN

        if( pj%fm%p(4) < 2.5 ) then
           call cgpLowExp(pj, massN, atomicN,  a, ntp)
        else
           call cfakeGH(pj, massN, atomicN,  a, ntp)
        endif
        end


      subroutine cgpLowExp(pj, massN, atomicN,  a, ntp)
!          basically Eg < 2.5 GeV.  use exp. data
      implicit none
#include  "Zptcl.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zmass.h"
#include  "Zcode.h"
      type(ptcl):: pj          ! input projectile photon
      integer,intent(in):: massN ! target A
      integer,intent(in):: atomicN ! target Z


      integer,intent(out):: ntp ! # of ptcls produced
      type(ptcl):: a(*)        ! output produced particles


      type(ptcl)::pic

      integer ic                ! target N charge
      logical fermim
      type(ptcl):: tgt
      type(ptcl)::  pjx
      integer icon
      integer jtype
      integer k

      if( massN >= 2 ) then
!           fix target charge (n or p)
         call cfxTgtChg(massN, atomicN, ic)
      else
         ic=atomicN
      endif
!          make target
      call cmkptc(knuc, regptcl, ic, tgt)
      fermim=(pj%fm%p(4) -pj%mass) .lt. Efermi
     *          .and. massN >= 2
      if(fermim) then
         call csampFermiM(tgt%fm) ! 4 mom. has been  set
!              boost the projectile into target
!              rest system (trs).
         call cbst1(1, tgt, pj,  pjx)
      else
         pjx = pj
!            rest target
         tgt%fm%p(1) = 0.
         tgt%fm%p(2) = 0.
         tgt%fm%p(3) = 0.
         tgt%fm%p(4) = tgt%mass
      endif
!             make cm ptcl
      call  cgeqm(pj, tgt, Cmsp, icon) ! not pjx
      if(icon /= 0  ) then
         write(0,*) ' cms cannot be formed in cgpLowExp'
         stop
      endif
!            fix collision type
      call cghCollType(pjx, jtype)
      if(jtype .eq. 0) then     ! will not happen
         ntp=1                  !   older version  0 and no product
         a(1)=pj                ! gamma 
      elseif(jtype .eq. 1) then
!           gp-->p+pi0 or gn-->n+pi0
!           'a' gets particles at target rest system
         call cg1pi0(pjx, ic, a, ntp)
      elseif(jtype. eq. 2) then
!           gp-->n+pi+ or gn-->p+pi-;  at target rest system
         call cg1pic(pjx, ic,  a, ntp)
      elseif(jtype .eq. 3) then
!           gp-->p pi+ pi- or gn --> n pi+ pi- at  CMS
         call cg2pi(ic,  a, ntp) 
      elseif(jtype .eq. 4) then
!                'a' gets CMS ptcls
         call cg3pi(ic, a, ntp)
      elseif(jtype .eq. 5) then ! will not come in our setting
!               vector meson type.  ptcls produced  in lab.
         call cfakeGH(pj, massN, atomicN, a, ntp)
      else
         write(0,*) ' strage jtype=',jtype, ' from cghCollType'
      endif
      if(fermim .and. jtype .le. 2) then
!            boost ptcls back to lab. 
         do   k=1, ntp
            call cibst1(k, tgt, a(k), a(k))
         enddo
      elseif(jtype .eq.  3 .or. jtype .eq. 4) then
!              now in cms. boost to lab
         do k =1, ntp
            call cibst1(k, Cmsp, a(k), a(k))
         enddo
      else
!          jtype =1 or 2 and fermin=F; then a is already in lab.
      endif
      end
!      ****************************************************************
!         fix g--->hadrons interaction type
!        jtype=1   gp-->p+pi0 or gn-->n+pi0
!   
!        jtype=2   gp-->n+pi+ or gn-->p+pi-
!                
!        jtype=3   gp --> p pi+ pi- pi0 or  gn n pi+pi- pi0
!
!        jtype=4   vector meson collision
!        jtype=5
!        jtype=0  no-production
!      ****************************************************************
      subroutine cghCollType(pj, jtype)
      implicit none
#include  "Zptcl.h"
       type(ptcl):: pj
       integer jtype

       real*8 egl, xs1, xs2, xs3, xs4, xso, xst, u
       real*8 xs
       if(pj%fm%p(4) .lt. 5.) then ! actually come here when < 2.5 GeV
!             log10(Eg/MeV); xs in micro barn
          egl=log10(pj%fm%p(4)) + 3
          call cgppi0(egl, xs1)
          call cgppip(egl, xs2)
          call cgppi2(egl, xs3)
          call cgppi3(egl, xs4)
       else
          xs1=0.
          xs2=0.
          xs3=0.
          xs4=0.
       endif
!            gp total x-section  xs in mb
       call cgpxs1(pj%fm%p(4),   xs)
       xs=xs*1000.              ! in micro barn
       xso=max(0.d0, xs-(xs1+xs2+xs3+xs4) ) ! other channel
       if(pj%fm%p(4) .lt. 2.5) xso=0.
       xst=xs1+xs2+xs3+xs4+xso
       if(xst .gt. 0.) then
          call rndc(u)
          if(u .lt. xs1/xst) then
!              gp-->p+pi0 or gn-->n+pi0
             jtype=1
          elseif(u .lt. (xs1+xs2)/xst) then
!                   gp-->n+pi+ or gn-->p+pi-
             jtype=2
          elseif(u .lt. (xs1+xs2+xs3)/xst) then
!                   gp-->p pi+ pi- or gn --> n pi+ pi-
             jtype=3
          elseif( u .lt.  (xs1+xs2+xs3+xs4)/xst) then
             jtype=4
          else
!                  vector meson collision
             jtype=5
          endif
       else
          jtype=0
       endif
       end
!          gn --> resonance production
       subroutine cg1pi0(pj, ic, a, ntp)
       implicit none
#include  "Zptcl.h"
#include  "Zmass.h"
#include  "Zcode.h"
#include  "Zevhnv.h"

       type(ptcl):: pj, a(*)
       integer ic, ntp
!

        real*8 cs, tmass
        type(ptcl):: eres
        save
!
        tmass=masp
!                   gp-->p+pi0 or gn-->n+pi0; sample cos of pi0 in cms
        call csPiAngOfPiN(Cmsp%mass, 1, 0, cs)
!          resonance energy in trs
        eres%fm%p(1) = 0.
        eres%fm%p(2) = 0.
        eres%fm%p(4) = pj%fm%p(4) + tmass
        eres%mass = Cmsp%mass
        eres%fm%p(3) = sqrt(eres%fm%p(4)**2 - eres%mass**2)
        call cmkptc(kpion, 0, 0, a(1))        
        call cmkptc(knuc, regptcl, ic, a(2))
        call c2bdcp(eres, a(1), cs, a(2))
!        call c2bdcp(Cmsp,  a(1), cs, a(2))
        ntp=2
        end
!       **************
        subroutine cg1pic(pj, ic, a, ntp)
       implicit none
#include  "Zptcl.h"
#include  "Zmass.h"
#include  "Zcode.h"
#include  "Zevhnv.h"

       type(ptcl):: pj, a(*)
       integer ic, ntp
!


        real*8 cs, tmass
        type(ptcl):: eres
        save

!       **************
!                   gp-->n+pi+ or gn-->p+pi-; sample cos of pi in cms
        tmass = masp
        call csPiAngOfPiN(Cmsp%mass, 0, 1, cs)
        eres%fm%p(4)=pj%fm%p(4) + tmass
        eres%mass = Cmsp%mass
        eres%fm%p(3) = sqrt(eres%fm%p(4)**2 - eres%mass**2)
        call cmkptc(kpion, 0, ic, a(1))
        call cmkptc(knuc, regptcl, (1-ic)/2, a(2))
        call c2bdcp(eres, a(1), cs, a(2))
!        call c2bdcp(Cmsp, a(1), cs, a(2))
        ntp=2
        end
!       **************
        subroutine cg2pi(ic, a, ntp)
!       **************
!          particles are produced in cms.
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnv.h"

       type(ptcl):: a(*)

       integer ic, ntp

        real*8 w
        integer icon
!                   gp-->p pi+ pi- or gn --> n pi+ pi-
       
       call cmkptc(knuc, regptcl, ic, a(1))
       call cmkptc(kpion, 0, 1, a(2))
       call cmkptc(kpion, 0, -1, a(3))

       call cnbdcy(3, Cmsp%mass, a,  0, w, icon)
       if(icon .eq. 1) then
          write(0, *)
     *    ' cnbdcy fails in gp-->p pi+ pi- ', 
     *    ' roots=',Cmsp%mass, ' icon=',icon
          ntp=0
       else
          ntp=3
       endif
       end
!       **************
       subroutine cg3pi(ic, a, ntp)
       implicit none
#include  "Zptcl.h"
#include  "Zmass.h"
#include  "Zcode.h"
#include  "Zevhnv.h"

       type(ptcl):: a(*)
       integer ic, ntp, icon
       real*8 w
!
!       **************
!                   gp-->p pi+ pi- pi0 or gn-> 3pi
!            in cms.
       call cmkptc(knuc, regptcl, ic, a(1))
       call cmkptc(kpion, 0, -1, a(2))
       call cmkptc(kpion, 0, 0, a(3))
       call cmkptc(kpion, 0, 1, a(4))
       call cnbdcy(4, Cmsp%mass, a,  0, w, icon)
       if(icon .eq. 1) then
           write(0,*) ' cnbdcy fails in gp--> p + 3pi ',
     *     ' roots=',  Cmsp%mass, ' icon=',icon
           ntp=0
       else
!          icon =2 comes here. no problem statistically.
!            few percent cases for mass=1.6 to 3 GeV happens to be icon=2
!            (icon = 2 means rejection after 20 trials due to wight problem)
          ntp=4
       endif
       end
!      ************************************************************
!         neutral  meson collision.
       subroutine cfakeGH(pj, massN, atomicN, a, ntp)
!       use modXsecMedia
       use modColInfo
       implicit none
#include  "Zcode.h"
#include  "Zptcl.h"   
#include  "Zevhnv.h"



       type(ptcl):: pj   ! input. photon
       integer,intent(in):: massN ! target  A
       integer,intent(in):: atomicN ! target Z
       type(ptcl):: a(*)  ! produced ptcls
       integer,intent(out)::ntp  ! # of ptcls produced

!
       type(ptcl):: vm
       integer jcon

       real(8)::u
       real(8)::xs
       integer nout
       integer::pichg
       type(ptcl):: pix

!           make pi+ or -    
       call rndc(u)
       
       pix = pj  

       if(ActiveMdl == "qgsjet2") then
          pichg = 0
          call cmkptc(kpion, 0, 0, pix)  ! can accept pi0
          call cadjm(pix, pix)  ! adjust momenutm
          call cxsecQGS(pix, massN,   xs )  
          call chAcol(pix, massN, atomicN, xs, a, ntp)
       elseif( ActiveMdl == "epos") then
          call cinelx(pix, massN, atomicN, xs)
          call chAcol(pj, massN, atomicN, xs, a, ntp)
       elseif (ActiveMdl /= "ad-hoc" ) then
          if(u < 0.5 ) then
             pichg = -1
          else
             pichg = 1
          endif
          call cmkptc(kpion, 0, pichg, pix)
          call cadjm(pix, pix)
!             some model(JAM) needs xs. for safety get xs 
          call cinelx(pix, massN, atomicN, xs)
          call chAcol(pix, massN, atomicN, xs, a, ntp)
          call cLeadingPiAfterCol(pix, a, ntp)
       else
          call cmkVectorMeson(pj, vm, jcon)
          if(jcon /= 0) then
             write(0,*) "cmkVectorMeson failed"
             ntp=1
             a(1) = pj
          else
             call cinelx(pix, massN, atomicN, xs)
             call chAcol(vm, massN, atomicN, xs, a, ntp)
             call cVecMesonAfterCol(vm, a, ntp, nout)
             ntp = nout
          endif
       endif
       
       end

      subroutine  cmkVectorMeson(pj,  vm, jcon)
      implicit none
#include  "Zptcl.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zmass.h"
#include  "Zcode.h"
      type(ptcl):: pj          ! input. photon
      type(ptcl):: vm          ! output
      integer,intent(out):: jcon !

      real(8):: p, alfa
!         fix vector meson (rho, omega, or phai)
!                              46  46         8 %
      call cfixVectorMeson(pj%fm%p(4), vm, jcon)
!     make vector meson proj.
      p=sqrt(pj%fm%p(4)**2 - vm%mass**2)
      alfa=p/pj%fm%p(4)
      vm%fm%p(1) = pj%fm%p(1)*alfa
      vm%fm%p(2) = pj%fm%p(2)*alfa
      vm%fm%p(3) = pj%fm%p(3)*alfa
      vm%fm%p(4) = pj%fm%p(4)
      end

      subroutine cVecMesonAfterCol(vm, a, nin, nout)
      implicit none
#include  "Zptcl.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zmass.h"
#include  "Zcode.h"
      type(ptcl):: vm          ! input. vector meson
      integer,intent(in):: nin  ! # of ptcls in a
      type(ptcl):: a(nin)      ! ptcls generated by col.
      integer,intent(out)::nout ! after vm treatment, # of ptcls in a
      type(ptcl):: b(10)
      integer i, nx, j
      real(8):: p, alfa
      nout  = nin
      do i = 1, nin
         if( a(i)%code == vm%code ) then
#if VECMESDECAY == 1
            call cvmdcy(a(i), b, nx)
            a(i) = b(1)
            do j = 2,  nx
               a(j+nout-1) = b(j)
            enddo
            nout = nout + nx -1
#else
            a(i)%code = kphoton
            a(i)%mass = 0.
            p=a(i)%fm%p(4)
            alfa=sqrt(dot_product( a(i)%fm%p(1:3),a(i)%fm%p(1:3)))
     *         /p
            a(i)%fm%p(1:3) = a(i:3)%fm%p(1)/alfa
#endif
         endif
      enddo
      end
      subroutine  cLeadingPiAfterCol(pix, a, ntp)
      implicit none
!         replace max energy pi with same type of pix
!      by pi0
#include  "Zptcl.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zmass.h"
#include  "Zcode.h"
      type(ptcl):: pix
      integer,intent(in):: ntp
      type(ptcl):: a(ntp)
      
      integer i, maxi
      real(8)::maxE
      maxE=-1.0
      maxi =0
      do i = 1, ntp
         if( pix%code == a(i)%code ) then
            if( pix%charge == a(i)%charge ) then
               if(maxE < a(i)%fm%p(4)) then
                  maxE = a(i)%fm%p(4)
                  maxi = i
               endif
            endif
         endif
      enddo
      if( maxi > 0 ) then
         call cmkptc(kpion, 0, 0, a(maxi))
         call cadjm(a(maxi),a(maxi))
      endif
      end
!      *****************************************
       subroutine cfixVectorMeson(e, vm, icon)
!      *****************************************
       implicit none

#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zmass.h"
       real*8 e
       type(ptcl):: vm
       integer icon
!       
       integer nc
       real*8 u, amass, w
!
       nc=0
!         *** until loop*** 
       do while (.true.)
          nc=nc+1
          call rndc(u)
          if(u .lt. .46) then
             call cmkptc(krho, 0, 0, vm)
             w=wrho
          elseif(u .lt. .92) then
             call cmkptc(komega, 0, 0, vm)
             w=womega
          else
             call cmkptc(kphi, 0, 0, vm)
             w=wphai
          endif
!              *** until loop*** 
          do while (.true.)
             call ksbwig(vm%mass, w, amass)
             if (amass .gt. vm%mass-w .and. amass .lt. vm%mass+w)
     *                           goto 10
          enddo
 10       continue
          if(e .le. amass) then
             icon=1
          else
             icon=0
          endif
          if  (icon .eq. 0 .or. nc .gt. 10)
     *                      goto 100
       enddo
       vm%mass = amass
 100   continue
       
       end
!      *****************************************************************
!          make decay of a vector meson
!      *****************************************************************
       subroutine cvmdcy(vm, a, np)
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
       type(ptcl):: vm, a(*)
       integer np
!
       if(vm%code .eq. krho) then
          call crhodc(vm, a, np)
       elseif(vm%code .eq. komega) then
          call comgdc(vm, a, np)
       elseif(vm%code .eq. kphi) then
          call cphidc(vm, a, np)
       endif
       end
