#include "ZcheckPHITS.h"
      subroutine cphits(pj0, ia0, iz0, sig, a, ntp)
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
      type(ptcl)::pj0           ! input. projectile.
      integer,intent(in)::ia0    ! target A (nucleon #)
      integer,intent(in)::iz0    ! target Z
      real(8),intent(in)::sig   ! cross-setion (mb) on this target
                         ! current  qmd00  dose not use but may  be
                         ! better to give it
	type(ptcl):: a(*)  ! output. produced particles
	      ! supoose Fe+Bi collsion.  If all Fe 
              ! nucleon of energy 4 GeV collide and produce
              !( 2pions + n  )*56 + 200 n; 3*56+200 = 380; max size of a.

	integer,intent(out):: ntp   ! total number of produced ptcls
	integer n, i, j

	integer code, subcode, charge
	integer ityp,ktyp
	real(8):: KEpn
        real(8):: u
        real(8):: totxs,  elaxs
        integer:: icon
        integer:: inela  ! used to select only inela event for Jam

        type(ptcl)::tgt, pj
        integer ia, iz
        logical,save:: exchange

        real(8),parameter::KEpnB=1.500  ! GeV/n (KE). below this pA or Ap
                                        ! Bertin is usd. 
        real(8),parameter::KEneucpro=1.0  ! n+p, use Bertin below this and Jqmd above this
                !        for  np use jqmd

        KEpn =( pj0%fm%p(4) - pj0%mass )
        if( pj0%code == kgnuc) then
           KEpn= KEpn/pj0%subcode
        endif


        if( pj0%code == kgnuc .and.  ia0 == 1 ) then
!             Ap or An collision. treat it as pA or nA 
           exchange = .true.
!             reverse the Pj and Traget;
!              first,  make tgt
           call cmkptc(knuc, -1, iz0, tgt)
           tgt%fm%p(1:3)= 0.
           tgt%fm%p(4) = tgt%mass
!              see tgt from pj0 and save it in pj
           call cbst1(1, pj0, tgt, pj)
!              make tgt
           ia = pj0%subcode
           iz = pj0%charge
        else
           exchange =.false.
           pj = pj0
           ia = ia0
           iz = iz0
        endif
           

	code = pj%code
	subcode  = pj%subcode
	charge = pj%charge


	call ccos2phits(code, subcode, charge, ityp, ktyp)
        if( code == knuc  ) then
           if( subcode /= antip ) then
              call rndc(u)
              call cphitsXs(pj, ia, iz, elaxs, totxs,icon)
              if(icon == -1) then
                 ntp = 0
              else 
                 if(u < elaxs/totxs .and. .not. exchange ) then
                     ! elastic. for Ap or An case, only inela. 

#if defined (CHECKPHITS)
                    write(0,*) ' entering elastic: pj=',
     *                  pj%code, pj%charge
#endif
                    call cnelas(pj, ia, iz, a, ntp)

#if defined (CHECKPHITS)
                    call cprintptcl(1,'after cnelas')
#endif
                 else

                    if( ( KEpn < KEpnB .and. ia > 1)  .or.
     *                  ( ia == 1 .and.  KEpn < KEneucpro) ) then

#if defined (CHECKPHITS)
                       write(0,'(a, 3i6, 1p,g12.4)')
     *                 ' entering cbertin pj code & KE(MeV)=',pj%code,
     *                  pj%subcode, pj%charge, (pj%fm%p(4)-pj%mass)*1.e3
                       write(0,*) ' target A,Z=', ia, iz
#endif

                       call cbertini(pj, ia, iz, sig, a, ntp)

#if defined (CHECKPHITS)
                       call cprintptcl(2, 'after cbertini 1')
#endif

                    else
                       call cjqmd(pj, ia, iz, sig, a, ntp)

#if defined (CHECKPHITS)
                       call cprintptcl(2, 'after cjqmd 2')
#endif

                    endif
                 endif
              endif
           else

#if defined (CHECKPHITS)
              write(0,'(a, 3i6, 1p,g12.4)')
     *             ' entering cbertin pj code & KE(MeV)=',pj%code,
     *             pj%subcode, pj%charge, (pj%fm%p(4)-pj%mass)*1.e3
              write(0,*) ' target A,Z=', ia, iz
#endif
              call cjqmd(pj, ia, iz, sig, a, ntp)

#if defined (CHECKPHITS)
              call cprintptcl(2, 'after cjqmd 1')
#endif
!////////////////
!           if( pj.code == 6 .and. pj.subcode == antip .and.
!     *        pj.charge == -1 .and. pj.fm.p(4) - pj.mass <= 0.)
!     *         then
!              write(0,*) '1  sig=',sig, ' ntp=',ntp
!              do i = 1, ntp
!                 write(0,*)  a(i).code, a(i).fm.p(4)
!              enddo
!           endif
!//////////////
           endif
        elseif( code == kpion ) then

#if defined (CHECKPHITS)
           write(0,*) ' pion enters cbertini'
           write(0,'(a, i3, 1p,g12.4)') ' chg  KE(MeV)=', pj%charge, 
     *         (pj%fm%p(4)-pj%mass)*1.e3
           write(0,*) ' target,A,Z=', ia, iz
#endif

           call cbertini(pj, ia, iz, sig, a, ntp)

#if defined (CHECKPHITS)
           call cprintptcl(2,'after cbertini 2')
#endif

        else

#if defined (CHECKPHITS)
           write(0,'(a, 3i6, 1p, g12.4)')
     *        ' Entering jqmd, pj & KE(MeV)=',pj%code,
     *         pj%subcode, pj%charge,
     *         (pj%fm%p(4)-pj%mass)*1.e3
           write(0,*) ' target,A,Z=', ia, iz
#endif

           call cjqmd(pj, ia, iz, sig, a, ntp)

#if defined (CHECKPHITS)
           call cprintptcl(2, 'after cjqmd 2')
#endif

!////////////////
!           if( pj.code == 6 .and. pj.subcode == antip .and.
!     *        pj.charge == -1 .and. pj.fm.p(4) - pj.mass <= 0.)
!     *         then
!              write(0,*) '2 sig=',sig, ' ntp=',ntp
!              do i = 1, ntp
!                 write(0,*)  a(i).code, a(i).fm.p(4)
!              enddo
!           endif
!//////////////
        endif
        if( exchange ) then
           !    move to reset system of pj0
           do i = 1, ntp
#if defined (CHECKPHITS)
              write(0,'(a,1p,4g12.4)') 'bef ibst: pj0 rest system', 
     *          a(i)%fm%p(1:3)*1000., 
     *          (a(i)%fm%p(4)-a(i)%mass)*1000.
#endif
              call cibstPol(i, pj0, a(i), a(i))  ! output canbe a(i)

#if defined (CHECKPHITS)
              write(0,'(a,4f7.1)') 'aft ibst;orig sys. ', 
     *          a(i)%fm%p(1:3)*1000., 
     *          (a(i)%fm%p(4)-a(i)%mass)*1000.
#endif

           enddo
        endif

        end subroutine
      subroutine cphitsXs(pj, ia, iz, elaxs, totxs, icon)
      use jqmd
      use bertini
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmass.h"

      type(ptcl):: pj  ! input projectile
      integer,intent(in)::ia,iz  ! target A,Z
      real(8), intent(out):: elaxs
      real(8),intent(out):: totxs   ! xs in mb
                                 ! for p,n total (inela + elaxs)
                                 ! for Heavy  inela, elaxs =
      integer,intent(out):: icon ! =0; oK =-1; N%G

      integer::incp
      real(8),save:: EkMeV, inelaxs
      real(8)::ap,zp,ep,at,zt,bmax

!            used to judge elastic or inela in cphits
      integer,save:: nstrange=0



      if(pj%code == knuc ) then
         EkMeV =( pj%fm%p(4) - pj%mass) *  1000.d0
         if(pj%charge == 1 ) then
            incp = 1
         elseif(pj%charge == 0 .and. pj%subcode == regptcl) then
            incp = 2
         else
            nstrange = nstrange + 1
            if( nstrange < 10) then
               write(0,*) 'in cphits: code, subc, charge=', pj%code,
     *             pj%subcode, pj%charge
               write(0,*) ' strange; neglect  '
            elseif( nstrange > 100 ) then
               write(0,*) ' too many strange ptcl in cphits'
               stop
            endif
            inelaxs = 0.
            elaxs = 0.
            icon = -1
            return
         endif
         call sigrc(incp, EkMeV, ia, iz, totxs, inelaxs, elaxs)
!             abvoe xs is in b.
      elseif(pj%code == kgnuc ) then
         ap = pj%subcode
         zp = pj%charge
         ep = pj%fm%p(4) * 1000.d0   ! total E in MeV
         at = ia
         zt = iz
!                use Shen's (or nasa) cross-section. (icrhi=0 in cprephists
!                  choose Shen's one (if it is  1-->Nasa XS)
!                elaxs = 0 always. bmax not used
         call sighi(ap,zp,ep,at,zt,inelaxs, elaxs,bmax)
         totxs =  inelaxs + elaxs
      else
         nstrange = nstrange + 1
         if( nstrange < 10 ) then
            write(0,*) 'in cphits: code, subc, charge=', pj%code,
     *        pj%subcode, pj%charge
            write(0,*) ' strange; neglected'
         elseif ( nstrange > 100) then
            write(0,*) ' too many strange ptcl in cphits'
            stop
         endif
         inelaxs = 0.
         elaxs = 0.
         icon = -1
      endif
      totxs  = totxs*1000.
      inelaxs = inelaxs *1000.
      elaxs = elaxs *1000.
      end subroutine
      subroutine cphitsADJnp(pj, ia, iz)
!     For pj= neutron and heavy ion target,
!    "grey" nucleons (knock-out nucleons)
!     are mostly neutrons in the phits model; this
!     is somewhat unnatural so, we re-assign the charge of
!     nucleons.
!     In the case of  pj=p, protons are too much, so we also
!     re-assign  the charge.  The total charge conservation is done
!     by adjusting the residual heavy ion charge which is placed
!     in the last of iclust array.
!     This routine should be called before the evaporation
!     routine 'nevap' is called,

!      nclst, iclst in
!      common /clustf/ nclst, iclust(nnn) 
!      may have effect.

!    This routine is activated when DoNPadjust=1
!    (see save statement below)
!    and will be de-activated by putting 0.
!
      implicit none 

#include "Zcode.h"
#include "Zptcl.h"
#include "Zmass.h"
#include "Zkfcode.h"
#include "Zevhnp.h"
      type(ptcl)::pj  ! projectile 
      integer,intent(in):: ia ! target A
      integer,intent(in):: iz ! target Z 
!              next is from Zevhnp.h

      integer nnn, nomp
#include "param00.inc"

      integer nclst, iclust
!        nclst   : total number of out going particles and nuclei    

      common /clustf/ nclst, iclust(nnn)
      integer jclust
      real(8):: qclust
      common /clustg/ jclust(0:7,nnn),  qclust(0:12,nnn)

      integer numpat(0:20)
      real(8):: rumpat(0:20)
      common /clustp/ rumpat, numpat

      integer jclusts
      integer nclsts, iclusts
      real(8):: qclusts
      common /clustt/ nclsts, iclusts(nnn)
      common /clustw/ jclusts(0:7,nnn),  qclusts(0:12,nnn)
!

      integer i, j
      real(8):: probp, u
!!!!!!!!!!!!
!      record /ptcl/ a(100)
!      integer ntp
!!!!!!!!!!!
!      should be  j=  iclust(xxx)
!                j = 0, nucleus
!                  = 1, proton 
!                  = 2, neutron
!                  = 3, pion   
!                  = 4, photon 
!                  = 5, kaon   
!                  = 6, muon   
!                  = 7, others 
!          jclust(j,:)         
!                j = 0, angular momentum  
!                  = 1, proton number   
!                  = 2, neutron number
!                  = 3, ip, see below 
!                  = 4,               
!                  = 5, charge        
!                  = 6, baryon number 
!                  = 7, kf code       

!        qclust(i,nclst)
!                i = 0, impact parameter 
!                  = 1, px (GeV/c)  
!                  = 2, py (GeV/c)  
!                  = 3, pz (GeV/c)  
!                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)
!                  = 5, rest mass (GeV)
!                  = 6, excitation energy (MeV)
!                  = 7, kinetic energy (MeV) 
!                  = 8, weight change 
!                  = 9, delay time 
!                  = 10, x-displace
!                  = 11, y-displace
!                  = 12, z-displace
!                          
!        numpat(i) : total number of out going particles or nuclei
!                i =  0, nuclei         
!                  =  1, proton      
!                  =  2, neutron    
!                  =  3, pi+ 
!                  =  4, pi0 
!                  =  5, pi- 
!                  =  6, mu+ 
!                  =  7, mu- 
!                  =  8, K+  
!                  =  9, K0  
!                  = 10, K-  
!                            
!                  = 11, other particles
!                                   
!                  = 12, electron 
!                  = 13, positron 
!                  = 14, photon  
!
!                  = 15, deuteron
!                  = 16, triton  
!                  = 17, 3He     
!                  = 18, Alpha   
!                  = 19, residual nucleus 
      logical dojob
      integer::sump, sumA
      integer:: HZ, HA  ! residual heavy ion'z Z,A
              !    HA is kept same but HZ may changen
      integer:: HZmin, HZmax   ! min/maximum acceptable residual Z
                 !  when HA is fixed.
      real(8):: Zs
      integer:: charge, subcode
      if( DoNPadjust ==  0 )  return   ! ********** 

#if defined (CHECKPHITS)
      if(sum( jclust(5,1:nclst) ) /= pj%charge + iz ) then
         write(*,*) ' ******** charge non consv '
         write(*,*) sum( jclust(5,1:nclst) ), ' input=', pj%charge+iz
         stop
      endif
#endif

      
      if( pj%code == knuc  .and.  ia > 10 .and. nclst > 0  ) then
!         if( ( (  pj.charge == 0 .and.  iclust(1) == 2 ) .or.
!     *       (  pj.charge == 1 .and.  iclust(1) == 1 ) )  .and.
!     *          qclust(7,1)/1000. > (pj.fm.p(4)-pj.mass)*0.95)  then
         if( qclust(7,1)/1000. > (pj%fm%p(4)-pj%mass)*0.95)  then
            return
         endif

           ! last one should be heavy remnant. require > He
         if( iclust(nclst) == 0 .and. jclust(5,nclst) > 2 ) then
            HA = jclust(6, nclst) ! remnant A; fixed
!            write(*,*) ' pn=',jclust(1,nclst)
!            write(*,*) ' nn=',jclust(2,nclst)
!            write(*,*) ' kf=',jclust(7,nclst)
            call cphitsFixMinMaxZ(HA, Zs, HZmin, HZmax)
            do j = 1, 3
               HZ = jclust(5, nclst) ! remnant Z; may be changed
!/////
!               write(*,*) '***HZ=',HZ
!/////////////
!               sump = jclust(5,1) ! # of p
               sump = 0
               sumA = 0
!               do i = 2, nclst-1 ! first one is leading nucleon, don't  touch
               do i = 1, nclst-1 
                  probp = float(iz-sump)/(ia-sumA) ! prob. of p
                           !  among nucleons  before HA is formed
                  if( iclust(i) == 1 .or. iclust(i) == 2 ) then
                     dojob = .true. ! n or p
                  elseif( iclust(i) == 0 ) then ! heavy ?
                     charge = jclust(5,i)  
                     subcode = jclust(6,i)  
                     if( subcode == 1 .and. charge <=1 ) then
                        dojob = .true. ! some strange ptcl
                                     ! treat it as p/n
                        iclust(i) = -charge+2
                     else       ! real heavy
                        sump = charge + sump
                        sumA = sumA + subcode
                        dojob = .false.
                     endif
                  else
                     ! pi etc do nothing
                     dojob = .false.
                  endif
                  if(dojob) then
                     call rndc(u)
                     if(u < probp) then
                        if( iclust(i) /= 1 ) then
                     !     use ..s as working array
                           iclusts(i) = 1
                           jclusts(5,i) = 1
                           qclusts(5,i) = masp
                           jclusts(1,i) = 1
                           jclusts(2,i) = 0
                           jclusts(7,i) = kfproton
                             ! current  n --> p
                           HZ = HZ - 1
!/////
!               write(*,*) 'i', i, ' - HZ=',HZ
!/////////////
                           sump = sump + 1
                        else
                           iclusts(i) = iclust(i)
                           jclusts(5,i) = jclust(5,i)
                           qclusts(5,i) = qclust(5,i)
                           jclusts(1,i) = jclust(1,i)
                           jclusts(2,i) = jclust(2,i)
                           jclusts(7,i) = jclust(7,i)
                           sump = sump + 1
                        endif
                        sumA = sumA + 1
                     else
                        if(iclust(i) /= 2 ) then
                           iclusts(i) = 2
                           jclusts(5,i) = 0
                           qclusts(5,i) = masn
                           jclusts(1,i) = 0
                           jclusts(2,i) = 1
                           jclusts(7,i) = kfneutron
                                ! current p-->n
                           HZ  = HZ  + 1
!/////
!               write(*,*) 'i ', i, '+ HZ=',HZ
!/////////////
                           sump = sump - 1
                        else
                           iclusts(i) = iclust(i)
                           jclusts(5,i) = jclust(5,i)
                           qclusts(5,i) = qclust(5,i)
                           jclusts(1,i) = jclust(1,i)
                           jclusts(2,i) = jclust(2,i)
                           jclusts(7,i) = jclust(7,i)
                           sump = sump + 1
                        endif
                        sumA = sumA + 1
                     endif
                  else
                       ! pion etc
                     iclusts(i) = iclust(i)
                     jclusts(5,i) = jclust(5,i)
                     qclusts(5,i) = qclust(5,i)
                     jclusts(1,i) = jclust(1,i)
                     jclusts(2,i) = jclust(2,i)
                     jclusts(7,i) = jclust(7,i)
                     sump = sump +  jclust(5,i)
                  endif
               enddo
               if( HZ >= HZmin .and. HZ <= HZmax) then
                     ! OK
                     ! move ...s to without s
!     
!                  write(*,*) 'OK j=',j, ' HZ=',HZ,
!     *             ' min max=',HZmin, HZmax
                 
                  do i = 1, nclst-1
!////////
!                    write(*,*)i,  iclust(i),"-->", iclusts(i)
!/////////
                     iclust(i) = iclusts(i)
                     jclust(5,i) = jclusts(5,i)
                     qclust(5,i) = qclusts(5,i)
                     jclust(1,i) = jclusts(1,i)
                     jclust(2,i) = jclusts(2,i)
                     jclust(7,i) = jclusts(7,i)
                  enddo
                  jclust(1, nclst) = HZ
                  jclust(2, nclst) = HA-HZ
                  jclust(5, nclst) = HZ  ! rest remnant charge
                  jclust(7, nclst) = HZ*1000000+HA ! kfcode ?
!                  call cphitsOut(1, pj, ia, iz, ntp, a)
!                      next will not be used in nevap
!                  numpat(0:2) = 0
!                  do i =1, nclst
!                     if(iclust(i) == 0) numpat(0) = numpat(0)+1
!                     if(iclust(i) == 1) numpat(1) = numpat(1)+1
!                     if(iclust(i) == 2) numpat(2) = numpat(2)+1
!                  enddo
                  exit
               endif
            enddo   ! after some trials, HZ did not satisfy
                    ! so use origianl one (kept in w/o s array) 
         endif
      endif
      end
      subroutine cphitsFixMinMaxZ(HA,  Zs, HZmin, HZmax)
      implicit none
      integer,intent(in):: HA  ! heavy ion A
      real(8),intent(out):: Zs  ! most stable Z for a give HA.
      integer,intent(out):: HZmin  ! minimum Z for HA
      integer,intent(out):: HZmax  ! maximum Z for HA

! From: http://www.th.phys.titech.ac.jp/~muto/lectures/INP02/INP02_chap03.pdf
      real(4),parameter::bsym = 23.29
      real(4),parameter::bcoul = 0.697
!      real(4),parameter::bsurf = 17.23
!      real(4),parameter::bvol = 15.56
      
! http://ja.wikipedia.org/wiki/核種の一覧
! however, we restrict to familier ones for Z=1 and Z=2
! and accept rather moderate range around the most stable
! line 
       !  Weizsacker.  /2  is missing in ~muto
      Zs = HA/(2.0 + bcoul*HA**0.666666/bsym/2)
      select case(HA)

      case(22:220)
         HZmax = Zs +  6.0*HA/140.
         HZmin = Zs -  2.5*HA/140.
      case(1:2)  
         HZmin = 1
         HZmax = 1
         Zs = 1
      case( 3 )
         HZmin = 1
         HZmax = 2
         Zs = 1.99
      case( 4 )
         HZmin = 2
         HZmax = 2
         Zs = 2
      case( 5 )
         HZmin = 3
         HZmax = 4
         Zs = 3.1
      case( 6 )
         HZmin = 3
         HZmax = 5
         Zs = 4.0
      case( 7 )
         HZmin = 3
         HZmax = 5
         Zs =  4.
      case( 8:10)
         HZmin = 4
         HZmax = 6
         Zs =5.
      case( 11:12)
         HZmin = 5
         HZmax = 7
         Zs = 6.
      case( 13:14)
         HZmin = 6
         HZmax = 8
         Zs = 7
      case( 15:17)
         HZmin = 7
         HZmax = 9
      case( 18:19)
         HZmin = 8
         HZmax = 10
      case( 20:21 )
         HZmin = 8
         HZmax = 11
      case default  
         write(0,*) ' strange HA=', HA, ' for cphitsFixMinMaxZ'
         stop
      end select
      end

