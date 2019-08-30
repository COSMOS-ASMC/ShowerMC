      subroutine cphitsADJnp(pj, ia, iz)
!     For pj= neutron and heavy ion target,
!    "grey" nucleons (knock-out nucleons)
!     are mostly neutrons in the phits model; this
!     is somewhat unnatural so, we re-assign the charge of
!     nucleons.
!     In the case of  pj=p, protons are too much, so we also
!     re-assign  the charge.  The total charge conservation is done
!     by adjusting the residual heavy ion charge which is placed
!     in the last of icust array.
!     This routine should be called before the evaporation
!     routine 'nevap' is called,

!      nclst, iclst in
!      common /clustf/ nclst, iclust(nnn) 
!      may have effect.

!    This routine is activated when DoNPadjust=1 (hardwired) 
!    (see save statement below)
!    and will be de-activated by putting 0.
!
      implicit none 

#include "Zcode.h"
#include "Zptcl.h"
#include "Zmass.h"
      type(ptcl)::pj  ! projectile 
      integer,intent(in):: ia ! target A
      integer,intent(in):: iz ! target Z 

      integer,save:: DoNPadjust=1  ! 0--> no adjust
!      include "JQMD/param00.inc"

      integer nclst, iclust
!      integer nclsts, iclusts
!        nclst   : total number of out going particles and nuclei      *
      common /clustf/ nclst, iclust(nnn)
      integer jclust
      real(8):: qclust
      common /clustg/ jclust(0:7,nnn),  qclust(0:12,nnn)
!
      integer jclusts
      real(8):: qclusts
      common /clustt/ nclsts, iclusts(nnn)
      common /clustw/ jclusts(0:7,nnn),  qclusts(0:12,nnn)
!

      integer i, j
      real(8):: probp, u
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
      integer::npnow, nnnow
      integer:: HZ, HA  ! residual heavy ion'z Z,A
              !    HA is kept same but HZ may changen
      integer:: HZmin, HZmax   ! min/maximum acceptable residual Z
                 !  when HA is fixed.
      real(8):: Zs
      integer:: charge, subcode

      if( DoNPadjust ==  0 )  return   ! ********** 
      dojob = .false.
      if( pj.code == knuc  &&  ia > 10 && nclst > 0  ) then
         if( iclust(nclst) == 0 .and. jclust(5,i) > 2 ) then
           ! last one should be heavy remnant.
           ! do below if it is  > He
            HA = jclust(6, nclst) ! remnant A
            call cphitsFixMinMaxZ(HA, HZ, Zs, HZmin, HZmax)
            do j = 1, 3
               HZ = jclust(5, nclst) ! remnant Z; may be changed
               sump = jclust(5,1) ! # of p
               do i = 2, nclst-1 ! first one is leading nucleon, don't  touch
                  probp = float(iz-sump)/ia ! prob. of p
                           !  among nucleons  before HA is formed
                  if( iclust(i) == 1 .or. iclust(i) == 2 ) then
                     dojob = .true. ! n or p
                  elseif( iclust(i) == 0 ) then ! heavy ?
                     charge = jclust(5,i)  
                     subcode = jclust(6,i)  
                     if( subcode == 1 .and. charge <=1 ) then
                        dojob = .true. ! some strange ptcl
                                     ! treat it as p/n
                     else       ! real heavy
                        sump = charge + sump
                     endif
                  endif
                  if(dojob) then
                     call rndc(u)
                     if(u < probp) then
                        if( iclust(i) /= 1 ) then
                     !     use ..s as working array
                           iclusts(i) = 1
                           jclusts(5,i) = 1
                           qclusts(5,i) = massp
                             ! current  n --> p
                           HZ = HZ + 1
                           sump = sump + 1
                        else
                           iclusts(i) = iclust(i)
                           jclusts(5,i) = jclust(5,i)
                           qclusts(5,i) = qclust(5,i)
                        endif
                     else
                        if(icust(i) /= 2 ) then
                           iclusts(i) = 2
                           jclusts(5,i) = 0
                           qclusts(5,i) = massn
                                ! current p-->n
                           HZ  = HZ  - 1
                           sump = sump - 1
                        else
                           iclusts(i) = iclust(i)
                           jclusts(5,i) = jclust(5,i)
                           qclusts(5,i) = qclust(5,i)
                        endif
                     endif
                  endif
               enddo
               if( HZ >= HZmin .and. HZ <= HZmax) then
                     ! OK
                     ! move ...s to without s
                  do i = 2, nclst-1
                     iclust(i) = iclusts(i)
                     jclust(5,i) = jclusts(5,i)
                     qclust(5,i) = qclusts(5,i)
                  enddo
                  jclust(5, nclst) = HZ  ! rest remnant charge
                  exit
               endif
            enddo   ! after some trials, HZ did not satisfy
                    ! so use origianl one (kept in w/o s array) 
         endif
      endif
      end
      subroutine cphitsFixMinMaxZ(HA, HZ, Zs, HZmin, HZmax)
      implicit none
      integer,intent(in):: HA  ! heavy ion A
      integer,intent(in):: HZ  ! heavy ion charge Z; not used now 
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
