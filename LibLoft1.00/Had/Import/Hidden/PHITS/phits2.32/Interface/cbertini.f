#include "ZcheckPHITS.h"
      subroutine cbertini(pj, ia, iz, sig, a, ntp)

      use modnevap
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
	type(ptcl):: pj  ! input. projectile   
        integer,intent(in):: ia  ! target A 
        integer,intent(in):: iz  ! target Z
        real(8),intent(in):: sig ! cross-section in mb on this target
                                 ! at presnt not used.
        type(ptcl):: a(*) ! output. generated ptlcs.
        integer,intent(out):: ntp   ! # of generated ptcls


	integer n

	integer ityp, ktyp
!/////////
        integer::j
!////////////
        real(8):: eein
	integer code, subcode, charge
!        real(8):: Tm  ! max recoil energy in lab
!        real(8):: T   ! recoil energy  in lab
!        real(8):: tetacm, tetalab ! angle in CM
!        real(8)::  pt

	integer,save::jcasc=1


*       jcasc   : =1 : bertini, =2 : isobert (2 should not be used)


*       ityp    : particle type of projectile                          *
*       eein    : energy of projectile (MeV)                           *
*       mmas    : mass of target                                       *
*       mchg    : charge of target                                     *
	code = pj%code
	subcode  = pj%subcode
	charge = pj%charge
	call ccos2phits(code, subcode, charge, ityp, ktyp)
!	eein = 500*56    total kinetic energy in MeV
	eein = pj%fm%p(4)-pj%mass
	eein = eein *1000.d0    ! MeV


!           sig = pi r**2*10 mb r in fm.
!	bmax0= sqrt(sig/3.141592/10.)


!     Tm = 4.0*mmas/(1.0+mmas)**2 * eein /1000.   ! GeV


        call cbertin(jcasc, ityp, eein, ia, iz)

#if defined (CHECKPHITS)
        call cprintptcl(1,'aft cbertin')
#endif
        call cphitsADJnp(pj, ia, iz)  ! may try to adjust # of n,p
        call nevap(0)

#if defined (CHECKPHITS)
        call cprintptcl(2,'aft cbertin+nevap')
#endif
        call cphitsOut(2,pj, ia, iz,  n, a)   ! n and a are output
	call crot3mom( pj, a, n ) ! rotate to the current  cooord

#if defined (CHECKPHITS)
        write(0,*) ' after crot3mom n=',n
        do j = 1, n
           write(0,'(4f8.1)')
     *         a(j)%fm%p(1:3)*1000., (a(j)%fm%p(4)-a(j)%mass)*1000.
        enddo
#endif

        ntp = n

!        do j=1, n
!           T = a(j).fm.p(4)- a(j).mass
!           tetacm = acos(1.-2*T/Tm)
!           pt = sqrt( a(j).fm.p(1)**2 + a(j).fm.p(2)**2)
!           tetalab = atan2( pt, a(j).fm.p(3) )
!           write(*,'(4i4, 1p, 7g12.4)')
!     *      j, a(j).code, a(j).subcode, a(j).charge, a(j).fm.p(:),
!     *      a(j).mass, tetacm, tetalab
!        enddo
!        write(*,*)
!        

*---- in common -------------------------------------------------------*
*                                                                      *
*        nclst   : total number of out going particles and nuclei      *
*                                                                      *
*        iclust(nclst)                                                 *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclust(i,nclst)                                               *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4,                                                *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                                                                      *
*        qclust(i,nclst)                                               *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, px (GeV/c)                                     *
*                  = 2, py (GeV/c)                                     *
*                  = 3, pz (GeV/c)                                     *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight change                                  *
*                  = 9, delay time                                     *
*                  = 10, x-displace                                    *
*                  = 11, y-displace                                    *
*                  = 12, z-displace                                    *
*                                                                      *
*        numpat(i) : total number of out going particles or nuclei     *
*                                                                      *
*                i =  0, nuclei                                        *
*                  =  1, proton                                        *
*                  =  2, neutron                                       *
*                  =  3, pi+                                           *
*                  =  4, pi0                                           *
*                  =  5, pi-                                           *
*                  =  6, mu+                                           *
*                  =  7, mu-                                           *
*                  =  8, K+                                            *
*                  =  9, K0                                            *
*                  = 10, K-                                            *
*                                                                      *
*                  = 11, other particles                               *
*                                                                      *
*                  = 12, electron                                      *
*                  = 13, positron                                      *
*                  = 14, photon                                        *
*                                                                      *
*                  = 15, deuteron                                      *
*                  = 16, triton                                        *
*                  = 17, 3He                                           *
*                  = 18, Alpha                                         *
*                  = 19, residual nucleus                              *

	end subroutine

      subroutine cbertin(jcasc, ityp, eein, ia, iz)
      use bertini
!ccc	implicit none
#include "param00.inc"
	! interface with jqmdin
      integer,intent(in):: jcasc ! =1 berting
      integer,intent(in):: ityp
      real(8),intent(in)::eein
      integer,intent(in):: ia, iz



      integer ncount

      common /clustf/ nclst, iclust(nnn)

      ncount = 0
      do while (ncount < 100)
         call bertin(jcasc, ityp, eein, ia, iz)
         if( nclst > 0 ) exit
         ncount = ncount + 1
      enddo
!//////
!      write(0,*) 'bertinLoop=', ncount
!//////////
      end


	

	
	



