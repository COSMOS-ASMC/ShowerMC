      subroutine subtestBert(pj, mmas, mchg, sig, a, ntp)
	use bertini
	implicit none
#include "Zptcl.h"
	type(ptcl):: pj, a(1000)
        integer ia, iz, ntp
	integer n, i, j
        real(8):: sig

	integer jcasc, ityp, mmas, mchg
	integer iargc
	integer narg, ldat
	character(20)::in
        real(8):: eein
        real(8):: Tm  ! max recoil energy in lab
        real(8):: T   ! recoil energy  in lab
        real(8):: tetacm, tetalab ! angle in CM
        real(8)::  pt
!       jcasc   : =1 : bertini, =2 : isobert, =3 : JAM                 *
!       ityp    : particle type of projectile                          *
!       eein    : energy of projectile (MeV)                           *
!       mmas    : mass of target                                       *
!       mchg    : charge of target                                     *


        
	   write(0,*)
     *    "jcasc  ityp  Ekproj(MeV)  (A, Z)target"
	   write(0,*) 'jcasc: 1--> Bertini model'
	   write(0,*) '       2--> Isobar model'
	   write(0,*) 'ityp: 1-->p, 2-->n,3-->pi+ 4-->pi0 5->pi-'
!           open(33,file="input") 
!           read(33,*) jcasc, ityp, eein,  mmas, mchg
           jcasc = 1
           ityp = 2
           eein = (pj.fm.p(4)-pj.mass)*1000.
           

        Tm = 4.0*mmas/(1.0+mmas)**2 * eein /1000.   ! GeV

!	call cprePhits
!///////////
!        call myPhitsDummy
!////////////
!        do i = 1, 100000
!        do i = 1, 30
           call bertin(jcasc,ityp,eein,mmas,mchg)
           call cphitsOut(n, a)
           if(n <= 0 ) cycle  ! n=-1 if probable. why ?


           do j=1, n
              T = a(j).fm.p(4)- a(j).mass
              tetacm = acos(1.-2*T/Tm)
              pt = sqrt( a(j).fm.p(1)**2 + a(j).fm.p(2)**2)
              tetalab = atan2( pt, a(j).fm.p(3) )
              write(*,'(4i4, 1p, 7g12.4)')
     *         j, a(j).code, a(j).subcode, a(j).charge, a(j).fm.p(:),
     *         a(j).mass, tetacm, tetalab
           enddo
           write(*,*)
 c          eein = eein*0.99
!        enddo

!---- in common -------------------------------------------------------*
!                                                                      *
!        nclst   : total number of out going particles and nuclei      *
!                                                                      *
!        iclust(nclst)                                                 *
!                                                                      *
!                i = 0, nucleus                                        *
!                  = 1, proton                                         *
!                  = 2, neutron                                        *
!                  = 3, pion                                           *
!                  = 4, photon                                         *
!                  = 5, kaon                                           *
!                  = 6, muon                                           *
!                  = 7, others                                         *
!                                                                      *
!        jclust(i,nclst)                                               *
!                                                                      *
!                i = 0, angular momentum                               *
!                  = 1, proton number                                  *
!                  = 2, neutron number                                 *
!                  = 3, ip, see below                                  *
!                  = 4,                                                *
!                  = 5, charge                                         *
!                  = 6, baryon number                                  *
!                  = 7, kf code                                        *
!                                                                      *
!        qclust(i,nclst)                                               *
!                                                                      *
!                i = 0, impact parameter                               *
!                  = 1, px (GeV/c)                                     *
!                  = 2, py (GeV/c)                                     *
!                  = 3, pz (GeV/c)                                     *
!                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
!                  = 5, rest mass (GeV)                                *
!                  = 6, excitation energy (MeV)                        *
!                  = 7, kinetic energy (MeV)                           *
!                  = 8, weight change                                  *
!                  = 9, delay time                                     *
!                  = 10, x-displace                                    *
!                  = 11, y-displace                                    *
!                  = 12, z-displace                                    *
!                                                                      *
!        numpat(i) : total number of out going particles or nuclei     *
!                                                                      *
!                i =  0, nuclei                                        *
!                  =  1, proton                                        *
!                  =  2, neutron                                       *
!                  =  3, pi+                                           *
!                  =  4, pi0                                           *
!                  =  5, pi-                                           *
!                  =  6, mu+                                           *
!                  =  7, mu-                                           *
!                  =  8, K+                                            *
!                  =  9, K0                                            *
!                  = 10, K-                                            *
!                                                                      *
!                  = 11, other particles                               *
!                                                                      *
!                  = 12, electron                                      *
!                  = 13, positron                                      *
!                  = 14, photon                                        *
!                                                                      *
!                  = 15, deuteron                                      *
!                  = 16, triton                                        *
!                  = 17, 3He                                           *
!                  = 18, Alpha                                         *
!                  = 19, residual nucleus                              *

	end subroutine

