      subroutine cphitsOut(id, pj, ia, iz,  n, a)
!       **** note  param00.inc by phits dose not accept implicit none
!            a-h,o-z: double.  i-m; integer
!            New variables here must be explicitly declared
!            without exception 
      implicit none

#include "Zcode.h"
#include "Zptcl.h"
      integer:: nnn, nomp
      include "JQMD/param00.inc"
      type(ptcl):: pj  ! input. projectile 
      integer,intent(in):: ia ! target A
      integer,intent(in):: iz ! target Z 
      integer,intent(in):: id ! 1--> output before nevap; 2-->output after nevap
      integer,intent(out)::n  ! total number of ptcls 
      type(ptcl):: a(*)       !   a(n). each ptcl property

      integer nclst, iclust
      integer nclsts, iclusts
*        nclst   : total number of out going particles and nuclei      *
      common /clustf/ nclst, iclust(nnn)
      integer jclust
      real(8):: qclust
      common /clustg/ jclust(0:7,nnn),  qclust(0:12,nnn)

      integer jclusts
      real(8):: qclusts
      common /clustt/ nclsts, iclusts(nnn)
      common /clustw/ jclusts(0:7,nnn),  qclusts(0:12,nnn)


      integer i, j
*      should be  j=  iclust(xxx)
*                j = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *

*          jclust(j,:)                                                 *
*                j = 0, angular momentum                               *
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
*                                                                      *
*                                                                      *
************************************************************************
      if( id == 1 ) then
         call cphitsOut0( nclst, iclust, jclust, qclust, n, a)
      elseif( id == 2 ) then
         call cphitsOut0( nclsts, iclusts, jclusts, qclusts, n, a)
      else
         write(0,*) 'strange id=',id, ' to cphitsOut'
         stop
      endif      
      end

      subroutine cphitsOut0( ncl, icl, jcl, qcl, n, a)
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
      integer nnn, nomp
      include "JQMD/param00.inc"
      integer,intent(in)::ncl  ! nclst or nclsts:total number of out
                               ! going particles and nuclei  
      integer,intent(in)::icl(nnn) ! iclust or iclusts
      integer,intent(in)::jcl(0:7, nnn)  ! jclst or jclsts
      real(8),intent(in):: qcl(0:12, nnn) ! qclust or qclusts
      integer,intent(out)::n  ! total number of ptcls put in a 
      type(ptcl):: a(*)       !   a(n). each ptcl property
      integer code, subcode, charge

      integer,save:: nstrange= 0
      integer kfcode, i

      n = 0
      do i=1, ncl
         if( icl(i) == 0 ) then
            ! nucleus  but sometimes  p/n comes out  or even
            !  strange ptcl  proton =2  neutron = -1 !
            charge = jcl(5,i)
            subcode = jcl(6,i)
            if( subcode == 1 ) then
               ! should be n or p; if |charge|>1, will be detected later
               code = knuc
               subcode = -1
            else 
               code = kgnuc
            endif
         else
            kfcode=jcl(7,i)
            call ckf2cos(kfcode, code, subcode, charge)
         endif
         if( code /= krare ) then
            if( abs(charge) > 1 .and. code == 6) then
            !  seems to appear for pi- proj.
            !      chg=2. i is normally
            !      last one. chg=3,4 is also seen.
               if( abs( qcl(5,i)- 0.9385d0) < 1.d-3 ) then
                        !  mass
                        !          T.E  - mass      - K.E (MeV) 
                  if( abs( qcl(4,i)-qcl(5,i)-qcl(7,i)/1.d3)
     *                  < 1.d-3 ) then
            !             charge >= 1 code =6 seems proton (since mass is p)
            !            and energy is mass + K.E so we assign it to p
                     charge = 1
                     goto 333
                  endif
               endif
               nstrange = nstrange + 1
               if(nstrange < 10 .or. nstrange >= 100 ) then
                  write(0,*) 'phits: strange ptcl'
                  write(0,*) 'code,sub,charge=', code, subcode,
     *                 charge
                  write(0,*)
     *                'i=',i, ' qcl(1~7,i)=', qcl(1:7, i)
                  write(0,*) ' nstrange=', nstrange
                  if( nstrange >= 100) then
                     write(0,*) 'too many '
                     stop
                  endif
               endif
               cycle  ! neglect this ptcl
            endif
 333        continue
            call cmkptc(code, subcode, charge, a(i))
            n = n + 1
            a(n)%fm%p(4) = qcl(7, i)/1000.d0 + a(n)%mass
                  ! (7,i) is  K.E in MeV; we don't use qcl(4,),qcl(5,)
                  !    since some diff. exists in mass definition
            a(n)%fm%p(1:3) = qcl(1:3, i) ! in GeV/c
!            a(n).fm.p(2) = qcl(2, i)
!            a(n).fm.p(3) = qcl(3, i)
         endif
      enddo
      end subroutine
!///////////////
      subroutine cprintptcl(id, cid)
      
#include "Zcode.h"
#include "Zptcl.h"
      integer:: nnn, nomp
      include "JQMD/param00.inc"
      integer,intent(in):: id ! 1--> output before nevap; 2-->output after nevap
      character(*),intent(in):: cid

      integer nclst, iclust
      integer nclsts, iclusts
*        nclst   : total number of out going particles and nuclei      *
      common /clustf/ nclst, iclust(nnn)
      integer jclust
      real(8):: qclust
      common /clustg/ jclust(0:7,nnn),  qclust(0:12,nnn)

      integer jclusts
      real(8):: qclusts
      common /clustt/ nclsts, iclusts(nnn)
      common /clustw/ jclusts(0:7,nnn),  qclusts(0:12,nnn)
      if(id == 1) then
         call cprintptcl0(cid, nclst, iclust, jclust, qclust)
      else
         call cprintptcl0(cid, nclsts, iclusts, jclusts, qclusts)
      endif
      end subroutine

      subroutine cprintptcl0(cid, ncl, icl, jcl, qcl)
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
      integer nnn, nomp
      include "JQMD/param00.inc"
      character(*),intent(in)::cid

      integer,intent(in)::ncl  ! nclst or nclsts:total number of out
                               ! going particles and nuclei  
      integer,intent(in)::icl(nnn) ! iclust or iclusts
      integer,intent(in)::jcl(0:7, nnn)  ! jclst or jclsts
      real(8),intent(in):: qcl(0:12, nnn) ! qclust or qclusts


      integer code, subcode, charge, kfcode
      integer i

      write(0,*) ' length  of cid=',len(cid), ' ncl=',ncl

      write(0,'(a,a)') ' id=', cid

      write(0,*)
     * '# code subc chg    Px       Py      Pz       KE (MeV)' 
      do i = 1, ncl
         kfcode=jcl(7,i)
!//////////////
         if(kfcode == 2) then
            if( jcl(1,i) == 0  ) then
               !  not proton
               if( jcl(5,i) == 0 ) then
                  ! no charge
                  if(jcl(6,i) == 2 ) then
                     !baryon #; 2 neutrons system
                     !  this will be decayed into two n. in nevap
                     !  so leave it
                     write(0,*) 'two neutron sytem'
                     write(0,*) ' will decay ==>2 n'
                     cycle
                  endif
               endif
            endif
            write(0,*) ' ncl=',ncl
            write(0,*) 'pr0 jcl=',jcl(:,i)
            stop
         endif
!////////////
         if( kfcode > 100000 ) then
            call cphits2cos(15, kfcode, code, subcode, charge)
         else
            call ckf2cos(kfcode, code, subcode, charge)
         endif
         write(0,'(i3, 2i4, 4g13.4)')
     *    code, subcode, charge, qcl(1:3,i)*1.d3, qcl(7,i)
      enddo
      end

      
      

