!----------------------------------------------------
!       cmkptc:  make  a particle
!        
!     implicit none
!     include '../Zptcl.h'
!     include '../Zcode.h'
!     type(ptcl)::  p
!     integer i
!     do i=1, klast
!        call cmkptc(i, 0, 0,  p)
!        write(*, *) p.mass, p.charge
!     enddo
!     end   
      subroutine cmkptc(code, subcode,  charge, p)
!             make a particle. 
!       code: integer. Input. Particle code defined by the Cosmos convention.
!    subcode: integer. Input. Particle subcode defined //
!                      It has meaning for k0. neutron, gamma.
!     charge: integer. Input. Charge of the particle.
!                             In case of heavy (alpha, etd) this should be
!                             1 or -1, indepndently of the real charge.
!                             -1 for anti-neucleus.
!          p: type ptcl. Output.
!                             Template particle is set.
!                  The attributes set are:
!                       px=undef  unchaged
!                       py=   //
!                       pz=   //
!                       e=    //
!                       mass=ptcl mass 
!                       code=ptcl code (same as input)
!                       subcode = ptcl sub code 
!                              This code is mainly used to identify
!                              particle/antiparticle.  If it is not
!                              important, or it is to be determined
!                              later, the user may give 0.
!
!                              This has meaning for the following
!                              particles. For other particles, 
!                              giving 0 is ok. It can be composed by
!                              'code' and 'charge'.
!----------------------------------------------------------------------
!                n           n~         k0s           k0l
! subcode
! defined      kneutron   kneutronb     k0s           k0l   
! in Zcode.h 
!----------------------------------------------------------------------
!           neutrino(e)  neutrino(mu)  neutrino(e)~   neutrino(mu)~
!
! subcode      regptcl            regptcl     antip          antip
!
!----------------------------------------------------------------------
!           direct gamma   brems gamma     d0          d0~
!
! subcode     kdirectg     kcasg          kd0          kdb
!
!----------------------------------------------------------------------
!                       charge=charge (if not heavy neuclus)
!                                     (charge * Z) (charge = 1, 0, -1)
!
!                             If subcode = 0 for  neutral partilces, this
!                             should be reset later, if they are
!                             not symmetric particle (k0, n, d0)
!              
!    
!                 
       implicit none
!----       include '../Zptcl.h'
#include  "Zptcl.h"
!----       include '../Zcode.h'
#include  "Zcode.h"
#include  "Zheavyp.h"
       type(ptcl):: p
       integer code, charge, subcode
!
          p%code = code
!          if(code .ge. kdeut .and. code .le. khvymax) then
          if(code .ge. kalfa .and. code .le. khvymax) then
              call cshvc(code, charge, p)
          elseif(code .eq. kdeuteron ) then
             p%code = kgnuc
             call cshvc(code, charge, p)
          elseif(code .eq. ktriton) then
             p%code = kgnuc
              call cshvc(code, charge, p)
          else
              p%charge = charge
          endif
          call csmass(code, subcode, charge, p)
          call cssubc(code, subcode, charge, p)
!           for heavy, we use only kgnuc here after (from v6.0)
          if(code .ge. kalfa .and. code .le. khvymax) then
             p%subcode = Code2massN(code)
             p%code = kgnuc
          endif   
      end
!     *******************************************************
      subroutine csmass(code, subcode, charge, p)
!          set  particle mass from ptcl code and charge.
!            code: Integer. Input. partcle code defined in COSMOS
!           charge:Integer. Input. partcle charge.
!                p:/ptcl/  Output.  p.mass will get partcle mass in GeV.
!                           For heavy neucleus, (massp + massn)/2*A
!                           is used.
      implicit none
!----      include '../Zptcl.h'
#include  "Zptcl.h"
!----      include '../Zcode.h'
#include  "Zcode.h"
!----      include '../Zmass.h'
#include  "Zmass.h"
!
       integer code, charge, subcode
       type(ptcl):: p
!
       real*8 x
       parameter (x = 1.d50)
       real*8 mass(0:klast, -1:1)
       character*8 id
       integer massn
       character*70  msg
       data 
     * mass(kphoton, :)/x, 0., x/,
     * mass(kelec,:)/masele,x, masele/,
     * mass(kmuon, :)/masmu,x,masmu/, 
     * mass(kpion, :)/maspic, maspi0, maspic/
       data
     * mass(kkaon, :)/maskc, mask0, maskc/,
     * mass(knuc, :)/masp, masn, masp/,
     * mass(kneue, :)/x, 0., x/,
     * mass(kneumu,:)/x, 0., x/,
     * mass(kneumu,:)/x, 0., x/,
     * mass(kneutau,:)/x, 0., x/,       
     * mass(knnb,  :)/x, masnnb, x/
       data
     * mass(kddb, :)/x,masddb,x/,
     * mass(kdmes,:)/masd, masd, masd/
     * mass(krho, :)/masrho,  masrho, masrho/,
     * mass(komega,:)/x,masomg, x/ 
     * mass(kphi,:)/x, masphi, x/
     * mass(keta, :)/x, maseta, x/
     * mass(ketap, :)/x, masetap, x/

       data
     * mass(ksigma, :)/massigmam, massigma0, massigmap/,
     * mass(kgzai,  :)/masgzaim, masgzai0, masgzaim/,
     * mass(klambda,:)/x, maslambda,x/,
     * mass(klambdac,:)/maslambdac, x, maslambdac/,
     * mass(krare, :)/0., 0., 0./,
     * mass(kgnuc, :)/x, x, x/
       data 
     * mass(kbomega, :)/masbomega, x, masbomega/,
     *      mass(kds,:) /masds, x, masds/,
     *      mass(kXic,:)/masXic, masXic0, masXic/,
     *      mass(ktau,:)/mastau, x, mastau/,        
     *      mass(komeC0,:)/x, masomC0, x/
     *      mass(kDelta,:)/masDelta,masDelta, masDelta/ 
!       if(code .ge. kdeut .and. code .le. khvymax) then
       if(code .ge. kalfa .and. code .le. khvymax) then
!                  get mass number
          call cghvm(code, massn)
          p%mass =( masn + masp)  * massn /2
       elseif(code .eq. kdeuteron) then
          p%mass = 1.875613d0
       elseif(code .eq. ktriton) then
          p%mass = 2.80891
       elseif(code  .eq. kgnuc) then
!             general nucleaus (A>1). subcode is A. very rough
!             binding energy. (Weizsacker-Bethe)
          p%mass = masn*(subcode-charge) + masp*charge
     *            -(15.68d-3*subcode-18.56d-3*(float(subcode))**0.6666
     *          -0.717d-3 * charge**2/(float(subcode))**0.33333)
       elseif(code .ge. 0 .and. code .le. klast) then
          p%mass = mass(code, charge)
          if(p%mass .eq. x) then
             call cgpid(code, id)
             write(0, *) ' error detected in csmass in cmkptc:'
             write(0, *)
     *            'input code, subcode, charge =',code, charge, subcode
             write(msg, *)
     *            ' charge=',charge,' invalid for csmass; code=',id
             call cerrorMsg(msg, 0)
          endif
       elseif( code .eq. klight ) then
          p%mass = 0.
       elseif( code .eq. kEdepo .or. code .eq. kchgPath ) then
!         energy deposit or charged ptcl streight path for
!         light emission; nothing to do
       else
            write(msg, *) ' code=',code,' invalid to csmass'
            call cerrorMsg(msg, 0)
       endif
      end
!     *******************************************************
      subroutine cssubc(code, subcode, charge, p)
!            set particle or anti particle subcode from 
!            ptcl code and charge.
!            code: Integer. Input. particle code defined in COSMOS
!          subcode: Integer. Input. paricle sub code //
!          charge:Integer. Input. partcle charge.
!             p: /ptcl/. Output. for most of particles,
!                        'ptcl' or 'antip' is set according to
!                        code and charge. For neutron, k0, gamma
!                        they are treated specially.
!                        for self conjugate particles, 0 is set.
!
      implicit none
!----      include '../Zptcl.h'
#include  "Zptcl.h"
!----      include '../Zcode.h'
#include  "Zcode.h"
#include  "Zheavyp.h"
!
       integer code, subcode, charge
       type(ptcl):: p
       character*70  msg
!
       if(code .ge. 1 .and. code .le. klast) then
!                   this should be consistent with regptcl/antip
!                   def. in Zcode.h
          if(code .eq. kphoton) then
             p%subcode = subcode
          elseif(code .eq. kelec .or. code .eq. kmuon ) then
             p%subcode = - charge * regptcl
          elseif(code .eq. kpion .or. code .eq. kkaon
     *            .or. code .eq. knuc) then
             p%subcode =  charge * regptcl
             if( code .eq. kkaon .and. charge .eq. 0 .and.
     *            subcode .ne. 0) then
                if(abs(subcode) .eq. k0s .or. 
     *               abs(subcode) .eq. k0l ) then
                   p%subcode = subcode
                else
                   write(msg,*) '1 strange subcode=', 
     *                  subcode,' to cssubc. code=', code
                   p%mass = -1.0
                   p%mass = sqrt(p%mass)
                   call cerrorMsg(msg, 0)
                endif
             elseif(code .eq. knuc .and. charge .eq. 0 
     *               .and.   subcode .ne. 0) then
                if(subcode .eq. kneutron .or.
     *               subcode .eq. kneutronb) then
                   p%subcode = subcode
                else
                   write(msg, *) '2 strange subcode=', 
     *                  subcode, ' to cssubc. code=', code
                   call cerrorMsg(msg, 0)
                endif
             endif
          elseif(code .eq. kdmes) then
             if(subcode .ne. 0 .and. charge .eq. 0)then
                if(subcode .eq. kd0 .or.
     *               subcode .eq. kd0b) then                       
                   p%subcode = subcode
                endif
             else
                p%subcode = charge * regptcl
             endif
!          elseif(code .ge. kdeut .and. code .le. khvymax) then
          elseif(code .ge. kalfa .and. code .le. khvymax) then
!             p.subcode = isign(1, charge) *regptcl; set A
             p%subcode = Code2massN(code)   ! mass #
          elseif(code == kdeuteron )  then
             p%subcode =  2
          elseif(code .eq. ktriton ) then
!             p.subcode = isign(1, charge) *regptcl
             p%subcode = 3   !   mass #
          elseif(code .eq. kgnuc) then
             p%subcode = subcode    ! mass #
          elseif(code .eq. kneumu .or. code .eq. kneue) then
             if(subcode .eq. regptcl .or.
     *            subcode .eq. antip .or.
     *            subcode .eq. 0  ) then
                p%subcode = subcode
             else
                write(msg, *) ' 3 strange subcode=', 
     *               subcode, ' to cssubc. code=', code
                call cerrorMsg(msg,  0)
             endif   
          elseif(code .ge. klambda .and.
     *            code .le. klast ) then
             p%subcode = subcode
          else      
             p%subcode = 0      ! should be fixed later
          endif     
       elseif( code .eq. klight .or. code .eq. kEdepo .or.
     *         code .eq. kchgPath ) then
!                not certain. 
          p%subcode = subcode
       elseif(code .eq. krare) then
          p%subcode = 0
       else     
          write(msg, *) ' code=',code,' invalid to cssubc'
          call cerrorMsg(msg, 0)
       endif
      end
!     ****************************************************
!           set heavy neucleus charge
      subroutine cshvc(code, charge, p)
!           code: Integer. Input.  ptcl code
!         charge: Integer. Input.  ptcl charge (1 or -1)
!                                  indicating only positive or
!                                  negative. True charge is
!                                  set here.
!              p: /ptcl/. Output. heavy neucleus charge 
!                           is set in p.charge
!
         implicit none
!----         include '../Zptcl.h'
#include  "Zptcl.h"
!----         include '../Zcode.h'
#include  "Zcode.h"
         integer code, charge
         type(ptcl):: p
         character*70  msg
!
!         integer zhvy(kdeut:khvymax)/1, 2, 4, 7, 12, 17, 26/
         integer zhvy(kalfa:khvymax)/2, 4, 7, 12, 17, 26/
!
!         if(code .ge. kdeut .and. code .le. khvymax ) then
         if(code .ge. kalfa .and. code .le. khvymax ) then
            p%charge =  zhvy(code) * isign(1, charge)
         elseif(code .eq. kdeuteron) then
            p%charge = 1
         elseif(code .eq. ktriton) then
            p%charge = 1
         else
            write(msg, *) 'error input code=',code,' to cshvc'
            call cerrorMsg(msg, 0)
         endif
       end
!     ***************************************************
!         get heavy neucleus mass number
       subroutine cghvm(code, massn)      
!         code: Integer input. ptcl code
!        massn: Integer  output.  mass number
         implicit none
!----         include '../Zcode.h'
#include  "Zcode.h"
#include  "Zheavyp.h"
         integer code, massn
         character*70  msg
!
!
!         if(code .ge. kdeut .and. code .le. khvymax) then
         if(code .ge. kalfa .and. code .le. khvymax) then
            massn = Code2massN(code)
         else
            write(msg, *) 'error input code=',code,' to cghvm'
            call cerrorMsg(msg, 0)
         endif
       end
!     ****************************************************
!           get particle id 
      subroutine cgpid(code, id)
!           get partilce id in character
!        code: Integer. Input.  particle code defined in COSMOS          
!          id: Character*8. Output. partcle id
         implicit none
!----         include '../Zcode.h'
#include  "Zcode.h"
         integer code
         character*8 id
!
         character*70  msg
         character*8 ida(klast)
         data ida(kphoton)/'photon'/, ida(keta)/'Eta'/,
     *        ida(kelec)/'Electron'/, ida(kmuon)/'Muon'/,
     *        ida(kpion)/'Pion'/,     ida(kkaon)/'Kaon'/,
     *        ida(knuc)/'Nucleon'/,   ida(kneue)/'Nue_e'/,
     *        ida(kneumu)/'Nue_mu'/,  ida(knnb)/'NN~'/,
     *        ida(kddb)/'DD~'/,        ida(kdmes)/'D_meson'/,
     *        ida(krho)/'Rho'/,       ida(komega)/'omega'/,
     *        ida(kphi)/'Phi'/,  ida(kgnuc)/'Nucleus'/,
     *        ida(kdeuteron) /'d'/, ida(ktriton)/'t'/
!     *        ida(kphi)/'Phi'/,  ida(kdeut)/'deuteron'/
!                heavy neucleus
         data ida(kalfa)/'Helium'/, ida(klibe)/'LiBeB'/,
     *        ida(kcno)/'CNO'/, ida(khvy)/'NaMgSi'/,
     *        ida(kvhvy)/'SClAr'/, ida(kiron)/'Fe'/,
     *        ida(keta+1)/'light'/, ida(keta+2)/'dE'/,
     *        ida(keta+2)/'cpath'/
         data ida(ksigma)/'sigma'/, ida(klambda)/'lambda'/,
     *   ida(kgzai)/'gzai'/, ida(klambdac)/'lambdac'/,
     *   ida(kbomega)/'Omega'/
         data ida(ktau)/'tau'/, ida(kneutau)/'neu_tau'/
         data ida(kds)/'Ds'/, ida(kXic)/'Xic'/
         data ida(komeC0)/'OmegaC0'/ 
         if(code .ge. 1 .and. code .le. klast)then
              id = ida(code)
         else
              write(msg, *) ' code=',code,' invalid to cgpid'
              call  cerrorMsg(msg,  0)
         endif
      end
!        ------------------------------------------
      subroutine cprptc(p, n)
!           print /ptcl/ strucuture; debug purpose
!      
!----      include '../Zptcl.h'
#include  "Zptcl.h"
      type(ptcl):: p(n)
!
      integer i, j, code
      character*8 id
      character*80 msg

!
      do i=1, n
         code = p(i)%code
         call cgpid(code, id)
         write(msg, *) ' ---------code=',p(i)%code, ' id=', id
         call cerrorMsg(msg, 1)
         write(0, *) ' 4 momentum=',(p(i)%fm%p(j),j=1, 4), ' mass=',
     *               p(i)%mass
!         call cerrorMsg(msg, 1)
         write(msg, *) ' charge=', p(i)%charge, ' subcode=',
     *    p(i)%subcode
         call cerrorMsg(msg, 1)
      enddo   
      end





