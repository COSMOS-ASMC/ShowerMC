!         hadron Air collision for inclusive dpmjet3 
        subroutine cincdpm3(pj, ia, iz, a, np)
        implicit none
#include  "Zmanagerp.h"
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
!c
        type (ptcl):: pj, a(*)
        integer np, ia, iz
        integer maxinclusive
        integer chg, code, subcode, nchild, hcode
        parameter (maxinclusive=500)

        integer kchild(maxinclusive), nextc(maxinclusive)
        real*8  pchild(5, maxinclusive), pin(5)
        integer i
        logical doinc
!
        code = pj%code
        chg = pj%charge
        subcode = pj%subcode
        if(code .eq. kgnuc) then
!           heavy so use original version
!            should not come here; because   for  inclusive treatment
!            superposition model has been already taken
!
!           call  cdpmjet(pj, ia, iz, a, np)
           write(0,*)
     *     'for incdpm3, heavy ptcl should not come here cincdpm3'
           stop 99999
        else
           doinc =  pj%fm%p(4) - pj%mass .gt. 0.2
           if(doinc ) then
              call cccode2hcode(pj, hcode)
              pin(1) = 0.
              pin(2) = 0.
              pin(3) = sqrt(pj%fm%p(4) - pj%mass**2)
              pin(4) = pj%fm%p(4)
              pin(5) = pj%mass
              call hadronint(hcode, 
     *          pin, nchild, kchild, pchild, nextc)

              if(nchild .gt. maxinclusive) then
                 call cerrorMsg(
     *           '# of ptcls by inclusive prod. exceeded limit', 1)
                 call cerrorMsg(
     *          'enlarge maxinclusive in cincdpm3.f', 0)
              endif
              do i = 1, nchild
                 call chcode2ccode(kchild(i), code, subcode, chg)
                 call cmkptc(code, subcode, chg,  a(i))
                 a(i)%fm%p(1) = pchild(1,i)
                 a(i)%fm%p(2) = pchild(2,i)
                 a(i)%fm%p(3) = pchild(3,i)
                 a(i)%fm%p(4) = pchild(4,i)
              enddo
              call crot3mom(pj, a, nchild)
              np = nchild
           else
              ActiveMdl ='byenergy'
!                 use fritiof 1.6 or nucrin
              call chALund(pj, ia, iz, a, np)
           endif
        endif
       end
!      ******************
       subroutine ciniincdpm3
       implicit none
#include  "Zevhnp.h"
!      *********************** for Inclusive init.
       real*8  fact_cross, fac_alpha
       common /crosssection_factor/fact_cross
       common /he_factor/fac_alpha
!      **********************
        character*160  incout
        logical first/.true./
        save    first


        if(first) then
!            read inclusive data table.
           incout = 'InclusiveInfo'
           call initvarip(InclusiveFile, incout)
           first = .false.
        endif
        fact_cross = 0.95
        fac_alpha = 1.77         ! not used in Cosmos
!                                                         
       end
!

      subroutine cccode2hcode(p,  hcode)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      type (ptcl):: p
!           convert cosmos code to honda code
      integer code, subcode, chg  ! input
      integer hcode      ! output. honda code
      character*80  msg

      code = p%code
      chg = p%charge
      subcode = p%subcode

      if(code .eq. knuc) then
         if(chg .eq. 1) then
            hcode = 10
         elseif(chg .eq. -1) then
            hcode = 11
         elseif(chg .eq. 0) then
            if(subcode .eq. antip) then
               hcode = 9
            else
               hcode = 8
            endif
         endif
      elseif(code .eq. kpion) then
         if( chg .eq. 1 ) then
            hcode = 12
         elseif(chg .eq. -1) then
            hcode = 13
         else
            hcode = 14
         endif
      elseif( code .eq. kkaon ) then
         if( chg .eq. 1 ) then
            hcode = 4
         elseif( chg .eq. -1) then
            hcode = 5
         else
            if( subcode .eq. k0l ) then
               hcode = 6
            else
               hcode = 7
            endif
         endif
      else
         call cerrorMsg('strange ptcl input to incdpm3', 1)
         write(msg, '("code=",i3)') code
         call cerrorMsg(msg, 0)
      endif
      end
      subroutine chcode2ccode(hcode, code, subcode, chg)
      implicit none
#include "Zcode.h"
!          hcode to code
      integer  hcode ! input
      integer  code, subcode, chg
      character*80  msg


      if(hcode .eq. 12) then
         code = kpion
         chg = 1
         subcode = antip
      elseif(hcode .eq. 13) then
         code = kpion
         chg = -1
         subcode = regptcl
      elseif(hcode .eq. 14) then
         code = kpion
         chg = 0
         subcode =regptcl
      elseif( hcode .eq. 10) then
         code = knuc
         chg = 1
         subcode = regptcl
      elseif(hcode .eq. 11 ) then
         code = knuc
         chg = -1
         subcode = antip
      elseif(hcode .eq. 8 ) then
         code = knuc
         chg = 0
         subcode = regptcl
      elseif(hcode .eq. 9 ) then
         code = knuc
         chg = 0
         subcode = antip
      elseif(hcode .eq. 4) then
         code = kkaon
         chg = 1
         subcode = antip
      elseif(hcode .eq. 5) then
         code = kkaon
         chg = -1
         subcode = regptcl
      elseif(hcode .eq. 6) then
         code = kkaon
         chg = 0
         subcode = k0l
      elseif(hcode .eq. 7) then
         code = kkaon
         chg = 0
         subcode = k0s
      elseif(hcode .eq. 3) then
         code = kphoton
         chg = 0
         subcode = regptcl
      else
         write(msg,
     *  '("ptcl code =", i3," from incdpm3 is strange")')  hcode
         call cerrorMsg(msg, 0)
      endif
         
      end
