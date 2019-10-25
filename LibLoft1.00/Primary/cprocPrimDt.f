!    ********************************************************
      subroutine cprocPrimDt(prm)
!        process primaries( examine  primary data in prm. and
!        make some computation and store the results in prm)
!        character data is coverted into lower case.
!       prm: /primaries/  Input/output
!
      implicit none

#include  "Zmanagerp.h"
#include  "Zptcl.h"
#include  "Zprimary.h"
!
      type(primaries):: prm

      character*70 msg
!
      integer i, icon
      character*16 temp1
!
      icon = 0
      prm%NoOfSamplings = EventsInTheRun ! counter, including discarded ones
                                        ! how many primaries generated
      EventNo = PrevEventNo             ! the same

      do i = 1, prm%no_of_comps
         prm%NoOfSampComp(i, 1) =0  ! all
         prm%NoOfSampComp(i, 2) =0  ! only for accepted ones
!             to lower case string
          temp1 = prm%each(i)%symb
          call c2lowerCase(temp1, prm%each(i)%symb)
          temp1 = prm%each(i)%eunit
          call c2lowerCase(temp1, prm%each(i)%eunit)
          temp1 = prm%each(i)%etype
          call c2lowerCase(temp1, prm%each(i)%etype)
          temp1 = prm%each(i)%diff_or_inte
          call c2lowerCase(temp1,  prm%each(i)%diff_or_inte)
!
          prm%each(i)%label = i           ! numbering
          call cexmPrimSymb(prm%each(i), icon)
          call cexmPrimEU(prm%each(i),icon)
          call cmkPrimSTbl(prm%each(i),icon)
      enddo
      if(icon .ne. 0) then
         write(msg, *) ' correct primary data table'
         call cerrorMsg(msg, 0)
      endif
      prm%cummInteFlux(1) = prm%each(1)%inte_value
      do i = 1, prm%no_of_comps-1
          prm%cummInteFlux(i+1) =   prm%cummInteFlux(i)
     *        + prm%each(i+1)%inte_value
      enddo
!       normalize 
      do i = 1, prm%no_of_comps
          prm%cummInteFlux(i) =  prm%cummInteFlux(i) /
     *          prm%cummInteFlux(prm%no_of_comps)
      enddo
      end
!     ***********************************
      subroutine cexmPrimSymb(each, icon)
!     ***********************************
!           examine a given primary symbol and if it is valid
!       one, set each.code, each.subcode, each.charge
!       if not, icon = 1 is given.
!
      implicit none

#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryc.h"
      type(component):: each
      integer icon
!
!
      character(16):: symb  !  12->16 (2019/Aug/29)
      integer k,  pap, massn, chgn
      character*70 msg
      integer:: kf

      symb = each%symb
      icon = 0
      k = index(symb, '~') 
      if(k .gt. 0) then
!              anti particle specification
            pap = 1
            symb = symb(1:k-1)  !drop ~
      else
            pap = 0
      endif
!      
      k = 1
      if(symb(1:3) .eq.  'iso')  then
!                  read A, Z
         read(symb(4:16), *) massn, chgn   !  12->16 2019/Aug/29
         if(massn .le. chgn) then
            write(msg,*) ' primary =',symb,
     *     ' invalid(becaus A=',massn,'<= Z=',chgn
            call cerrorMsg(msg, 0)
         endif
         symb = 'iso'          
      elseif(symb(1:3) .eq.  'pdg')  then
!                  read PDG code   ! no heavy primary
         read(symb(4:16), *) kf !  12->16 2019/Aug/29  now heavy
                                ! 1ry can be managed by pdg   
         call ckf2cos(kf, each%code, each%subcode, each%charge)
         symb = 'pdg'           ! erase kf code
         return !!!!!!!
      endif


      
      do while ( k .le. NoOfSymbols )
         if(PrimaryIdTbl(k)%symb .eq. symb) then
            each%code = PrimaryIdTbl(k)%code
            each%subcode = PrimaryIdTbl(k)%subcode
            each%charge = PrimaryIdTbl(k)%charge
            k = NoOfSymbols +1
         endif
         k = k+1
      enddo

      if(k .ne. NoOfSymbols + 2) then
         write(msg, *) each%label, '-th primary component=',
     *           each%symb,' invalid'
         call cerrorMsg(msg, 1)
         icon = 1
      else
         if(symb .eq. 'iso') then
            each%subcode = massn
            each%charge = chgn
         endif
         if(pap .ne. 0) then
!               anti particle
            each%charge  = - each%charge
            if( symb .ne. 'iso' ) then
               each%subcode = - each%subcode               
            endif
         endif   
      endif
      end
!     *********************************
      subroutine cexmPrimEU(each, icon)
!     *********************************
!           examine a given primary energy unit and if it is valid
!       one, set each.togev
!       if not, icon = 1 is given.
!
      implicit none
#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryc.h"
      type(component):: each
      integer  icon
!     
!
      character(16):: symb  !  12-> 16  2019/Aug/29
      integer k
      character*70 msg

      symb = each%eunit
!      
      k = 1
      do while ( k .le. maxErgUnit)
         if(ErgUnitTbl(k)%symb .eq. symb) then
             each%togev = ErgUnitTbl(k)%togev
             k = maxErgUnit +1
         endif
         k = k+1
      enddo
      if(k .ne. maxErgUnit + 2) then
          write(msg, *) each%label, '-th primary energy unit symbol=',
     *    each%eunit,' invalid'
          call cerrorMsg(msg, 1)
          icon = 1
      endif
!             set emin and emax
      each%emin = each%energy(1)
      each%emax = each%energy(each%no_of_seg + 1)
      end
