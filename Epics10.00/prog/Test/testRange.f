#include "ZcosmosBD.h"
#include "ZepicsBD.h"
      program main
      use modepRange
      implicit none
#include "ZepTrackv.h"
#include "Zmass.h"
#include "Zcode.h"
#include "Zcnfig.h"
!
!          test  range
! 
       type(epmedia):: mediax
       type(ptcl)::  aPtcl
      integer  ncomp      
      integer io, subcode, charge, code, nr
      real*8 E,  dedx, w, betag, Ek, PperN, dedxfull
      character*20 epifile, config
      character(256):: input 
      character(20):: b(6)
      integer icon, i
      real(8):: glen, Rcm, Ekmax
      io = 10



      read(*,'(a)') input
      write(0,*) ' input =', trim(input)
      call ksplit(input, 20, 5, b,  nr)
      if( nr /= 5 ) then
         write(0,*) ' wrong number of input = ', nr
         stop
      endif
      
      epifile = b(1)
      config = b(2)
    
      write(0,*) 'epicsfile path is ', trim(epifile)
      write(0,*) 'config path  is ', trim(config)

      call epprmr(epifile)
      call epcmp1  ! compute some (EminGsave etc)
      call eprcnf(config)

      read(b(3), *)  code
      read(b(4), *)  subcode
      read(b(5), *)  charge
      write(0,*) ' code, subcode, charge=',code, subcode, charge

     

      call epqncp(ncomp)

!           make incident
      call cmkptc(code, subcode, charge, aPtcl)


!      call epRangeAlloc(noOfmedia)   ! this has been done in eprcnf
!
!      do i=1, NoOfMedia
!         call epRangeMkTbl(Media(i), i)
!      enddo


      do i = 1, NoOfMedia
         write(*,'(a)') 
     *  "#  media#  Ek/n  R(g/cm2) R(cm)  gbeta  Ek  P/n  Pt  Et"//
     *  " name "
         if(code > 2 ) then
            Ek =  (sqrt( (aPtcl%mass*0.007)**2 + aPtcl%mass**2)
     *              - aPtcl%mass )
               ! upto  gbeta = 2.2
            Ekmax =(sqrt( (aPtcl%mass*2.2)**2 + aPtcl%mass**2)
     *              - aPtcl%mass )
            if( code == 9) then
               Ek = Ek/aPtcl%subcode
               Ekmax = Ekmax/aPtcl%subcode
            endif
         else
            Ek = 1.0d-6         ! 1 keV /n  Ek is k.e/n
            Ekmax = 20.d-3     ! 20 MeV
         endif
         do while (Ek < Ekmax)
           if( aPtcl%code == kgnuc )then
              E =  aPtcl%mass + Ek*aPtcl%subcode
           else
              E =  aPtcl%mass + Ek
           endif
        
           aPtcl%fm%p(4) = E
           call epGetRange(i, Media(i), aPtcl, glen, Rcm)
           if( aPtcl%mass == 0.) then
              betag = 0.
           else
              betag = sqrt(( E/aPtcl%mass)**2 -1)
           endif

           PperN = sqrt(E**2-aPtcl%mass**2)
           if(aPtcl%code == kgnuc ) then
              PperN = PperN/aPtcl%subcode 
           endif
           write(*,'(i3,1p,8g12.3,2x,a)') i,
     *       Ek,  glen, Rcm, betag, E-aPtcl%mass, PperN,
     *       sqrt(E**2-aPtcl%mass**2), E,  Media(i)%name
     *       
           Ek = Ek*10.**0.03
        enddo
      enddo
        write(*,'(a)') 
     *  "#  media#  Ek/n  R(g/cm2) R(cm)  gbeta  Ek  P/n  Pt  Et"//
     *  " name "
      end

