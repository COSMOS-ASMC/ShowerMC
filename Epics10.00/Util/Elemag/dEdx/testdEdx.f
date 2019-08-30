#include "ZcosmosBD.h"
#include "ZepicsBD.h"
      implicit none
!           next v9.162 needs next for Media
!  #include "ZepTrackp.h"   
#include "ZepTrackv.h"
#include "Zmass.h"
#include "Zcode.h"
!
!          test  dE/dx formula
! 
       type(epmedia):: mediax
       type(ptcl)::  aPtcl
      
      integer  subcode, charge, code, nr
      integer,save:: cmpn=1
      real*8 E,  dedx, w, betag, Ek, PperN, dedxfull
      character*20 epifile, config
      character(256):: input 
      character(20):: b(5) ! 6-->5 v9.162
      real(8)::minErg(5)  ! v9.162


      read(*,'(a)') input
      write(0,*) ' input =', trim(input)
      call ksplit(input, 20, 5, b,  nr) !  6 --> 5 v9.162
      if( nr /= 5 ) then     ! 6-->5 v9.162
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

      Cn = cmpn ! only one comp is assumed. This common variable
              ! must be set  v9.162
      call epSetEmin(cmpn)   ! component dependent Emin,  
      call epqCn2Media(cmpn, mediax)  
      write(0,*) ' media name is ' , mediax%name
      write(0,*) ' mediax%sh%w0=', mediax%sh%w0, 
     *  ' MeV. should be the same as the next one'
      call epqEmin(Cn, minErg)   !  inquire RecoilKEmin
      write(0,*) ' Knockon min K.E =', minErg(3)*1.e6, 'keV'
!           make incident
      call cmkptc(code, subcode, charge, aPtcl)

      Ek = 300d-6 ! 100 keV /n  Ek is k.e/n
      write(0,'(a,a)')
     *  '# beta*g  dE/dx (Rest Full,MeV/(g/cm2))   p(GeV) ',
     *  '   Et(GeV)   K.Et(GeV)   K.E/n(GeV)   p/n(GeV/c)'
      write(0,'(a,a)') 
     *  '#   1             2      3                  4    ',
     *  '     5          6          7            8 '   
      do while ( Ek .lt. 10000.)
         if(aPtcl%code == kgnuc) then
            E= Ek*aPtcl%subcode + aPtcl%mass
         else
            E = Ek + aPtcl%mass
         endif
         aPtcl%fm%p(4) = E
         if(aPtcl%code .eq. kelec) then
            call epdedxe(mediax, aPtcl, dedx, dedxfull)
         elseif(aPtcl%code == kgnuc) then
            call epdedxhvy(mediax, aPtcl, dedx, dedxfull )
         else
            call epdedxNone(mediax, aPtcl, dedx, dedxfull)
         endif
         betag = sqrt(( E/aPtcl%mass)**2 -1)
!          betag = p/m =sqrt( (Ek+m)**2-m**2)/m
!                      =sqrt(  Ek*(Ek+2m) )/m
         PperN = sqrt(E**2-aPtcl%mass**2)
         if(aPtcl%code == kgnuc ) then
            PperN = PperN/aPtcl%subcode 
         endif
         write(*,'(1p,8g13.4)')
     *       betag, dedx*1000., dedxfull*1000.,
     *       sqrt(E**2-aPtcl%mass**2), E,
     *       E-aPtcl%mass, Ek,  PperN
         Ek = Ek*10.**0.01 
      enddo
      end
