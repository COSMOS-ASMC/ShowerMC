!     ******************************************************************
!     *                                                                *
!     * chookHybAS. There are two routines here
!     *     chookHybAS: is the interface when a component A.S has been
!     *                made from an electron.  
!     *     chookHybAS2: is the interface which may be called at the 
!     *                 end of one event generation for air shower business.
!     *   Note: The former is called from the system.
!     *         The latter must be called in your chookEnEvent, if
!     *         necessary.  These are template rouitnes.  You must/may
!     *         modify them.  The latter may be included in your
!     *         chookEnEvent routine directly. The latter name 
!     *         can be another one.
!     *                                                        
!     ******************************************************************
!
!
!
      subroutine chookHybAS(el, never)
      implicit none

!      This routine is called at the end of the cobAS.f in Tracking/AS/,
!      that is, when a component A.S is made from an electron.  If the
!      user has nothing to do, give never=1,  then the
!      routine will not be called again.
!      If you have something to do (say, business for 
!      generating air fluorescence light), do it here.
!      Give never=0 in such a case.

#include "Ztrack.h"
! #include "Zmagfield.h"      
#include "Zobs.h"
#include "Zobsv.h"
#include "Zelemagp.h"
#include "ZmediaLoft.h"
      type(track):: el    ! input. an electron produced  component A%S
      integer never        ! input /output.  give 1 if you don't need
                           ! the routine else give 0.   

!     The following will be the typical stuff you may want to use in
!     this routine. (i=1, NoOfASSites; index for observation depths)
!
!   CompASNe(i):  component A.S size produced by the input electron.
!                 For depths where this value is 0, 
!                 avoid doing something here.  
!   CompASAge(i): age of component A.S produced by the input electron.
!                 If this value is 2.0, the A.S is assumed to be very
!                 old and the CompASNe(i) is 0.  You should skip
!                 treating deeper depths.
!

      real*8  zobas, zp
      real*8  elog, eno, age
      real*8  cvh2temp, tk  ! temperature in Kelvin
      real*8  dedx  ! to store <dE/dx> 
      real*8  cvh2den, rho  ! density of air in kg/m^3

      integer xsite
!
!     **********
      never = 1        ! change this to 0 if you need this routine
!     **********

      zp = el%pos%depth      ! starting vertical depth of
                             ! the component electron (kg/m^3)
!
!        get average dE/dx for every depth.
!


      do   xsite = 1, NoOfASSites
         age = CompASAge(xsite) 
         eno = CompASNe(xsite)


         if(age .eq. 2.) then
!            store 0 or ... in your own array 

            goto 100
         endif
         if(eno .gt. 0.) then
!      
!            you may do some business. Say  generate fluorescence light.
!            you may need some array to store the quantities
!            you compute here. (presumably in your own common block).
!        Following is typical quantities you may need for such a
!        computation
!             temerature in Kelvin of 'site'
            tk = cvh2temp(ASObsSites(xsite)%pos%height)
!             density of air in kg/m^3;  multiply 10^-3 to get it in
!             g/cm^3.
            rho = cvh2den(ASObsSites(xsite)%pos%height)  
            call cavedEdx(CompASNe(xsite), CompASAge(xsite), dedx)
!          ****** dedx >0 and in GeV/(kg/m^2). To convert it to
!          ****** MeV/(g/cm^2).  Multiply 100. 

!             vertical depth of site(  kg/m2)
            zobas=ASObsSites(xsite)%pos%depth
            
!              log10 of elecrton energy  in terms of critical energy
            elog = log10(el%p%fm%p(4)/Media(MediaNo)%Ecrit) 
            
         endif
      enddo
  100 continue
      end
!
      subroutine chookHybAS2
      implicit none
!     
!           You may utilize this routine for computing, say,
!      air fluo. light.  for a given air shower.  This should
!      be called from chookEnEvent routine.
!
!
#include "Ztrack.h"
! #include "Zmagfield.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"



      real*8  eno, age,  zobas, nmu
      real*8  muonno(maxNoOfASSites)
      real*8  cvh2temp, tk  ! temperature in Kelvin
      real*8  dedx  ! to store <dE/dx> 
      real*8  eth
      real*8  cvh2den, rho  ! density in  kg/m^3.
!
      data eth/1.0/   ! dummy Emu threshold

      integer xsite
!
!      ******************** Below:  not usable for a while
!                                 (as of May/10,'97)
!                Nmu (E>Eth)
       call cgetNmu(eth, muonno)
!      ********************

      if(ObserveAS) then
!
         do   xsite = 1, NoOfASSites
            age = ASObsSites(xsite)%esize
            eno =  ASObsSites(xsite)%age
            nmu = muonno(xsite)
!             Ne or Nmu > 0
            if(eno .gt. 0. .or. nmu .gt. 0.)  then
!      
!            you may do some business. Say  generate fluorescence light.
!            you may need some array to store the quantities
!            you compute here. 
!        Following is typical quantities you may need for such a
!        computation
!             temerature in Kelvin of 'site'
               tk = cvh2temp(ASObsSites(xsite)%pos%height)
!                 in kg/m^3; x 10^-3 in g/cm^3
               rho = cvh2den(ASObsSites(xsite)%pos%height)

!                get average <dE/dx>at site            
               call cavedEdx(ASObsSites(xsite)%esize, 
     *          ASObsSites(xsite)%age,  dedx)
!      ******* dedx >0 and in GeV/(kg/m^2). To convert it to
!      ******* MeV/(g/cm^2).  Multiply 100. 

!             vertical depth of site(  kg/m2)
               zobas=ASObsSites(xsite)%pos%depth
            endif
         enddo
      endif
      end









