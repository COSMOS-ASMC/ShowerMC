      subroutine epGetEffZA(media)
      implicit none
#include  "Zglobalc.h"
#include  "ZbasicCnst.h"
#include  "Zmedia.h"

       type(epmedia):: media  ! input.  should contain basic media info. 
                           ! output. other info is computed here.

      integer i, nc, j
      real*8 sum, Lrad, Lradp, f, cs
!
!     &&&&&&&&&&&&&&&&&&&&&&&&&&         
!       some trial to get Zeff so that it (as a sigle virtual atom)
!       gives the same r.l obtained for complex material.
!       This gives almost the same Zeff as obtained in this prog.
!       (normally within 1 %).
!&      record /media/mediax
!&      common /Zmedia/ mediax
!&      real*8 newz
!&      integer icon
!&      real*8 epsolverl
!&      external epsolverl
!     &&&&&&&&&&&&&&&&&&&&&&&&&&&

      nc = media%noOfElem
      sum = 0.
      do i = 1, nc
         sum = sum + media%No(i) 
      enddo
      media%sumNo = sum
! &&&&&&&&&& from  v7.0 we use normalized No(i).  and this leads to  the 
!      change of the  cross-section (in mb) (new one  is  old one / sum(old No)).
!      and also change of mbto* constant.  However, there is no change
!      in the sampling table and actual simulations.
!
      do i = 1, nc
         media%OrigNo(i) = media%No(i)    ! 
         media%No(i) = media%No(i)/sum
      enddo
      sum = 0.
      media%A = 0.
      media%Z = 0.
      media%Z2 = 0.
      media%ZZ1 = 0.
      media%MoliereExpb = 0.
      media%Z1_3rd =0.
      media%Z2_3rd =0.

      do i = 1, nc
         media%A = media%A + media%No(i)* media%elem(i)%A
         media%Z = media%Z + media%No(i)* media%elem(i)%Z
         media%Z2 = media%Z2 + media%No(i)* media%elem(i)%Z**2
         media%ZZ1 = media%ZZ1 + 
     *     media%No(i)* media%elem(i)%Z*(media%elem(i)%Z + 1.)
!            for beta=1 and z=1; Moliere's  exp(b) = MoliereExb*t
         media%MoliereExpb = media%MoliereExpb +
     *       media%No(i)*media%elem(i)%Z**0.333333 *
     *        (media%elem(i)%Z+1.d0)/
     *              (1.0d0 + 3.327d0*(media%elem(i)%Z/137.d0)**2)
         media%Z1_3rd =  media%Z1_3rd 
     *          +  media%No(i)* media%elem(i)%Z**(1./3.)
         media%Z2_3rd =  media%Z2_3rd 
     *          +  media%No(i)* media%elem(i)%Z**(2./3.)
      enddo
!                only for e+/e- 
      media%MoliereForXc2 =0.6011 * media%ZZ1/media%A

      media%MoliereExpb = 6702. * media%MoliereExpb/media%A
!        we simply deinfe  Z^(2/3) by <Z>^(2/3) (not <Z^2/3>
!        also  Z^(1/3)
!      media.Z1_3rd = media.Z**(1./3.)
!      media.Z2_3rd = media.Z**(2./3.)

      do i = 1, nc
         media%w(i) = media%No(i)* media%elem(i)%A/media%A
         media%npercm3(i) =N0/media%elem(i)%A *media%w(i)*media%rho
      enddo

      sum = 0.
!
!         Since many of cross-sections depend primarily on
!         Z**2/A, we take effective Z**2/A first and  
!         multiply it by Aeff to get effective Z**2 and
!         further get effective Z by sqrt(Zeff**2)
!         For the cross-sections which depend on Z/A,
!         we use effective Z/A directly.
      sum = 0.
      do i =1, nc
         sum = sum + media%w(i)*media%elem(i)%A
      enddo
      media%Aeff = sum    ! not to be used alone

!         get effective Z**2/A
      sum = 0
      do i = 1, nc
         sum = sum + media%w(i)*media%elem(i)%Z**2/
     *         media%elem(i)%A
      enddo
      media%Z2byAeff = sum   ! OK  = sum(niZi^2)/sum(niAi)

      media%Z2eff = media%Z2byAeff * media%Aeff  ! not to be used alone
      media%Zeff = sqrt(media%Z2eff)         ! not to be used  

!  &&&&&&&&&&&&&&&
!      write(*,*) ' meida.No', ( media.No(i), i=1,nc )
!      write(*,*) ' meida.w ', ( media.w(i), i=1,nc )
!      write(*,*) ' A,Z, Aeff, Zeff=', media.A, media.Z, media.Aeff,
!     *   media.Zeff
! &&&&&&&&&&&&&&&&
      media%LogZ = log(media%Zeff)  ! not used anywhere
      media%Zeff3 = media%Zeff**(1.d0/3.d0)   ! not to be used

      sum = 0.
      do i = 1,nc
         sum  =  sum +
     *             media%w(i)*media%elem(i)%Z/media%elem(i)%A
      enddo
!         =  sum(niZi)/sum(niAi)  OK
      if( media%format .eq. 2 .and. media%ZbyAeff .ne. -100.) then
         if(abs(sum-media%ZbyAeff)/media%ZbyAeff .gt. 0.01 ) then
            write(0,*) "weighted <Z/A> inp=",media%ZbyAeff,
     *      " computed inside =", sum
            write(0,*) "former will be used"
         endif
      else
         media%ZbyAeff = sum
      endif
!            for photo electric effect
      sum = 0.
      do i = 1, nc
         sum = sum + 
     *        media%w(i)*media%elem(i)%Z**5/media%elem(i)%A
      enddo
      media%Z5byAeff = sum   ! ok


!         compute radiation length X0g(g/cm2) and  X0 (cm)
      call epX0(media)
      media%X0m = media%X0/100.d0  ! 1 r.l in m
!         g/cm^2 to cm  for const den.
      media%gtocm = media%X0/media%X0g
      media%X0kg = media%X0g *10.d0 !  X0 in kg/m2
!        kg/m2 to m    for const den.    X0 1e-2/ (X0g10)
      media%kgtom = media%gtocm*1.0d-3 !
!     Air:    X0-36.5g/cm2 ~ 30000 cm =300 m  gtocm~ 30000/36.5
!                    kgtom ~ 3e4/36.5/1e3~ 30/36.5 
!             X0kg =  365 kg/m2 -300 m    ; 1kg/m2  ~ 300/365 m  ok

!           this may be used when summing each
!           element contribution
      media%mbtoPgrm = 1.d-27 *N0/media%A

!     mfp = 1.0d0/( media%mbtoPkgrm * media%xs)   in kg/m2
!          if mfp is 100 g/cm2 it is 1000 kg/m2 so next coef is 0.1
      media%mbtoPkgrm = media%mbtoPgrm * 0.1d0

      media%mbtoPcm =  media%mbtoPgrm * media%rho
      media%mbtoPX0 =  media%mbtoPgrm *media%X0g

!         <Z/A> must be used explicitly for crosssections with Z
      media%basearea = pir02 * media%ZbyAeff * N0*1.d-27 * media%X0g
!
!           this may be used when we approximate
!      a compound or molecule as one element
!        wrong: ! not used anywhere
!
      media%mbtoPgrm2 = 1.d-27*N0/media%Aeff   
      media%mbtoPcm2 =  media%mbtoPgrm2 * media%rho
      media%mbtoPX02 =  media%mbtoPgrm2 * media%X0g

! set approx Ecrit;  from PDG.  For Air,
!     the Ecrit is used only for rough calculation of Moliere
!     unit and hybrid A.S size calculation.
!     For the rough estimation of casade spread,
!     the scattering constant, Es, is also used .
!     Es=19.3 MeV is said to be better than tradistional Es=21 MeV
!     but could lead to bit too small spread if we use
!     the value below.  (Air Ec=81 MeV in older versions. and Es=19.3
!     MeV. ).  The pdg value for Air (gas), ~86.6 MeV will be obtained.
!     In old days, Ec=76, or 84.2 MeV was standard in A.S community. But
!     81 is recommended by Linsley (some communication with him)
!     Probably 81 may be better. To use it the use must
!     change Media(MediaNo)%Ecrit 
!      
!      For  muons, they are not yet used.       
      if( media%gasF == 0 ) then
!     solid
         media%Ecrit = 0.61/(media%Z + 1.24)
         media%Ecritmu = 5700.0/(media%Z+1.47)**0.838
      else
!         gas
         media%Ecrit = 0.71/(media%Z + 0.92)
         media%Ecritmu = 7980.0/(media%Z+2.03)**0.879         
      endif


!  &&&&&&&&&&&&&&&&&&&&&&&&
!        get Zeff for compound/molecule so that Zeff and Aeff as a
!        single atom gives  the same X0 as above.
!&      mediax = media
!&      write(*,*) ' X0=', media.X0, ' X0g=',media.X0g
!&      write(*,*) ' Zeff=',media.Zeff, ' Aeff=',media.Aeff
!&      call kbchop(epsolverl, 1.d0,  140.d0, 0.001d0, newz, icon)
!&       write(*,*) ' newz=',newz
!  &&&&&&&&&&&&&&&&&&&&&&&
!       ------------LPM -----------;
!        these are not atom base. but media base.
!      media.s1 = (media.Zeff**(1.d0/3.d0)/183.0d0)**2
      media%s1 = (media%Z1_3rd/183.0d0)**2
      media%logs1 = log(media%s1)
!       -----------complete screening cross-sec. coeff.

!        const used in the complete screening cross-sec.
!      f(y)dy ~dy (4/3(1-y)/y + y)C1 + (1-y)/y C2
!
      media%cScrC1 = 0.
      media%cScrC2 = 0.
!

      do i = 1, media%noOfElem
         call epGetLrad(media%elem(i)%Z, Lrad, Lradp, f)
         media%cScrC1 =  media%cScrC1  + 
     *     ( media%elem(i)%Z**2 *(Lrad-f) +
     *            media%elem(i)%Z*Lradp) * media%No(i)
         media%cScrC2 = media%cScrC2  + 
     *      (media%elem(i)%Z**2  + media%elem(i)%Z)/9.d0
     *       * media%No(i)
      enddo
!       f(y) ~ (4/3C1+ C2)(1-y)/y + yC1; main term
      media%cScrMain = 4.*media%cScrC1/3.d0 + media%cScrC2
!          plasma energy in GeV
      media%wp = 28.816d-9 *sqrt(media%rho*media%Z/media%A)
!        no. of ingredient /cm^3
      media%nd = media%rho*N0/media%A
!
!           make cross-section table: for 10 mb to 100 mb
!
      media%ndensity = 1.d-27*media%rho*N0/media%A
      end

      subroutine epResetEcrit(io, name, newV, oldV, icon)
!       to reset Crittical energy of a given media with "name"
      implicit none
#include "ZmediaLoft.h"
      integer,intent(in):: io   ! output message device #
           !  0-->  some message is put as standard Fortran error message  
           !  6-->  some message is put as sysout.
           ! >0 --> assume logical device is open with that number
           !  <0 --> no message is put,  but see next
      character(*),intent(in):: name ! media name such as "Air"
          ! if media with "name" is not found eroor message is
          ! 
      real(8),intent(in):: newV ! new critical energy (GeV)
       ! new value is set to media%Ecrit
      real(8),intent(out):: oldV ! E crit  so far defined. (GeV)
      integer,intent(out):: icon ! 0 if ok.  -1 if some error
      integer:: i

      icon = -1
      do i = 1, NoOfMedia
         if( trim(Media(i)%name) == trim(name) ) then
            oldV = Media(i)%Ecrit
            Media(i)%Ecrit = newV
            icon = 0
            exit
         endif
      enddo
      if(icon == -1 ) then
         write(0,*)
     *     'epResetEcrit: input media name ', trim(name),
     *     ' not found in Meida: availabe median names are:'
         do i = 1, NoOfMedia
            write(0,*) Media(i)%name
         enddo
         call cerrorMsg('none of them matches',0)            
         stop
      elseif( io >= 0 ) then
         write(io,
     *        '("Critical Energy for media ",a," was reset to")')
     *        trim(name)
         write(io, '(f6.1, "MeV from old value ", f6.1," MeV")')
     *        newV*1.0d3, oldV*1.0d3
      endif
      end
