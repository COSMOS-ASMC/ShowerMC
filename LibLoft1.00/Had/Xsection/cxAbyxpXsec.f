!     ratio of SxA/Sxp as a funciton of Sxp and A
!     (inelastic cross section) 
!     Need table in Cosmos/Data/DPM/sigxAbysigxp
!                       ...   /QGS/ ..
!                       ...   /EPOS/ ..
!      if xp < 15, use cxAbyxpXsecOld else use table
!
      function cxAbyxpXsec(xpin, Ain) result(ratio)
      implicit none
#include "Zevhnp.h"
#include "Zmanagerp.h"
      real(8),intent(in):: xpin   ! xp inelastic cross section. (mb)
                        ! x may be any. (pi,K,p, gamma, nu )
      real(8),intent(in):: Ain    ! target mass #. May be non integer.
                            ! A = 2 ~ 210
      real(8):: ratio           !   SigxA/sigxp

      logical,save:: first=.true.
      character(50)::fileDPM="$LIBLOFT/Data/DPM/sigxAbysigxp"
      character(50)::fileQGS="$LIBLOFT/Data/QGS/sigxAbysigxp"
      character(50)::fileEPOS="$LIBLOFT/Data/EPOS/sigxAbysigxp"
      character(50)::file
      integer::icon

      real(4),allocatable,save::XSratio(:,:)
      real(4),allocatable,save::XSa(:)
      real(4),allocatable,save::Aa(:)
      real(4)::xp, A, ratioS,  error
      real(4):: XSv, Av, XSr, Ar

      integer,save::iXS, iA      
      integer::iXSc, iAc
      real(8):: cxAbyxpXsecOld


      real(8):: xsxa
      if( first ) then
       !  read table sigxAbysigxp 
         if( SxAbySxpOpt == 1 ) then
            file =fileQGS
         elseif(SxAbySxpOpt == 2) then
            file =fileDPM
         elseif(SxAbySxpOpt == 3) then
            file = fileEPOS
         elseif(SxAbySxpOpt == 4) then
           !  use old one
            goto 200
         elseif( SxAbySxpOpt == 0) then
            ! use default
            file =fileQGS
         else
            write(0,*) ' SxAbySxpOpt =',SxAbySxpOpt,' invalid'
            write(0,*) ' it must be 1 ~4'
            stop
         endif
         call copenf(TempDev, file, icon)
         if(icon /= 0) then
            write(0,*) file, " could not be opened "
            stop
         endif
         XSv =0.
         Av = 0.
         iXS = 0
         iA = 0
         do while(.true.)
            read(TempDev,*, end=100) XSr, Ar 
            if( XSr > XSv ) then
               iXS = iXS + 1
               XSv = XSr
            endif
            if( Ar > Av ) then
               iA = iA + 1 
               Av = Ar
            endif
         enddo
 100     continue
         rewind TempDev
         allocate( XSratio(iXS, iA) )
         allocate( XSa(iXS) )
         allocate( Aa(iA) )
         do iXSc = 1, iXS
            do iAc = 1, iA
               read(TempDev,*)
     *              XSa(iXSc), Aa(iAc), XSratio(iXSc, iAc)
            enddo
         enddo
         close(TempDev)
 200     continue
         first = .false.
      endif
      if( Ain == 1.0 ) then
         ratio = 1.
         return  !******************
      endif

      if( xpin < 15.d0 .or. SxAbySxpOpt == 4 ) then
         ratio = cxAbyxpXsecOld(xpin, Ain)
      else
         xp = xpin
         A = Ain
         call  kpolintp2S(XSa, 1, 0., Aa, 1, 0.,
     *   XSratio, iXS, iXS, iA,  2, 2, xp, A, ratioS, error)
         ratio = ratioS
         if( ratio < 1. ) then
            ratio = cxAbyxpXsecOld(xpin, Ain)
         endif
      endif
      end   function cxAbyxpXsec

      function cxAbyxpXsecOld(xpin, a) result(ratio)
!       This give old value; for large xpin and small Ain
!       ration seems to be little bit smaller than other
!       code. 
      implicit none
      real(8),intent(in):: xpin ! see cxAbyxpXsec
      real(8),intent(in):: a ! see cxAbyxpXsec
      real(8):: ratio   !  xA/xp cross-section ratio
      real(8):: cinelCosByPdg

      real ca, cb, cc
      if(xpin .lt. 15.) then
         ca = .93812*a**1.0215
         cb = .65000E-02* a**0.5 + (((-.37436E-10* a +.18937E-07)*a 
     *      -.33424E-05)* a + .21646E-03)*a -.31385E-02
         ratio = ca/(1.0 + cb* xpin)
      elseif(xpin .lt. 35.) then
         cb = (0.17859E-05*a-.48592)*a +  .66909E-02  +
     *        0.48968 * a**0.99875       
         ca = (0.67395E-03*a +   .93859 )*a + .26890  +
     *        .34987E-01 * a**1.4525   
         cc = 0.925
         ratio =  ca/(1.0 + cb* xpin**cc)
      else
         ca = 80.000 * a**0.17829  + (((-.52469E-07*a +.30395E-04)*a
     *        -.69823E-02)*a +1.0517)*a +16.427
         if(a .lt. 20.) then
            cb  = 90.029* a**(-1.0219) + 
     *         (.94351E-02*a -.45131)*a+6.4772
         else
            cb = 130.79*a**(-1.0) + 
     *        (-.89803E-05*a +.52418E-02)*a -1.2073
         endif
         ratio = ca* log(1.0 + xpin/cb)/xpin
        endif
        ratio = ratio/cinelCosByPdg(a)
        ratio =  max(1.d0, ratio)
        end
