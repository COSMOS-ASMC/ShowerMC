!          compute  K-N formula
!         ds/dcos or ds/dEg' (mb) or other normalization
      implicit none
#include "Zmedia.h"
#include "Zmass.h"

       type(epmedia):: media
      integer  io, icon
      real*8 Ee,  Eg, cosg,  Egp,  v, dsdv, dsdcos
      real*8 Egpmin, tprob, path, Nc
      character*120 file
      character(len=8):: name
      integer norm

      io = 10
      write(0,*)
     *  'Enter Eg(GeV), normalization, media name'
      write(0,*) 
     *"   where norm: 1->/r.l; 2->No norm; 3->/(g/cm2);"//
     * "  4->/cm; 5->area"
      read(*, *) Eg, norm, name
      file ="$EPICSTOP/Data/Media/"//trim(name)

      call copenf(io, file, icon)
      call epReadMTbl(io, media)
      close(io)
      call epGetEffZA(media)
      Egpmin = Eg/(1.+2.0*Eg/masele)
      Egp = Egpmin
      call epcompp(media, Eg, tprob, path)
      tprob = tprob/media%mbtoPX0   ! in mb. to be consistent with other
                               ! treatment in BremPair
      if(norm .eq. 5) then
         Nc = 1./tprob
      elseif(norm .eq. 1) then
         Nc=media%mbtoPX0  ! prob/X0                                            
      elseif( norm .eq.  2) then
         Nc = 1
      elseif(norm .eq.  3) then
         Nc = media%mbtoPgrm
      elseif(norm .eq. 4) then
         Nc = media%mbtoPcm
      else
         call cerrorMsg('input error for norm',0)
      endif

      do while (.true.)
         call epcomptonFunc(media, 1, Eg, Egp, cosg, dsdv)
         dsdv = dsdv*Nc
         dsdcos = dsdv *Eg/masele /( 1.+Eg/masele*(1.-cosg))**2
         v = Egp/Eg
         write(*,'(1p,8g12.4)') 
     *    v, dsdv, cosg, dsdcos, Eg, Egp, Nc, tprob
         Egp = Egp*10.**0.02
         if(Egp .gt. Eg) exit
      enddo
      end


      subroutine epcomptonFunc(media, j, Eg, Egp, cosg, f)
      implicit none
!          ds/dcos or ds/dEg' (mb)
#include "Zglobalc.h"
#include "Zmedia.h"
#include "ZbasicCnst.h"
#include "Zmass.h"

       type(epmedia)::  media ! input
      integer j ! input. 1: ==>  f = ds/dv  v=Egp/Eg  (mb)
                !        2: ==>  f = ds/dcos  (mb)    

      real*8  Eg ! input. incident photon energy (GeV)
      real*8  Egp ! input. if j=1.  scattered photon energy (GeV)
                  ! output if j=2.   //
      real*8  cosg ! input if j=2.  //      angle
                   ! output if j=1 ///

      real*8  f   ! output.
      real*8 dsdomega, norm

      if(j .eq. 1 ) then
!      Egp is given
         cosg = 1. - (Eg/Egp-1.)*masele/Eg
         norm =Eg/masele/(1.+Eg/masele*(1.-cosg))**2
      else if( j .eq. 2) then
         Egp = Eg/(1. +Eg/masele*(1.-cosg) )
         norm = 1.
      else
         write(0,*) ' error j=',j, 'in drawCompFunc.f'
      endif
      dsdomega = 
     * pir02 * (Egp/Eg)**2 *
     *( Eg/Egp + Egp/Eg -(1.-cosg**2))*media%Z 
      
      f = dsdomega/norm
      end
