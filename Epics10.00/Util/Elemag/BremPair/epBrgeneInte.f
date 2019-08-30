!           generic brems function to be used when integrating
!     it 

      real*8 function epBrgenex(x)
      implicit none
!       generic brems function ds/dx in mb where 
!       x = Eg/Ee.   
!
!  This must be called after epBPGeneini has been called
! and then  epBrgenePreInte  is called.
!
!       Three functions 
!              epBrSfs  in the low energy region 
!              epBremS  in the partial screening region  
!                          There are three same name  functions in epBPfuncX.f
!                          (X=1,2,3).  One of them must be selected at the
!                          compilation time.
!              epCompScrBr  in the 
!                          copmlete screeing region
!                   are used.  In each region, if the LPM works
!              epBremSH is employed
!       This is the interface routine to call epBregne and can be
!       used when integrating the brems function.
!       The flow is as follows. 
!           1)  select one of the 3 functions depending on Ee
!           2)  if LPMeffect=t & Ee >  media.cnst.BremEemaxL 
!                    see if x is LPM region. (x<Xc)
!           3)  if 2) is yes
!                  compute LPM func at Xc and get normalization
!                  const to the selected function*normfactor
!                  (normfactor is normalization value of each
!                   selected func.)
!                  and compute LPMfunc(x)*normalization const
! 
!            
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"
!

      real*8 x  !  input Eg/Ee. 

      real*8 epBrgene
      real*8 Ee

      Ee= Eeme * masele
      epBrgenex = epBrgene(media, force, Ee, x)
      end function epBrgenex

      subroutine epBrgenePreInte(mediain, forcein, Eemein)
!    this must be called before  epBrgenex
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"
       type(epmedia)::  mediain  ! input media
      real(8),intent(in):: Eemein  ! Ee/me
      character(8),intent(in):: forcein
      media = mediain
      Eeme = Eemein
      force = forcein
      end
