!  test program of cxAbyxpXsec.  
#include "ZcosmosBD.h"
      program main
      implicit none
 !      next is for enabling change of the default (=1) of
 !         SxAbySxpOpt 
#include "Zevhnp.h"  
#include "Zmass.h"

      real(8):: xp, A, ratio, ratio200
      real(8):: cxAbyxpXsec
      real(8)::  cPDGsigmaInepA
      real(8)::  Ek, p, xs
      real(8):: AA(30)
      integer:: i
      write(0,*) ' Enter SxAbySxpOpt (1~4) '
      write(0,*) ' 1-->QGS like xsec(A)/xsec(p)'
      write(0,*) ' 2-->DPM like xsec(A)/xsec(p)'
      write(0,*) ' 3-->EPOS like xsec(A)/xsec(p) (A<=56)'
      write(0,*) ' 4-->old formula'
      read(*,*) SxAbySxpOpt
      write(0,*) " enter A's (max of 30) with / at last"
      AA(:) = 0.
      read(*,*) AA(:)

      Ek =200.
      p =sqrt( (Ek+masp)**2 - masp**2 )
      call cppInelaXs(p, xs)

      
      do xp= 15., 110, 3.
         do i= 1, 30
            A = AA(i)
            if(A <=1) exit
            ratio200 =cPDGsigmaInepA(A)/xs
            ratio = cxAbyxpXsec(xp, A)
            write(*,*) "E ", xp, A, ratio
            write(*,*) "200 ",  xs, A, ratio200
         enddo
      enddo
      end
