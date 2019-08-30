!     ---------------------------------------------------------------------
!          cgettauc:  get (mean life time) x (light velocity)  of
!      a given particle 
!
      subroutine cgetctau(proj, ctau)

      implicit none

#include  "Zglobalc.h"
#include  "Zlife.h"
#include  "Zptcl.h"
#include  "Zcode.h"
!     **************************************************
      type(ptcl)::proj
      real*8 ctau
!
!
      real*8 tcpi, tck, tcmu, tcks, tckl, tcpi0, tcd0, tcdc, 
     *       tceta, tcsigma0, tcsigmap, tcsigmam, tcgzai0,
     *       tcgzaim, tcbomega, tclambda, tclambdac

      parameter (tcpi=t0pi*c, tck=t0k*c,
     *        tcmu=t0mu*c, tcks=t0ks*c,
     *        tckl=t0kl*c, tcsigma0=t0sigma0*c,
     *        tcsigmap = t0sigmap*c, tcsigmam=t0sigmam*c,
     *        tcgzai0 = t0gzai0 * c, tcgzaim=t0gzaim*c,
     *        tclambda =t0lambda*c, tclambdac=t0lambdac*c,
     *        tcbomega = t0bomega*c,
     *        tcpi0=t0pi0*c, tcd0=t0d0*c,
     *        tcdc =t0dc*c, tceta= t0eta*c)

      real(8),parameter:: tcds=t0ds*c, tctau=t0tau*c, tcXic=t0Xic*c,
     *     tcomeC0 = t0omeC0*c, tcetap= t0etap*c, tcXic0=t0Xic0*c
!     *     tcrho = t0rho*c,
!     *     tcomega = t0omega*c, tcphi=t0phi*c 
            
      
      if(proj%code .eq. kpion) then
          if(proj%charge .ne. 0) then
             ctau = tcpi
          else
             ctau = tcpi0
          endif
      elseif(proj%code .eq. kkaon) then
         if(proj%charge .ne. 0) then
            ctau = tck
         elseif(abs(proj%subcode) .eq. k0l) then
            ctau = tckl
         elseif(abs(proj%subcode) .eq. k0s) then
            ctau = tcks
         else
            write(0,*) ' K0 subcode=', proj%subcode, ' stragne'
            stop
         endif
      elseif(proj%code .eq. kmuon) then
         ctau = tcmu
      elseif(proj%code == knuc ) then
         ctau = Infty
      elseif(proj%code .eq. knnb) then
         ctau =0.
      elseif(proj%code .eq. kddb) then
         ctau = 0.
      elseif(proj%code .eq. kdmes) then
         if(proj%charge .eq. 0) then
            ctau = tcd0
         else
            ctau = tcdc
         endif
      elseif(proj%code .eq.  ksigma) then
         if(proj%charge .eq. 0) then
            ctau = tcsigma0
         elseif(proj%charge .eq. 1) then
            if(proj%subcode .eq. regptcl) then
               ctau = tcsigmap
            else
               ctau = tcsigmam
            endif
         else
            if(proj%subcode .eq. regptcl) then
               ctau = tcsigmam
            else
               ctau = tcsigmap
            endif
         endif
      elseif(proj%code .eq. kgzai) then
         if(proj%charge .eq. 0) then
            ctau = tcgzai0
         else
            ctau = tcgzaim
         endif
      elseif(proj%code .eq. klambda ) then
         ctau = tclambda
      elseif(proj%code .eq. klambdac) then
         ctau = tclambdac
      elseif(proj%code .eq. kbomega) then
         ctau = tcbomega
      elseif(proj%code .eq. krho) then
         ctau = 0.
      elseif(proj%code .eq. kphi) then
         ctau =0.
      elseif(proj%code .eq. komega) then
         ctau = 0.
      elseif(proj%code .eq. keta ) then
         ctau = tceta
      elseif(proj%code .eq. kds  ) then
         ctau = tcds
      elseif( proj%code .eq. ketap  ) then
         ctau = tcetap
      elseif( proj%code .eq. kDelta  ) then
         ctau = 0
      elseif( proj%code .eq. kXic  ) then
         if( proj%charge == 0) then
            ctau = tcXic0
         else
            ctau = tcXic
         endif
      elseif( proj%code .eq. komeC0  ) then
         ctau = tcomeC0
      elseif( proj%code .eq. ktau  ) then
         ctau = tctau
      else
         ctau = Infty
      endif
      end
