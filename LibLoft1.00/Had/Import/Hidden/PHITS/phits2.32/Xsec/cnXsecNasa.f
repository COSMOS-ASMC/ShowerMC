      subroutine cnXsecNasa(ap0,zp0,ep0,at0,zt0,signe,sigel,bmax)
c This is to follow the  nasa paper exactly for n projectile case
c  but the result is diff. from the paper
*        calculates reaction cross-sections of nucleus-nucleus         *
*                                                                      *
*        this parametrization is taken from                            *
*        NIM B 129 (1997)  11 by R.K. Tripathi, et.al.,,               *
*                                                                      *
*     input:                                                           *
*        ap     : mass number of projectile                            *
*        zp     : charge number of projectile                          *
*        ep    : incident nucleus total energy (MeV)                  *
*        at     : mass number of target                                *
*        zt     : charge number of target                              *
*                                                                      *
*     output:                                                          *
*       signe   : reaction cross-section (b)                           *
*       sigel   : elastic cross-section (b)                            *
*       bmax    : correspond impact parameter (fm)                     *
*                                                                      *
************************************************************************

      implicit none

      real*8 signe, pi, r0, deltaE, B, Ecm, R, rp, rt, rmsp, rmst,
     $     radius, S, CE, rmsc,  Ek, rap, rat,  D, Rc, fact,
     $     G, T1, SL, X1, Xm, sigel, bmax, Ep, Dg, rca, rac
      real(8):: ap13, at13
      real*8 Ap, At, Zp, Zt
      real*8 Ap0, At0, Zp0, Zt0, ep0
c///////////
c      real(8):: mn
c      parameter(mn = 939)
c/////////
      integer i

      parameter (pi    = 3.14159265358979d0)
      parameter (r0    = 1.1)   ! (fm)

      Ep = ep0
      ap = ap0
      zp = zp0
      at = at0
      zt = zt0
      if(ap  /= 1.0  .or.  zp  /= 0.) then
         write(0,*) ' cnXsecNasa: this is for neutron '
         write(0,*) ' but ap=',ap, ' zp=',zp
         stop
      endif
c----------------------------------------------------------------------

      signe = 0.0
      sigel = 0.0
      bmax  = 0.0

      Ek   = Ep / ap
      Ecm  = Ep * at / ( ap + at )   ! non rela; up to Ecm=2GeV
c///////
c      Ecm = sqrt( ap**2 + at**2 + 2*at*(Ek+ap*mn)/mn)*mn -
c     *   mn*(ap+at)
c///////////
      rmsp  = radius(ap,zp)
      rmst  = radius(at,zt)

      fact = sqrt(5./3.) ! = 1.29

      rp  = fact * rmsp
      rt  = fact * rmst
      rca = fact * 2.471

      rap = 3. * ap  / ( 4. * pi * rp**3. )
      rat = 3. * at  / ( 4. * pi * rt**3. )
      rac = 3. * 12. / ( 4. * pi * rca**3. )


c---------------------------------------------------- D selection start
      

c     neutron + X systems
      ap13 = 1.
      at13 = at**(1./3.)
      S=at13/(1.+at13/ap13)
      if(at >= 11. .and. at <= 40.) then
         T1 =30.
      else
         T1 = 40.
      endif

      if( at < 200.) then
         X1 = 2.83 - 3.1e-2 * at + 1.7e-4 * at * at
         if(at < 12.0) then
            SL = 0.6
         elseif( at == 12.0) then
            SL = 1.6
         else
            SL = 1.0
         endif
         Xm =1. - X1* exp(-Ek/X1/SL)
      else
         Xm = (1.- 0.3*exp( -(Ek-1.0)/15.0)) *
     *        (1. -exp(0.9- Ek))
         SL = -1.   ! not used
         X1 = -1.    ! //
      endif

      Dg = 0.538/( rap + rat)
      if( at <= 40. ) then
         D = Dg- 1.5*(at-2*zt)/at + 
     *        0.25/(1 + exp( (Ek-170.)/100. ))
      elseif( at < 60. ) then
         D = Dg-1.5*(at-2*zt)/at
      elseif( zt > 82.) then
         D = Dg - zt/(at-zt)
      else
c            all others
         D = Dg
      endif
      CE = D* (1.0 - exp( -Ek/T1)) 
     *        - 0.292*exp( -Ek/792.0)*
     *          cos(0.229*Ek**0.453)
      deltaE = 1.85 * S + ( 0.16 * S / Ecm **(1./3.) ) - CE 
c     * +  0.91*(at - 2 * zt ) * zp / ( at * ap )    ! since zp =0


      signe = pi * r0 * r0 * ( ap13 + at13 + deltaE )**2
     *  * Xm
      if( signe .le. 0.) then
         signe=0
      endif
c////////////
c      write(0,'(1p,9g12.4)' )
c     * Ek, T1, D, CE, deltaE, S, SL, Xm, X1
c////////////////

      signe = signe * 10        ! convert the unit from [fm**2] to [mb]
      signe = signe / 1000      ! convert the unit from [mb] to [b]

      bmax =  r0 * ( ap13 + at13 + deltaE )
      end subroutine
