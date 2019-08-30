************************************************************************
*                                                                      *
      subroutine nasa(ap0,zp0,ep0,at0,zt0,signe,sigel,bmax)
*                                                                      *
*        calculates reaction cross-sections of nucleus-nucleus         *
*                                                                      *
*        this parametrization is taken from                            *
*        NIM B 117 (1996) 347 by R.K. Tripathi, et.al.,                *
*        NIM B 129 (1997)  11 by R.K. Tripathi, et.al.,,               *
*        and                                                           *
*        NIM B 155 (1999) 349 by R.K. Tripathi, et.al.                 *
*                                                                      *
*        coded by H. Iwase, on 2004/02/06                              *
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
     $     radius, S, CE, rmsc, rca, E, rap, rat, rac, D, Rc, fact,
     $     G, T1, SL, X1, Xm, sigel, bmax, ep

      real*8 Ap, At, Zp, Zt, Ac
      real*8 Ap0, At0, Zp0, Zt0, ep0
      integer i

      parameter (pi    = 3.14159265358979d0)
      parameter (r0    = 1.1)   ! (fm)
      parameter (Ac    = 12.)

c----------------------------------------------------------------------
c     this routine requires that at > ap
c----------------------------------------------------------------------

      if( ap0 .gt. at0 ) then

         ep = ep0 / ap0 * at0
         ap = at0
         zp = zt0
         at = ap0
         zt = zp0

      else

         ep = ep0
         ap = ap0
         zp = zp0
         at = at0
         zt = zt0

      end if

c----------------------------------------------------------------------

      signe = 0.0
      sigel = 0.0
      bmax  = 0.0

      E    = Ep / Ap
      Ecm  = Ep * At / ( Ap + At )

      rmsp  = radius(Ap,Zp)
      rmst  = radius(At,Zt)

      fact = sqrt(5./3.) ! = 1.29

      rp  = fact * rmsp
      rt  = fact * rmst
      rca = fact * 2.471

      rap = 3. * Ap  / ( 4. * pi * rp**3. )
      rat = 3. * At  / ( 4. * pi * rt**3. )
      rac = 3. * 12. / ( 4. * pi * rca**3. )

c---------------------------------------------------- D selection start

      if(ap.eq.1.and.zp.eq.0)then
c     neutron + X systems
         T1 = 18.
         D  = 1.85 + ( 0.16 ) / ( 1. + exp( (500.-E) / 200. ))

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.4.and.zt.eq.2)then
c     proton + alpha
         T1 = 40.
         D  = 2.05

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.3.and.zt.eq.2)then
c     proton + 3He
         T1 = 58.
         D  = 1.70

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.6.and.zt.eq.3)then
c     proton + 6Li
         T1 = 40.
         D  = 2.05

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.7.and.zt.eq.3)then
c     proton + 7Li
         T1 = 37.
         D  = 2.15

      elseif((ap.eq.1.and.zp.eq.1).and.At.le.7)then
c     proton + light nucleus
         T1 = 23.
         D  = 1.85 + ( 0.16 ) / ( 1. + exp( (500.-E) / 200. ))

      elseif(ap.eq.1.and.zp.eq.1)then
c     proton + others
         T1 = 40
         D  = 2.05

      elseif(ap.eq.2.and.zp.eq.1.and.at.eq.4.and.zt.eq.2)then
c     deuteron + alpha
         T1 = 23
         D  = 1.65 + ( 0.22 ) / ( 1. + exp( (500.-E) / 200.  ))

      elseif(ap.eq.2.and.zp.eq.1)then ! deuteron + X systems
         T1 = 23.
         D  = 1.65 + ( 0.1 ) / ( 1. + exp( (500.-E) / 200.  ))

      elseif(ap.eq.3.and.zp.eq.2 .or. at.eq.3.and.zt.eq.2)then
         T1 = 40.
         D  = 1.55


      elseif((ap.eq.4.and.zp.eq.2).and.(At.eq.4.and.Zt.eq.2))then
c     alhpa + alpha
         T1 = 40.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.4))then
c     alpha + Be
         T1 = 25.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.7))then
c     alpha + N
         T1 = 40.
         G  = 500.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.13))then
c     alpha + Al
         T1 = 25.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.26))then
c     alpha + Fe
         T1 = 40.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif(ap.eq.4.and.zp.eq.2 .or. at.eq.4.and.zt.eq.2)then
c     alpha + other systems
         T1 = 40.
         G  = 75.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif(zp.eq.3 .or. zt.eq.3)then
c     Li + X
         T1 = 40.
         D = 1.75 * ( rap + rat ) / ( rac + rac ) /3

      else
c     others
         T1 = 40.
         D = 1.75 * ( rap + rat ) / ( rac + rac )

      endif

c---------------------------------------------------- D selection end


c----------------------------------------------------------------------

      CE = D * ( 1 - exp( -E / T1) )  - 0.292 * exp( -E / 792 ) *
     $     cos( 0.229 * E**(0.453) )

      S = ( Ap**(1./3.) * At**(1./3.) ) / ( Ap**(1./3.) + At**(1./3.) )

      deltaE = 1.85 * S + ( 0.16 * S / Ecm **(1./3.) ) - CE + 0.91 * (
     $     At - 2 * Zt ) * Zp / ( At * Ap )


      R = rp + rt + 1.2 * ( Ap**(1./3.) + At**(1./3.) ) / ( Ecm**(1./3.)
     $     )


      B = 1.44 * Zp * Zt / R

c----------------------------------------------------------------------


c--------------------------------------------------- Rc selection start

      if((ap.eq.1.and.zp.eq.1).and.(at.eq.2.and.zt.eq.1))then
         Rc = 13.5              ! p + d

      elseif((ap.eq.1.and.zp.eq.1).and.(at.eq.3.and.zt.eq.2))then
         Rc = 21                ! p + 3He

      elseif((ap.eq.1.and.zp.eq.1).and.(at.eq.4.and.zt.eq.2))then
         Rc = 27                ! p + 4He

      elseif((ap.eq.1.and.zp.eq.1).and.(zt.eq.3))then
         Rc = 2.2               ! p + Li

      elseif((ap.eq.1.and.zp.eq.1).and.(zt.eq.6))then
         Rc = 3.5               ! p + C

      elseif((ap.eq.2.and.zp.eq.1).and.(at.eq.2.and.zt.eq.1))then
         Rc = 13.5              ! d + d

      elseif((ap.eq.2.and.zp.eq.1).and.(at.eq.4.and.zt.eq.2))then
         Rc = 13.5              ! d + 4He

      elseif((ap.eq.2.and.zp.eq.1).and.(zt.eq.6))then
         Rc = 6.0               ! d + C

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.73))then
         Rc = 0.6               ! alpha + Ta

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.79))then
         Rc = 0.6               ! alpha + Au

      else
         Rc = 1.0
      endif

c------------------------------------------------- Rc selection end


c---------------------------------------------- Xm calculation start

      if (ap.eq.1.and.zp.eq.0)then ! neutron + X systems

         SL = 1.2 + 1.6 * ( 1.0 - exp( -E/15. ) )


         if((at.eq.4.and.zt.eq.2) )then
c        n + alpha
            X1 = 5.2

         else
c        others
            X1 = 2.83 - 3.1E-2 * At + 1.7E-4 * At * At

         endif

         Xm = 1. - X1 * exp( -E / ( X1 * SL ) )

      else
         Xm = 1.0
      endif

c---------------------------------------------- Xm calculation end


c------------------------------------------- signe calculation start

      signe = pi * r0 * r0 * ( Ap**(1./3.) + At**(1./3.) + deltaE )**(2.
     $     ) * ( 1 - Rc * B / Ecm ) * Xm

      if(signe.le.0) then
         signe=0
         return
      endif

      signe = signe * 10        ! convert the unit from [fm**2] to [mb]
      signe = signe / 1000      ! convert the unit from [mb] to [b]

      bmax =  r0 * ( Ap**(1./3.) + At**(1./3.) + deltaE )
     $        * sqrt(1. -  Rc * B / Ecm)

c------------------------------------------- signe calculation end

      end subroutine


c----------------------------------------------------------------------
      function radius(a,z)
c----------------------------------------------------------------------
c Purpose   : to obtain the "r_rms,i" values
c References: Atomic Data adn Nuclear Data Tables 36, 495-536 (1987)
c             and
c             NIM-B 152(1999)425-431
c             by H.Iwase Fri Feb  6 2004
c----------------------------------------------------------------------

      implicit none

      real*8 rms(300,2)
      real*8 a, z, radius
      integer i
      integer irnm

c----------------------------------------------------------------------
c     these data are refered from
c     Atomic Data adn Nuclear Data Tables 36, 495-536 (1987)
c     data from page 503 to 510 are used
c----------------------------------------------------------------------

       data irnm / 135 /

       data rms(  1,1), rms(  1,2) /     1, 0.3407/
       data rms(  2,1), rms(  2,2) /  1001, 0.8507/
       data rms(  3,1), rms(  3,2) /  1002, 2.1055/
       data rms(  4,1), rms(  4,2) /  1003, 1.7200/
       data rms(  5,1), rms(  5,2) /  2003, 1.8990/
       data rms(  6,1), rms(  6,2) /  2004, 1.6810/
       data rms(  7,1), rms(  7,2) /  3006, 2.5567/
       data rms(  8,1), rms(  8,2) /  3007, 2.4000/
       data rms(  9,1), rms(  9,2) /  4009, 2.5095/
       data rms( 10,1), rms( 10,2) /  5010, 2.4500/
       data rms( 11,1), rms( 11,2) /  5011, 2.3950/
       data rms( 12,1), rms( 12,2) /  6012, 2.4690/
       data rms( 13,1), rms( 13,2) /  6013, 2.4400/
       data rms( 14,1), rms( 14,2) /  6014, 2.5600/
       data rms( 15,1), rms( 15,2) /  7014, 2.5480/
       data rms( 16,1), rms( 16,2) /  7015, 2.6537/
       data rms( 17,1), rms( 17,2) /  8016, 2.7283/
       data rms( 18,1), rms( 18,2) /  8017, 2.6620/
       data rms( 19,1), rms( 19,2) /  8018, 2.7270/
       data rms( 20,1), rms( 20,2) /  9019, 2.9000/
       data rms( 21,1), rms( 21,2) / 10020, 3.0120/
       data rms( 22,1), rms( 22,2) / 10022, 2.9690/
       data rms( 23,1), rms( 23,2) / 11023, 2.9400/
       data rms( 24,1), rms( 24,2) / 12024, 3.0467/
       data rms( 25,1), rms( 25,2) / 12025, 3.0565/
       data rms( 26,1), rms( 26,2) / 12026, 3.0600/
       data rms( 27,1), rms( 27,2) / 13027, 3.0483/
       data rms( 28,1), rms( 28,2) / 14028, 3.1140/
       data rms( 29,1), rms( 29,2) / 14029, 3.1045/
       data rms( 30,1), rms( 30,2) / 14030, 3.1760/
       data rms( 31,1), rms( 31,2) / 15031, 3.1880/
       data rms( 32,1), rms( 32,2) / 16032, 3.2423/
       data rms( 33,1), rms( 33,2) / 16034, 3.2810/
       data rms( 34,1), rms( 34,2) / 17035, 3.3880/
       data rms( 35,1), rms( 35,2) / 16036, 3.2780/
       data rms( 36,1), rms( 36,2) / 18036, 3.3270/
       data rms( 37,1), rms( 37,2) / 17037, 3.3840/
       data rms( 38,1), rms( 38,2) / 19039, 3.4040/
       data rms( 39,1), rms( 39,2) / 18040, 3.4320/
       data rms( 40,1), rms( 40,2) / 20040, 3.4703/
       data rms( 41,1), rms( 41,2) / 20048, 3.4605/
       data rms( 42,1), rms( 42,2) / 22048, 3.6550/
       data rms( 43,1), rms( 43,2) / 22050, 3.5730/
       data rms( 44,1), rms( 44,2) / 24050, 3.6690/
       data rms( 45,1), rms( 45,2) / 23051, 3.5975/
       data rms( 46,1), rms( 46,2) / 24052, 3.6467/
       data rms( 47,1), rms( 47,2) / 24053, 3.7260/
       data rms( 48,1), rms( 48,2) / 24054, 3.7127/
       data rms( 49,1), rms( 49,2) / 26054, 3.6957/
       data rms( 50,1), rms( 50,2) / 25055, 3.6800/
       data rms( 51,1), rms( 51,2) / 26056, 3.7503/
       data rms( 52,1), rms( 52,2) / 26058, 3.7750/
       data rms( 53,1), rms( 53,2) / 28058, 3.7683/
       data rms( 54,1), rms( 54,2) / 27059, 3.8130/
       data rms( 55,1), rms( 55,2) / 28060, 3.7953/
       data rms( 56,1), rms( 56,2) / 28061, 3.8060/
       data rms( 57,1), rms( 57,2) / 28062, 3.8263/
       data rms( 58,1), rms( 58,2) / 29063, 3.9187/
       data rms( 59,1), rms( 59,2) / 28064, 3.8673/
       data rms( 60,1), rms( 60,2) / 30064, 3.9370/
       data rms( 61,1), rms( 61,2) / 29065, 3.9440/
       data rms( 62,1), rms( 62,2) / 30066, 3.9580/
       data rms( 63,1), rms( 63,2) / 30068, 3.9667/
       data rms( 64,1), rms( 64,2) / 30070, 4.0077/
       data rms( 65,1), rms( 65,2) / 32070, 4.0565/
       data rms( 66,1), rms( 66,2) / 32072, 4.0550/
       data rms( 67,1), rms( 67,2) / 32074, 4.0750/
       data rms( 68,1), rms( 68,2) / 32076, 4.0810/
       data rms( 69,1), rms( 69,2) / 38088, 4.2060/
       data rms( 70,1), rms( 70,2) / 39089, 4.2500/
       data rms( 71,1), rms( 71,2) / 40090, 4.2707/
       data rms( 72,1), rms( 72,2) / 40091, 4.3090/
       data rms( 73,1), rms( 73,2) / 40092, 4.2970/
       data rms( 74,1), rms( 74,2) / 42092, 4.3047/
       data rms( 75,1), rms( 75,2) / 41093, 4.3205/
       data rms( 76,1), rms( 76,2) / 40094, 4.3235/
       data rms( 77,1), rms( 77,2) / 42094, 4.3340/
       data rms( 78,1), rms( 78,2) / 40096, 4.3960/
       data rms( 79,1), rms( 79,2) / 42096, 4.3640/
       data rms( 80,1), rms( 80,2) / 42098, 4.3880/
       data rms( 81,1), rms( 81,2) / 42100, 4.4300/
       data rms( 82,1), rms( 82,2) / 46104, 4.4370/
       data rms( 83,1), rms( 83,2) / 46106, 4.4670/
       data rms( 84,1), rms( 84,2) / 46108, 4.5240/
       data rms( 85,1), rms( 85,2) / 46110, 4.5900/
       data rms( 86,1), rms( 86,2) / 48110, 4.5780/
       data rms( 87,1), rms( 87,2) / 48112, 4.6080/
       data rms( 88,1), rms( 88,2) / 50112, 4.6205/
       data rms( 89,1), rms( 89,2) / 48114, 4.6305/
       data rms( 90,1), rms( 90,2) / 50114, 4.6020/
       data rms( 91,1), rms( 91,2) / 49115, 4.6460/
       data rms( 92,1), rms( 92,2) / 48116, 4.6390/
       data rms( 93,1), rms( 93,2) / 50116, 4.6240/
       data rms( 94,1), rms( 94,2) / 50117, 4.6250/
       data rms( 95,1), rms( 95,2) / 50118, 4.6630/
       data rms( 96,1), rms( 96,2) / 50119, 4.6390/
       data rms( 97,1), rms( 97,2) / 50120, 4.6430/
       data rms( 98,1), rms( 98,2) / 50122, 4.6580/
       data rms( 99,1), rms( 99,2) / 51122, 4.6300/
       data rms(100,1), rms(100,2) / 50124, 4.6807/
       data rms(101,1), rms(101,2) / 56138, 4.8360/
       data rms(102,1), rms(102,2) / 57139, 4.8500/
       data rms(103,1), rms(103,2) / 60142, 4.9253/
       data rms(104,1), rms(104,2) / 60144, 4.9260/
       data rms(105,1), rms(105,2) / 62144, 4.9470/
       data rms(106,1), rms(106,2) / 60146, 4.9815/
       data rms(107,1), rms(107,2) / 60148, 5.0020/
       data rms(108,1), rms(108,2) / 62148, 4.9890/
       data rms(109,1), rms(109,2) / 60150, 5.0037/
       data rms(110,1), rms(110,2) / 62150, 5.0450/
       data rms(111,1), rms(111,2) / 62152, 5.0947/
       data rms(112,1), rms(112,2) / 62154, 5.1259/
       data rms(113,1), rms(113,2) / 64154, 5.1240/
       data rms(114,1), rms(114,2) / 64156, 5.0680/
       data rms(115,1), rms(115,2) / 64158, 5.1720/
       data rms(116,1), rms(116,2) / 67165, 5.2100/
       data rms(117,1), rms(117,2) / 68166, 5.2593/
       data rms(118,1), rms(118,2) / 70174, 5.4100/
       data rms(119,1), rms(119,2) / 71175, 5.3700/
       data rms(120,1), rms(120,2) / 70176, 5.3790/
       data rms(121,1), rms(121,2) / 73181, 5.4800/
       data rms(122,1), rms(122,2) / 74184, 5.4200/
       data rms(123,1), rms(123,2) / 74186, 5.4000/
       data rms(124,1), rms(124,2) / 76192, 5.4130/
       data rms(125,1), rms(125,2) / 78196, 5.3800/
       data rms(126,1), rms(126,2) / 79197, 5.3000/
       data rms(127,1), rms(127,2) / 81203, 5.4630/
       data rms(128,1), rms(128,2) / 82204, 5.4790/
       data rms(129,1), rms(129,2) / 81205, 5.4745/
       data rms(130,1), rms(130,2) / 82206, 5.4963/
       data rms(131,1), rms(131,2) / 82207, 5.5050/
       data rms(132,1), rms(132,2) / 82208, 5.5016/
       data rms(133,1), rms(133,2) / 83209, 5.5167/
       data rms(134,1), rms(134,2) / 90232, 5.7087/
       data rms(135,1), rms(135,2) / 92238, 5.8470/

c----------------------------------------------------------------------

      do i = 1, irnm

         if( nint( z*1000 + a ) .eq. nint( rms(i,1) ) ) then

            radius = rms(i,2)
            return

         end if

      end do

c----------------------------------------------------------------------
c     if there is no data, a fitted values is used.
c     the formula is cited from NIM-B 152(1999)425-431
c----------------------------------------------------------------------

      radius = ( 0.84 * a**(1./3.) + 0.55 )

*-----------------------------------------------------------------------

      return
      end function
      include "cnXsecNasa.f"
