      function kftp(kf)
*                                                                      *
*        get ityp from kf code                                         *
*                                                                      *
*     input:                                                           *
*         kf    : kf code                                              *
*                                                                      *
*     output:                                                          *
*        kftp   : ityp particle id                                    *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

         if( kf .eq. 2212 ) then
               kftp = 1
         else if( kf .eq. 2112 ) then
               kftp = 2
         else if( kf .eq. 211 ) then
               kftp = 3
         else if( kf .eq. 111 ) then
               kftp = 4
         else if( kf .eq. -211 ) then
               kftp = 5
         else if( kf .eq. -13 ) then
               kftp = 6
         else if( kf .eq. 13 ) then
               kftp = 7
         else if( kf .eq. 321 ) then
               kftp = 8
         else if( kf .eq. 311 ) then
               kftp = 9
         else if( kf .eq. -321 ) then
               kftp = 10
         else if( kf .eq. 11 ) then
               kftp = 12
         else if( kf .eq. -11 ) then
               kftp = 13
         else if( kf .eq. 22 ) then
               kftp = 14
         else if( kf .eq. 1000002 ) then
               kftp = 15
         else if( kf .eq. 1000003 ) then
               kftp = 16
         else if( kf .eq. 2000003 ) then
               kftp = 17
         else if( kf .eq. 2000004 ) then
               kftp = 18
         else if( abs(kf) .gt. 1000000 ) then
               kftp = 19
         else
               kftp = 11
         end if

*-----------------------------------------------------------------------

      return
       end function

      function ibryf(ityp,ktyp)
*                                                                      *
*        get baryon number of particles                                *
*                                                                      *
*     input:                                                           *
*        ityp   : particle id                                          *
*        ktyp   : kf code of particle                                  *
*                                                                      *
*     output:                                                          *
*       ibryf   : baryon number of particle                            *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)
!    next is from original Jam; only kchg is used
#include "jam2.inc"   

*-----------------------------------------------------------------------

            if( ityp .eq. 1 .or. ityp .eq. 2 ) then

                  ibryf = 1

            else if( ityp .eq. 11 ) then

                  ibryf = kchg(jamcomp(ktyp),6) / 3
     &                  * isign(1,ktyp)

            else if( ityp .eq. 15 ) then

                  ibryf = 2

            else if( ityp .eq. 16 .or. ityp .eq. 17 ) then

                  ibryf = 3

            else if( ityp .eq. 18 ) then

                  ibryf = 4

            else if( ityp .eq. 19 .and. ktyp .gt. 0 ) then

                  ibryf = ktyp - ktyp / 1000000 * 1000000

            else

                  ibryf = 0

            end if

*-----------------------------------------------------------------------

      return
       end function


************************************************************************
*                                                                      *
      function ichgf(ityp,ktyp)
*                                                                      *
*        get charge of particles                                       *
*                                                                      *
*     input:                                                           *
*        ityp   : particle id                                          *
*        ktyp   : kf code of particle                                  *
*                                                                      *
*     output:                                                          *
*       ichgf   : charge of particle                                   *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

#include "jam2.inc"

*-----------------------------------------------------------------------

            if( ityp .eq. 2 .or. ityp .eq. 4 .or. ityp .eq. 9 .or.
     &          ityp .eq. 14 ) then

                  ichgf = 0

            else if( ityp .eq.  1 .or. ityp .eq. 3 .or.
     &               ityp .eq.  6 .or. ityp .eq. 8 .or.
     &               ityp .eq. 13 ) then

                  ichgf = 1

            else if( ityp .eq.  5 .or. ityp .eq.  7 .or.
     &               ityp .eq. 10 .or. ityp .eq. 12 ) then

                  ichgf = -1

            else if( ityp .eq. 11 ) then

                  ichgf = kchg(jamcomp(ktyp),1) / 3
     &                  * isign(1,ktyp)

            else if( ityp .eq. 15 .or. ityp .eq. 16 ) then

                  ichgf = 1

            else if( ityp .eq. 17 .or. ityp .eq. 18 ) then

                  ichgf = 2

            else if( ityp .eq. 19 ) then

                  ichgf = abs(ktyp) / 1000000

            else

                  ichgf = 0

            end if

*-----------------------------------------------------------------------

      return
       end function
