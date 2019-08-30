************************************************************************
*                                                                      *
      function rmtyp(ityp,ktyp)
*                                                                      *
*        get mass (MeV) of particles                                   *
*                                                                      *
*     input:                                                           *
*        ityp   : particle id                                          *
*        ktyp   : kf code of particle                                  *
*                                                                      *
*     output:                                                          *
*       rmtyp   : mas (MeV) of particle                                *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

      include 'jam2.inc'
      parameter ( rpmass = 938.27, rnmass = 939.58 )

*-----------------------------------------------------------------------

      dimension rmas(20)

      data rmas/938.27,939.58,139.6,135.0,139.6,105.7,105.7,
     &          493.6,497.7,493.6,0.0,0.511,0.511,
     &          7*0.0/

*-----------------------------------------------------------------------

                  rmtyp = rmas(ityp)

            if( ityp .eq. 11 ) then

                  kc    = jamcomp(ktyp)
                  rmtyp = pmas(kc,1) * 1000.0

            else if( ityp .eq. 15 ) then

                  iz = 1
                  in = 1
                  be = bindeg(iz,in)

                  rmtyp = iz * rpmass + in * rnmass - be

            else if( ityp .eq. 16 ) then

                  iz = 1
                  in = 2
                  be = bindeg(iz,in)

                  rmtyp = iz * rpmass + in * rnmass - be

            else if( ityp .eq. 17 ) then

                  iz = 2
                  in = 1
                  be = bindeg(iz,in)

                  rmtyp = iz * rpmass + in * rnmass - be

            else if( ityp .eq. 18 ) then

                  iz = 2
                  in = 2
                  be = bindeg(iz,in)

                  rmtyp = iz * rpmass + in * rnmass - be

            else if( ityp .eq. 19 .and. ktyp .gt. 0 ) then

                  iz = ktyp / 1000000
                  in = ktyp - ktyp / 1000000 * 1000000 - iz
                  be = bindeg(iz,in)

                  rmtyp = iz * rpmass + in * rnmass - be

            end if

*-----------------------------------------------------------------------

      return
       end function


************************************************************************
*                                                                      *
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

      include 'jam2.inc'

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

      include 'jam2.inc'

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


************************************************************************
*                                                                      *
      function kfft(ityp1)
*                                                                      *
*        get kf code from ityp                                         *
*                                                                      *
*     input:                                                           *
*        ityp   : particle id                                          *
*                                                                      *
*     output:                                                          *
*       kfft    : kf code                                              *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

            if( ityp1 .eq. 1 ) then
               kf1 = 2212
            else if( ityp1 .eq. 2 ) then
               kf1 = 2112
            else if( ityp1 .eq. 3 ) then
               kf1 = 211
            else if( ityp1 .eq. 4 ) then
               kf1 = 111
            else if( ityp1 .eq. 5 ) then
               kf1 = -211
            else if( ityp1 .eq. 6 ) then
               kf1 = -13
            else if( ityp1 .eq. 7 ) then
               kf1 =  13
            else if( ityp1 .eq. 8 ) then
               kf1 =  321
            else if( ityp1 .eq. 9 ) then
               kf1 =  311
            else if( ityp1 .eq. 10 ) then
               kf1 = -321
            else if( ityp1 .eq. 12 ) then
               kf1 =  11
            else if( ityp1 .eq. 13 ) then
               kf1 =  -11
            else if( ityp1 .eq. 14 ) then
               kf1 =  22
            else if( ityp1 .eq. 15 ) then
               kf1 =  1000002
            else if( ityp1 .eq. 16 ) then
               kf1 =  1000003
            else if( ityp1 .eq. 17 ) then
               kf1 =  2000003
            else if( ityp1 .eq. 18 ) then
               kf1 =  2000004
            else
               kf1 = 0
            end if

               kfft = kf1

*-----------------------------------------------------------------------

      return
       end function


************************************************************************
*                                                                      *
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


************************************************************************
*                                                                      *
      function ipatf(ityp,ktyp)
*                                                                      *
*        get iclusts code from particle id and kf code                 *
*                                                                      *
*     input:                                                           *
*        ityp   : particle id                                          *
*        ktyp   : kf code of particle                                  *
*                                                                      *
*     output:                                                          *
*       ipatf   : iclust code                                          *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

            if( ityp .eq. 0 ) then
               ipatf = 0
            else if( ityp .ge. 15 .and. ityp .le. 18 ) then
               ipatf = 0
            else if( ityp .eq. 1 ) then
               ipatf = 1
            else if( ityp .eq. 2 ) then
               ipatf = 2
            else if( ityp .ge. 3 .and. ityp .le. 5 ) then
               ipatf = 3
            else if( ityp .ge. 6 .and. ityp .le. 7 ) then
               ipatf = 6
            else if( ityp .ge. 8 .and. ityp .le. 10 ) then
               ipatf = 5
            else if( ityp .eq. 11 .and. ktyp .eq. 22 ) then
               ipatf = 4
            else if( ityp .eq. 11 ) then
               ipatf = 7
            else if( ityp .ge. 12 .and. ityp .le. 13 ) then
               ipatf = 7
            else if( ityp .eq. 14 ) then
               ipatf = 4
            else
               ipatf = 0
            end if

*-----------------------------------------------------------------------

      return
       end function


************************************************************************
*                                                                      *
      subroutine iptch(ipt,ityp,ktyp,ipat,ichg,ibry,iprt,inut,rmss)
*                                                                      *
*        get iclusts code, ityp, ktyp, ipatf, ichgf, ibryf             *
*            from particle id of MCNPX                                 *
*                                                                      *
*     input:                                                           *
*        ipt    : particle id of MCNPX                                 *
*                                                                      *
*     output:                                                          *
*        ityp   : particle id                                          *
*        ktyp   : kf code of particle                                  *
*        ipat   : iclust code                                          *
*        ichg   : charge                                               *
*        ibry   : barion number                                        *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

               ityp = -1

         if( ipt .eq. 1 ) then

               ityp = 2
               ktyp = 2112
               ipat = 2
               ichg = 0
               ibry = 1
               iprt = 0
               inut = 1
               rmss = rmtyp(ityp,ktyp)

         else if( ipt .eq. 2 ) then

               ityp = 14
               ktyp = 22
               ipat = 4
               ichg = 0
               ibry = 0
               iprt = 0
               inut = 0
               rmss = rmtyp(ityp,ktyp)

         else if( ipt .eq. 3 ) then

            if( unirn(dummy) .gt. 0.5 ) then

               ityp = 12
               ktyp = 11
               ipat = 7
               ichg = -1
               ibry = 0
               iprt = 0
               inut = 0
               rmss = rmtyp(ityp,ktyp)

            else

               ityp = 13
               ktyp = -11
               ipat = 7
               ichg = 1
               ibry = 0
               iprt = 0
               inut = 0
               rmss = rmtyp(ityp,ktyp)

            end if

         else if( ipt .eq. 9 ) then

               ityp = 1
               ktyp = 2212
               ipat = 1
               ichg = 1
               ibry = 1
               iprt = 1
               inut = 0
               rmss = rmtyp(ityp,ktyp)

         else if( ipt .eq. 20 ) then

            if( unirn(dummy) .gt. 0.5 ) then

               ityp = 3
               ktyp = 211
               ipat = 3
               ichg = 1
               ibry = 0
               iprt = 0
               inut = 0
               rmss = rmtyp(ityp,ktyp)

            else

               ityp = 5
               ktyp = -211
               ipat = 3
               ichg = 1
               ibry = 0
               iprt = 0
               inut = 0
               rmss = rmtyp(ityp,ktyp)

            end if

         else if( ipt .eq. 21 ) then

               ityp = 4
               ktyp = 111
               ipat = 3
               ichg = 0
               ibry = 0
               iprt = 0
               inut = 0
               rmss = rmtyp(ityp,ktyp)

         else if( ipt .eq. 31 ) then

               ityp = 15
               ktyp = 1000002
               ipat = 0
               ichg = 1
               ibry = 2
               iprt = 1
               inut = 1
               rmss = rmtyp(ityp,ktyp)

         else if( ipt .eq. 32 ) then

               ityp = 16
               ktyp = 1000003
               ipat = 0
               ichg = 1
               ibry = 3
               iprt = 1
               inut = 2
               rmss = rmtyp(ityp,ktyp)

         else if( ipt .eq. 33 ) then

               ityp = 17
               ktyp = 2000003
               ipat = 0
               ichg = 2
               ibry = 3
               iprt = 2
               inut = 1
               rmss = rmtyp(ityp,ktyp)

         else if( ipt .eq. 34 ) then

               ityp = 18
               ktyp = 2000004
               ipat = 0
               ichg = 2
               ibry = 4
               iprt = 2
               inut = 2
               rmss = rmtyp(ityp,ktyp)

         end if

*-----------------------------------------------------------------------

      return
       end subroutine


