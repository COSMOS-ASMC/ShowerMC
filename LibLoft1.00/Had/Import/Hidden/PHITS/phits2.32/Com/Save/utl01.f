************************************************************************
*                                                                      *
      subroutine cputime(ic)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to detect cpu time                                      *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision(a-h, o-z)

      parameter ( icpum = 30 )

*-----------------------------------------------------------------------

      common /cputim/ stime(30), cputm(30)
      common /paraj/ mstz(200), parz(200)

*-----------------------------------------------------------------------

      dimension estime(icpum), ecputm(icpum)
      dimension ietm(icpum)

      data  ietm  / icpum * 0 /
      data  stime / icpum * 0.0 /
      data  cputm / icpum * 0.0 /

      data  estime / icpum * 0.0 /
      data  ecputm / icpum * 0.0 /

      save ietm
      save estime, ecputm

*-----------------------------------------------------------------------

      if( mstz(44) .eq. 0 .and. ic .ne. 1 ) return

*-----------------------------------------------------------------------

         if( ic .lt. 1 .or. ic .gt. icpum ) then

            write(*,*) ' **** Error at cputime(ic) : ic is too large'
            write(*,*) ' ==========================================='
            call parastop( 777 )

         end if

*-----------------------------------------------------------------------

         if( ietm(ic) .eq. 0 ) then

            estime(ic) = sect_a(0.0d0)
            ecputm(ic) = cput_a(0.0d0)

            ietm (ic) = 1

         else

            stime(ic) = stime(ic) + sect_a(estime(ic))
            cputm(ic) = cputm(ic) + cput_a(ecputm(ic))

            ietm (ic) = 0

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine timex(time3,ic)
*                                                                      *
*        Last Revised:     2000 06 24                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to get cpu time and print                               *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common/inout/in,io

      data istart / 0 /
      data time0 / 0.0d0 /
      save time0
      save time1

*-----------------------------------------------------------------------

            if( istart .eq. 0 ) then

               time0  = cput_a( time0 )

               istart = 1
               return

            end if

*-----------------------------------------------------------------------

            time2 = cput_a( time0 )
            time3 = time2 - time1
            time1 = time2

            ih = int( time3 / 3600. )
            rm = time3 - dble(ih) * 3600.
            im = int( rm / 60. )
            ts = rm - dble(im) * 60.

      if( ic .eq. 0 ) then

         if( ih .ge. 1 ) then

            write(io,'(''          cpu time = '',
     &                 i3,'' h.'',i3,'' m.'',f6.2,'' s.'')')
     &                 ih, im, ts

         else if( im .ge. 1 ) then

            write(io,'(''          cpu time = '',
     &                 i3,'' m.'',f6.2,'' s.'')')
     &                 im, ts

         else

            write(io,'(''          cpu time = '',
     &                 f7.3,'' s.'')')
     &                 ts

         end if

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function pcmsr(a,b,c)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine the cm momentum from energy and mass       *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              a           : sqrt of s                                 *
*              b           : mass b                                    *
*              c           : mass c                                    *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      if( a .gt. b + c ) then

         pcmsr = sqrt( ( a**2 - ( b + c )**2 )
     &               * ( a**2 - ( b - c )**2 ) )
     &               / ( 2.0 * a )

      else

         pcmsr = 0.0

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine chname(ic,mmas,mchg,chau)
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give Nucleus name as character string                *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              ic          : =0 normal, =1 error                       *
*              mmas        : mass of Nucleus                           *
*              mchg        : charge of Nucleus                         *
*              chau        : output of the name as character           *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      character       chau*8
      character       nucl*3
      character       nmas*3
      character       num(0:9)*1

      character       elmnt(104)*3

*-----------------------------------------------------------------------

      data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

      data num / '0','1','2','3','4','5','6','7','8','9' /

*-----------------------------------------------------------------------
*        nucleon and pion
*-----------------------------------------------------------------------

                  ic = 0

            if( mchg .le. 104 .and.
     &          mchg .gt. 0 .and. mmas .gt. 1 ) then

                  nucl(1:3) = elmnt( mchg )

            else if( mchg .le. 104 .and.
     &               mchg .gt. 0 .and. mmas .eq. 0 ) then

                  nucl(1:3) = elmnt( mchg )

            else

               if( mchg .eq. 1 .and. mmas .eq. 1 ) then

                  chau = 'p'

               else if( mchg .eq. 0 .and. mmas .eq. 1 ) then

                  chau = 'n'

               else if( mchg .eq. 1 .and. mmas .eq. 2 ) then

                  chau = 'd'

               else if( mchg .eq. 1 .and. mmas .eq. -1 ) then

                  chau = 'pi+'

               else if( mchg .eq. 0 .and. mmas .eq. -1 ) then

                  chau = 'pi0'

               else if( mchg .eq. -1 .and. mmas .eq. -1 ) then

                  chau = 'pi-'

               else

                  ic = 1

               end if

                  return

            end if

*-----------------------------------------------------------------------
*        Nucleus
*-----------------------------------------------------------------------

            if( mmas .gt. 0 ) then

                  i1   = mod( mmas,       10 )
                  i10  = mod( mmas / 10,  10 )
                  i100 = mod( mmas / 100, 10 )

               if( mmas .lt. 10 ) then

                  nmas(1:1) = num(i1)
                  lmas      = 1

               else if( mmas .lt. 100 ) then

                  nmas(1:2) = num(i10)//num(i1)
                  lmas      = 2

               else

                  nmas(1:3) = num(i100)//num(i10)//num(i1)
                  lmas      = 3

               end if

                  chau  = nmas(1:lmas)//nucl

            else

                  chau  = nucl

            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine kfcname(inkf,m,sname)
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to transfer kf code to character                        *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      character*(*)   sname

*-----------------------------------------------------------------------

            write(sname,'(i9)') inkf

            i = 0
            k = 0

  100       i = i + 1

  200       if( i .eq. m .or. sname(i:i) .ne. ' ' ) then

               if( k .gt. 0 ) then

                  do l = m, m - k + 1, -1

                     sname(l:l) = ' '

                  end do

               end if

               return

            end if

            if( sname(i:i) .eq. ' ' ) then

                  k = k + 1

               do j = i, m - 1

                  sname(j:j) = sname(j+1:j+1)

               end do

               goto 200

            end if

*-----------------------------------------------------------------------

      return
      end

