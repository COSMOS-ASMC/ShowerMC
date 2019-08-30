************************************************************************
*                                                                      *
      function sect_a(stime)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to get elapse time                                      *
*              use function : secnds(stime)                            *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              stime    : starting time                                *
*              sect_a   : elapse time from starting time               *
*                                                                      *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

      real*4 secnds

*-----------------------------------------------------------------------

            sect_a = dble( secnds( real( stime ) ) )
c           sect_a = dble( secnds( ) )

*-----------------------------------------------------------------------

      return
       end function

************************************************************************
*                                                                      *
      function cput_a(cputm)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to get cpu time                                         *
*              use function : dtime(tarray)                            *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              cputm    : starting time                                *
*              cput_a   : cpu time from starting time                  *
*                                                                      *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

      real*4 tarray
      real*4 dtime

*-----------------------------------------------------------------------

      dimension tarray(2)

      data cpuall / 0.0 /
      save cpuall

*-----------------------------------------------------------------------

            cpuall = cpuall + dble( dtime(tarray) )
            cput_a = cpuall - cputm

*-----------------------------------------------------------------------


      return
       end function


************************************************************************
*                                                                      *
      subroutine date_a_time(iyer,imon,iday,ihor,imin,isec)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to get current date                                     *
*              use subroutine idate and itime                          *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              iyer     : 4-digit year                                 *
*              imon     : month                                        *
*              iday     : day                                          *
*              ihor     : hour                                         *
*              imin     : minite                                       *
*              isec     : secont                                       *
*                                                                      *
*                                                                      *
************************************************************************

      integer   inum(3)

*-----------------------------------------------------------------------
!     gfortran on Mac fail to comile. says too many argmuents
!     don't know on Linux.
!            call idate(inum(1),inum(2),inum(3))
!     put dummy values
      inum(3)=17
      inum(2) = 1
      inum(1) = 1
!      
            if( inum(3) .lt. 50 ) then

                  inum(3) = 2000 + inum(3)

            else if( inum(3) .lt. 1000 ) then

                  inum(3) = 1900 + inum(3)

            end if

            iyer = inum(3)
            imon = inum(1)
            iday = inum(2)

C for SUN ----------------------------------
C           imon = inum(2)
C           iday = inum(1)
C-------------------------------------------

            call itime(inum)

            ihor = inum(1)
            imin = inum(2)
            isec = inum(3)

*-----------------------------------------------------------------------
*     for f90
*-----------------------------------------------------------------------

c     dimension ivalue(8)
c     character*12 dum(3)

*-----------------------------------------------------------------------

c           call date_and_time(dum(1),dum(2),dum(3),ivalue)

c           iyer = ivalue(1)
c           imon = ivalue(2)
c           iday = ivalue(3)
c           ihor = ivalue(5)
c           imin = ivalue(6)
c           isec = ivalue(7)

*-----------------------------------------------------------------------

      return
       end subroutine


