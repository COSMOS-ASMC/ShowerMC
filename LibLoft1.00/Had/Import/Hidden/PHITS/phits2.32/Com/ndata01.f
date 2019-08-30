************************************************************************
*                                                                      *
      function weitn(nz,nn)
*                                                                      *
*                                                                      *
*        Last Revised:     2002 02 15                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give weight (1.e-24 g)                               *
*              Experimental Data                                       *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              nz          : proton  number                            *
*              nn          : neutron number                            *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      data rpms / 938.27d0 /
      data rnms / 939.58d0 /

*-----------------------------------------------------------------------

      weitn = nz * rpms + nn * rnms - bindeg(nn,nz)

      weitn = weitn * 1.78266270d-3

*-----------------------------------------------------------------------

      return
       end function


************************************************************************
*                                                                      *
      function bindeg(nz,nn)
*                                                                      *
*                                                                      *
*        Last Revised:     2002 02 15                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give binding energy (MeV)                            *
*              Experimental Data                                       *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              nz          : proton  number                            *
*              nn          : neutron number                            *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param01.inc'

*-----------------------------------------------------------------------

      parameter( bvol=15.56, bsur=17.23, bsym=46.57, coul=1.2 )

*-----------------------------------------------------------------------

      common /masbe0/ be0(0:nnuc)
      common /masid0/ inuc0(0:maxpt,0:maxnt)

*-----------------------------------------------------------------------

               bindeg = 0.0

               mm = nz + nn

*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------

            if( nz .lt. 0 .or. nn .lt. 0 .or.
     &          mm .le. 1 ) then

               return

*-----------------------------------------------------------------------
*        Normal Nucleus
*-----------------------------------------------------------------------
                 !!!KK  if check mode, ifort evaluates incu0
  !    even if nz> maxpt and  results in error. so we 
  !    split next
!!!            else if( nz .le. maxpt .and. nn .le. maxnt .and.
!!!     &               inuc0(nz,nn) .ne. 0 ) then
            else if( nz .le. maxpt .and. nn .le. maxnt) then
               if( inuc0(nz,nn) .ne. 0 ) then
                  bindeg = be0( inuc0(nz,nn) )
               endif
*-----------------------------------------------------------------------
*        Super Heavy Nuclei or out of table
*        Use Liquid Drop Mass Formula
*-----------------------------------------------------------------------
            endif  !!! KK

            if( bindeg == 0.) then  !!! KK

               a  = dble( mm )
               rc = 1.24 * a**(1./3.)

               bindeg = ( bvol * a
     &                - bsur * a**(2./3.)
     &                - bsym / 2. * ( nn - nz ) * ( nn - nz ) / a
     &                - 3./5. * nz * coul * nz * coul / rc )

            end if

*-----------------------------------------------------------------------

      return
       end function


************************************************************************
*                                                                      *
      subroutine masstbl
*                                                                      *
*        Last Revised:     2002/02/15                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to read mass table from block data                      *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param01.inc'

*-----------------------------------------------------------------------

      common /masmx0/ imax0

      common /maspn0/ nzz0(0:nnuc), nnn0(0:nnuc), jgs0(0:nnuc)
      common /masbe0/ be0(0:nnuc)
      common /masid0/ inuc0(0:maxpt,0:maxnt)

*-----------------------------------------------------------------------

               do iz = 0, maxpt
               do in = 0, maxnt

                  inuc0(iz,in) = 0

               end do
               end do

*-----------------------------------------------------------------------
*           Normal Nucleus data from block data
*-----------------------------------------------------------------------

               do i = 1, imax0

                  iz = nzz0(i)
                  in = nnn0(i)

                  inuc0(iz,in) = i

               end do

*-----------------------------------------------------------------------

      return
       end subroutine
