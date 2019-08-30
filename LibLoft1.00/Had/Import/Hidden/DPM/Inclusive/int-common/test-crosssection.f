      program test_crs_engel
      implicit doubleprecision(a-h,o-z)
      character line*80
      common /crosssection_factor/fact_cross
      common /prod_ratio/ prodratio, crs_e
      common /he_factor/fac_alpha
      fac_alpha = 1.77d0

c      if(iargc().lt.1) then
c         stop 'test-crs kin'
c      end if
c      call getarg(1, line)
c      read(line, *) kin

      fact_cross = 1
      do i=1, 7*20+1
         ek = 10.d0**(real(i-1)/20 - 1.d0)
         crs_inel = crossint(10, ek)
         crs_prd = crs_inel * prodratio
         
         crs_he = crossint(21, 4*ek)
         write(*,'(1p,5E13.5)') ek,crs_inel,crs_prd,crs_e
c         write(*,'(1p,5E13.5)') ek,crs_add(ek),crs_prd,crs_e
c         write(*,*) kin, ek, crs
      enddo
      end

      function crs_add(ek)
      implicit doubleprecision (a-h,o-z)
      elk = log10(ek)

c      crs_add = 305.*(1.-1.2*(elk-0.1)**2)
      crs_add = 310.*(1.-0.9*(elk-0.15)**2)
      end

c         crs_he = crossint(21, 4*ek)
c         write(*,'(1p,5E13.5)') ek, crs_inel, crs_prd, crs_e
c         write(*,*) kin, ek, crs
c      enddo
c      end
