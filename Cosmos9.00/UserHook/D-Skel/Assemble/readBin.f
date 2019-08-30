      include '../FleshHist/Zprivate1.h'
      character*120 file
      real*8  Et, wx, wy, wz
      integer EventNo
      integer*2  code

      read(*,'(a)') file
      write(0,*) ' file name is '
      write(0,*) file
      open(fnodat, file=file, form="unformatted")
      read(fnodat)
     *        EventNo, code,   Et,  wx, wy, wz
      write(*,'("i ",  i3,  i4, g13.4,3f11.7)')
     *        EventNo, code,   Et,  wx, wy, wz
      
      do while(.true.)
         read(fnodat,end=100) bufc, buf
         do i = 1, bufc
#if KeepWeight != yes
           write(*,
     *      '(6i3, 1pE11.3, 0p,f6.1,1p2E11.3,0p, 2f8.4,f10.6)')
     *       buf(i).ldep,  buf(i).code,  buf(i).subcode,
     *       buf(i).charge, buf(i).ridx, buf(i).faiidx,
     *       buf(i).rinmu, buf(i).fai, buf(i).Ek,
     *       buf(i).t, buf(i).wx, buf(i).wy, buf(i).wz
#else
           write(*,
     *    '(6i3, 1pE11.3, 0p,f6.1,1p2E11.3,0p, 2f8.4,f10.6,1pE11.3')
     *       buf(i).ldep,  buf(i).code,  buf(i).subcode,
     *       buf(i).charge, buf(i).ridx, buf(i).faiidx,
     *       buf(i).rinmu, buf(i).fai, buf(i).Ek,
     *       buf(i).t, buf(i).wx, buf(i).wy, buf(i).wz,
     *       buf(i).wgt
#endif

         enddo
      enddo
 100  continue
      end
