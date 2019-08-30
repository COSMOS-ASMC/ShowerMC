!          don't touch next
#define FNODATDEF 33
!       convert ascii output file made by FleshHist into binary file
!   input; filename is from stdin; let it be, say, xyz.dat 
!   output: fiename for binary is xyz.bdat 
!    echo xyz.dat |  ascii2bin 
#include "../FleshHist/Zprivate1.h"
#include "../FleshHist/Zprivate3.h"
      character*120 file, fileo
      real*8  Et, wx, wy, wz
      character*3 id
      integer EventNo
      integer*2  code
      integer loc, fnodat2
      fnodat2=34

      read(*,'(a)') file
      write(0,*) ' file name is '
      write(0,*) file
      open(fnodat, file=file, form="formatted",action="read")

!      form output format
      loc = index(file, ".dat")
      fileo=file(1:loc)//"bdat"
      write(0,*) ' output filename='
      write(0,*) fileo
      open(fnodat2, file=fileo, form="unformatted")


      read(fnodat,*)
     *      id,  EventNo, code,   Et,  wx, wy, wz
      write(0,*) " id etc=",
     *      id,  EventNo, code,   Et,  wx, wy, wz
      write(fnodat2)
     *        EventNo, code,   Et,  wx, wy, wz
      
      bufc = 0
      do while(.true.)
         i = bufc + 1
#if KeepWeight != yes
         read(fnodat, *, end=100) 
     *        buf(i).ldep,  buf(i).code,  buf(i).subcode,
     *        buf(i).charge, buf(i).ridx, buf(i).faiidx,
     *        buf(i).rinmu, buf(i).fai, buf(i).Ek,
     *        buf(i).t, buf(i).wx, buf(i).wy, buf(i).wz
#else
         read(fnodat, *, end=100) 
     *        buf(i).ldep,  buf(i).code,  buf(i).subcode,
     *        buf(i).charge, buf(i).ridx, buf(i).faiidx,
     *        buf(i).rinmu, buf(i).fai, buf(i).Ek,
     *        buf(i).t, buf(i).wx, buf(i).wy, buf(i).wz,
     *        buf(i).wgt
#endif
         bufc = i
         if(bufc .eq. bufsize) then
            write(fnodat2) bufc, buf
            bufc=0
         endif
      enddo
 100  continue
      if(bufc .gt.  0) then
         write(fnodat2) bufc, buf
         bufc=0
      endif
      end
