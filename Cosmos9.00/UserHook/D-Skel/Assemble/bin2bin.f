#define  FNODATDEF 33
#include "../FleshHist/Zprivate1.h"
#include "../FleshHist/Zprivate3.h"
!
!       convert .bdat made by FlehsHist into another .bdat
!       with a shorter record length. 
!       input:  
!         number of particles to be put in a record  (< 100000)
!         input filename 
!       echo  '10  "xxx.bdat"' | bin2binPCLinuxIFC 
!      output:  filelname automatically made  to be xxx.bdat2
!
      character*120 file, file2
      real*8  Et, wx, wy, wz
      integer EventNo
      integer*2  code
      integer fnodat2, nptcl, i,j,k,l

      fnodat2=fnodat+1

      read(*, *) nptcl, file
      write(0,*) ' file name is '
      write(0,*) nptcl, file
      file2=trim(file)//"2"
      open(fnodat, file=file, form="unformatted")
      open(fnodat2, file=file2,form="unformatted")
      read(fnodat)
     *        EventNo, code,   Et,  wx, wy, wz
      write(fnodat2)
     *        EventNo, code,   Et,  wx, wy, wz
      
      do while(.true.)
         read(fnodat,end=100) bufc, buf
         do i = 1, bufc, nptcl
            j = min(bufc, i+nptcl-1)
            k = j-i+1
            write(fnodat2)  k, (buf(l), l= i, j)
         enddo
      enddo
 100  continue
      end
