!
!         make cnst and  sampling table for Piar
!         into data statement
!
      implicit none
#include "Zmedia.h"
#include "Zmass.h"

       type(epmedia):: media
      integer  io
      character*20 file

      io = 10

      call  cerrorMsg(
     * "media file path( "//
     * "such as '../../Data/Media/Pb')",1)
      file = '../../Data/Media/Pb'
      read(*, *) file
      call cerrorMsg(file,  1)

      open(io, file=file, action='read' ) 
      call epReadTab(io, media)
!        basic consts
      call epwtmedia(media)

!         print brem const; Seltzer region
      call epwtBrCnstS(media%cnst)
!         at low energy
      call epwtBrCnst(media%cnst)
!         other const.
      write(*,*) 'c   consts for pair at complete screening energy'
      write(*,*) '     real*8 cScrMain'
      write(*,*) '     real*8 cScrC1'

      write(*,*) '     data cScrMain/',media%cScrMain, '/'
      write(*,*) '     data cScrC1/',  media%cScrC1, '/' 
!
!         at high energy
      call epwtBrCnstH(media%cnst)
!
      call epwt1dTbl('Brem at Seltzer table region: TX',
     *   'BrTXS', media%tbl%BrTXS, media%cnst%BrTXTS)

      call epwt2dTbl('Brem at Selzer table region; Region A', 
     *  'BrSTSA', media%tbl%BrSTSA, media%cnst%BrUszSA,
     *   media%cnst%BrES)

      call epwt2dTbl('Brem at Seltzer table regtion: Region B', 
     *  'PrSTSB', media%tbl%BrSTSB, media%cnst%BrUszSB,
     *   media%cnst%BrES)

!
      call epwt1dTbl('Brem at low energy: TX',
     *   'BrTXL', media%tbl%BrTXL, media%cnst%BremTXTL)

      call epwt2dTbl('Brem at low energy; Region A', 
     *  'BrSTLA', media%tbl%BrSTLA, media%cnst%BremUszLA,
     *   media%cnst%BremEsize)

      call epwt2dTbl('Brem at low energy: Region B', 
     *  'PrSTLB', media%tbl%BrSTLB, media%cnst%BremUszLB,
     *   media%cnst%BremEsize)



      call epwt1dTbl('Brem at high energy: TX',
     *   'BrTXH', media%tbl%BrTXH, media%cnst%BrneH)

      call epwt2dTbl('Brem at high energy; Region A',
     *  'BrSTHA', media%tbl%BrSTHA, media%cnst%Brnu1H,
     *   media%cnst%BrneH2)

      call epwt2dTbl('Brem at high energy; Region B',
     *  'BrSTHB', media%tbl%BrSTHB, media%cnst%Brnu2H,
     *   media%cnst%BrneH2)

      end
