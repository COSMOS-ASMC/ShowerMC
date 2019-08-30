!
!         make cnst and  sampling table for Brem
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
!         general consts for the media
      call epwtmedia(media)

!         print pair const. 
      call epwtPrCnst(media%cnst)

!
      call epwtPrCnstH(media%cnst)
!
      call epwt1dTbl('Pair at Low E: TX', 'PrTXL',
     * media%tbl%PrTXL, media%cnst%PairTXTL)

      call epwt2dTbl('Pair at Low E: Region A', 
     *  'PrSTLA', media%tbl%PrSTLA, media%cnst%PairUszLA,
     *   media%cnst%PairEsize)

      call epwt2dTbl('Pair at Low E: Region B', 
     *  'PrSTLB', media%tbl%PrSTLB, media%cnst%PairUszLB,
     *   media%cnst%PairEsize)

!      call epwt1dTbl('Pair at High E: TX', 'PrTXH',
!     * media.tbl.PrTXH, media.cnst.PrneH)

!      call epwt2dTbl('Pair at High E:',
!     *  'PrSTH', media.tbl.PrSTH, media.cnst.Prnu1H,
!     *   media.cnst.PrneH)

      end
