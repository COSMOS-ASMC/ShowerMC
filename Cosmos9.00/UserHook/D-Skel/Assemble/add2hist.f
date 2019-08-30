      implicit none
!   Suppose  we have n histograms
!     h1, h2, h3, ....hn
!   1) cp h1 h0
!   2) h0 + h2 --> h; mv  h h0
!   3) h0 + h3 --> h; mv  h h0
!   ..
!   n) h0 + hn --> h; mv  h h0
! 
!   This program add two histogram; h0 + hx--> h
!            environmental variable
!   file h0: HISTFILE0
!   file hx: HISTFILEX
!   file h:  HISTFILET
!
      include "../../Hist/Z90histc.h"
      include "../../Hist/Z90histo.h"
      include "../../Hist/Z90hist1.h"
      include "../../Hist/Z90hist2.h"
      include "../../Hist/Z90hist3.h"
      type(histogram1) h10, h1x, h1t
      type(histogram2) h20, h2x, h2t
      type(histogram3) h30, h3x, h3t


      integer fn0, fnx, fnt
      integer kgetenv2

      integer leng, i
      integer icon0, iconx, icont
      character*128 hist0, histx, histt
      character*6 histid0, histidx, oldhist
      
      fn0 = 2
      fnx = 3
      fnt = 4
      leng = kgetenv2("HISTFILE0", hist0)
      call copenfw2(fn0, hist0, 2, icon0)
      if(icon0 .ne. 1)  then
         write(0,*) hist0(1:leng)
         if( icon0 .eq. 0) then
            write(0,*) 'not exists'
         else
            write(0,*) ' cannot be opened '
         endif
         write(0,*) ' icon=',icon0
         stop 9999
      else
         write(0,*)  hist0(1:leng), ' opened'
      endif
      leng = kgetenv2("HISTFILEX", histx)
      call copenfw2(fnx, histx, 2, iconx)
      if(iconx .ne. 1)  then
         write(0,*) histx(1:leng)
         if( iconx .eq. 0) then
            write(0,*) 'not exists'
         else
            write(0,*) ' cannot be opened '
         endif
         write(0,*) ' icon=',iconx
         stop 9999
      else
         write(0,*)  histx(1:leng), ' opened'
      endif

      leng = kgetenv2("HISTFILET", histt)
      call copenfw2(fnt, histt, 2, icont)
      if(icont .ne. 0)  then
         write(0,*) histt(1:leng)
         write(0,*) ' cannot be opened '
         write(0,*) ' icon=',icont
         stop 9999
      else
         write(0,*)  histt(1:leng), ' opened'
      endif

      leng = kgetenv2("OLDHIST",  oldhist)
      if(oldhist .eq. "yes") then
         write(0,*) ' old histogram format is assumed'
         call kwhistfmt(.true.) ! old format
      else
         write(0,*) ' new histogram format is assumed'
      endif
      do while(.true.)
         read( fn0, end=1000 ) histid0
         read( fnx ) histidx
         if(histid0 .ne. histidx) then
            write(0,*) histid0, histidx, ' differ'
            stop 9876
         endif
         if( histid0 .eq. '#hist1' ) then
            call kwhistr(h10, fn0, icon0)
            call kwhistr(h1x, fnx, iconx)
            call kwhista(h10, h1x, h1t)

            call kwhistw(h1t, fnt)
            call kwhistd(h10)
            call kwhistd(h1x)
            call kwhistd(h1t)
         elseif(histid0 .eq. '#hist2' ) then
            call kwhistr2(h20, fn0, icon0)
            call kwhistr2(h2x, fnx, iconx)
            call kwhista2(h20, h2x, h2t)
            call kwhistw2(h2t, fnt)
            call kwhistd2(h20)
            call kwhistd2(h2x)
            call kwhistd2(h2t)
         elseif(histid0 .eq. '#hist3' ) then
            call kwhistr3(h30, fn0, icon0)
            call kwhistr3(h3x, fnx, iconx)
            call kwhista3(h30, h3x, h3t)
            call kwhistw3(h3t, fnt)
            call kwhistd3(h30)
            call kwhistd3(h3x)
            call kwhistd3(h3t)
         else
            write(0,*) 'histid=', histid0, ' invalid'
            stop 9000
         endif
      enddo
 1000 continue
      write(0,*) 'all data  processed '
      end


            
