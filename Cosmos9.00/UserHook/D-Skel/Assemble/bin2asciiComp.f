      implicit none
!       This progaram reads binary histogram file and ascii
!       hybrid data  and output it to stdout.  
!       For hist file id is modified to be consistent with
!       each assembled shower. (This is the difference from 
!       bin2ascii.f
      integer ndepth
      parameter (ndepth= 30)
      integer fn1
      real*8  ASdep(ndepth),  munit(ndepth)
      real*8  Esize0(ndepth),
     *     age0(ndepth),  cogdep0(ndepth),
     *     SEloss0(ndepth),
     *     Ng0(ndepth), Ne0(ndepth), Nmu0(ndepth),
     *     cog0
      integer EvNo0

      common /Zbin2ascii/ 
     *    ASDep, Esize0, age0,
     *    cogdep0,  SEloss0, munit,
     *    Ng0, Ne0, Nmu0, cog0,
     *    fn1, EvNo0


      include "../../Hist/Z90histc.f"
      include "../../Hist/Z90hist.f"
      include "../../Hist/Z90hist2.f"
      include "../../Hist/Z90hist3.f"
      type(histogram1) h10
      type(histogram2) h20
      type(histogram3) h30


      integer fn0
      integer kgetenv2
      integer EventNo
      integer leng, i, idx
      integer icon0
      character*120 hist0, hyb0
      character*6 histid0
      

      fn0 = 2
      fn1 = 3
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

      leng = kgetenv2("HYBFILE0", hyb0)
      call copenfw2(fn1, hyb0, 1, icon0)
      if(icon0 .ne. 1)  then
         write(0,*) hyb0(1:leng)
         if( icon0 .eq. 0) then
            write(0,*) 'not exists'
         else
            write(0,*) ' cannot be opened '
         endif
         write(0,*) ' icon=',icon0
         stop 9999
      else
         write(0,*)  hyb0(1:leng), ' opened'
      endif

      EvNo0 =0 
      do while(.true.)
         read( fn0, end=1000 ) histid0
 100     continue
         if( histid0 .eq. '#hist1' ) then
           call kwhistr(h10, fn0, icon0)
           read( h10%id(1:21), '(3x, i4, i3)') EventNo,idx
           do while (EventNo .ne. EvNo0) 
              call get1hyb
           enddo
           write( h10%id(1:21),
     *      '(" # ", i4, i3, f5.2, f5.2)') EventNo,
     *      idx, age0(idx), ASdep(idx)/cog0
           call kwhists(h10, 0.0)
           call kwhistpr(h10, -6) !  to stdout
           call kwhistd(h10)
        elseif(histid0 .eq. '#hist2' ) then
           call kwhistr2(h20, fn0, icon0)

           read( h20%id(1:21), '(3x, i4, i3)') EventNo, idx
           do while (EventNo .ne. EvNo0) 
              call get1hyb
           enddo
           write( h20%id(1:21), '(" # ", i4, i3, f5.2, f5.2)') EventNo,
     *      idx, age0(idx), ASdep(idx)/cog0

           call kwhists2(h20, 0.0)
           call kwhistpr2(h20, -6)
           call kwhistd2(h20)
        elseif(histid0 .eq. '#hist3' ) then
           call kwhistr3(h30, fn0, icon0)

           read( h30%id(1:21), '(3x, i4, i3)') EventNo, idx
           do while (EventNo .ne. EvNo0) 
              call get1hyb
           enddo
           write( h30%id(1:21), '(" # ", i4, i3, f5.2, f5.2)') EventNo,
     *      idx, age0(idx), ASdep(idx)/cog0

           call kwhists3(h30, 0.0)
           call kwhistpr3(h30, -6)
           call kwhistd3(h30)
        else
           write(0,*) 'histid=', histid0, ' invalid'
           stop 9000
        endif
      enddo
 1000 continue
      write(0,*) 'all events processed '
      end
!     **********************
      subroutine get1hyb
      implicit none
      character*120 input0
      integer i, klena

      integer ndepth
      parameter (ndepth= 30)
      integer fn1
      real*8  ASdep(ndepth),  munit(ndepth)
      real*8  Esize0(ndepth),
     *     age0(ndepth),  cogdep0(ndepth),
     *     SEloss0(ndepth),
     *     Ng0(ndepth), Ne0(ndepth), Nmu0(ndepth),
     *     cog0
      integer EvNo0

      common /Zbin2ascii/
     *    ASDep, Esize0, age0,
     *    cogdep0,  SEloss0, munit,
     *    Ng0, Ne0, Nmu0, cog0,
     *    fn1, EvNo0
      input0 = "x" 
      do while (input0 .ne. " ")
         read( fn1 ,'(a)') input0
         if(input0 .ne. " ")  then
            read(input0(1:klena(input0)), *)
     *        EvNo0, i, ASDep(i), Esize0(i), age0(i),
     *        cogdep0(i), SEloss0(i),
     *         munit(i), Ng0(i), Ne0(i), Nmu0(i), cog0
         endif
      enddo
      end

            
