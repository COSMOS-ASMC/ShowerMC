      program main
      use modHistogram1
      use modHistogram2
      use modHistogram3
      implicit none
!       binary hist file:  binary output file made by programs
!       in UserHook/Hist
!       This program is intended to read binary hist file 
!       and hybrid data file to modify the hybrid data information
!       in the binary hist file, and then re-make the binary
!       hist file.
!  Usage: compile:  make -f bin2bin.mk
!         execution:
!            set environmental variable HISTFILE0 to be
!            a binary hist file.
!            you MUST give another env. var.
!            HYBFILE0 to be a hybrid data output file
!            made by a program in UserHook/DisPara/FleshHist.
!        bin2bin$ARCH 
!            temp.hist will be created.
!

      type(histogram1) h10
      type(histogram2) h20
      type(histogram3) h30


      integer kgetenv2
      integer leng, i, fn0, fout
      integer icon0, iconhyb
      character*120 hist0, histout
      character*6 histid0, oldhist
      real normf
      data normf/-1.0/

      fn0 =  2
!     fn1 =  3 may be defined later for hybrid data
      fout = 4
      call kwhistso( 2 ) ! binary write 
      leng = kgetenv2("HISTFILE0", hist0)
      call copenfw2(fn0, hist0, 2, icon0)
      if(icon0 .ne. 1)  then
         write(0,*) "File specified by HISTFILE0 "
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
      leng = kgetenv2("HISTFILE1", histout)
      if(leng .le. 0) then
         write(0,*)
     *       "Output file name is not given by the HISTFILE1"
         write(0,*) "Set Env. Var. HISTFILE1"
         stop 1234
      endif
      call copenfw2(fout, histout(1:leng), 2, icon0)
      if(icon0 .ne. 0)  then
         write(0,*) "For binary output "//histout(1:leng)
         write(0,*) "is specified but old one seems to exist"
         write(0,*) "delete or mv that file beforehand"
        stop 9999
      else
         write(0,*)
     *     histout(1:leng)//' will be created for binary hist"'
      endif

!c?      call openhyb(iconhyb)
      leng = kgetenv2("OLDHIST", oldhist)
      if( leng .gt. 0  .and. oldhist .eq. "yes") then
         write(0,*) 'Old hist format assumed'
         call kwhistfmt(.true.)
      else
         call kwhistfmt(.false.)
      endif

      do while(.true.)
         read( fn0, end=1000 ) histid0
 100     continue
         if( histid0 .eq. '#hist1' ) then
           call kwhistr(h10, fn0, icon0)

           if(iconhyb .eq. 1)  then
!              call mergehyb1(h10)
           endif
           call kwhists(h10, normf)
           call kwhistp(h10,  fout) 
           call kwhistd(h10)
        elseif(histid0 .eq. '#hist2' ) then
           call kwhistr2(h20, fn0, icon0)
           if(iconhyb .eq. 1) then
              call mergehyb2(h20)
           endif
           call kwhists2(h20, normf)
           call kwhistp2(h20, fout)
           call kwhistd2(h20)
        elseif(histid0 .eq. '#hist3' ) then
           call kwhistr3(h30, fn0, icon0)
           if(iconhyb .eq. 1) then
              call mergehyb3(h30)
           endif
           call kwhists3(h30, normf)
           call kwhistp3(h30, fout)
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
      subroutine get1hyb( rew ) 
      implicit none
      logical rew
      character*128 input0
      integer i, klena

      integer ndepth
      parameter (ndepth= 50)
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

      if(rew) rewind fn1

      input0 = "x" 
!////////////
!      write(0,*) ' while'
!////////
      do while (input0(1:10) .ne. " ")
          input0=" "
          read( fn1 ,'(a)') input0
!////////////
!          write(0,*) ' input0=',input0
!//////////////
          if(input0(1:10) .ne. " ")  then
            read(input0(1:klena(input0)), *)
     *        EvNo0, i, ASDep(i), Esize0(i), age0(i),
     *        cogdep0(i), SEloss0(i),
     *         munit(i), Ng0(i), Ne0(i), Nmu0(i), cog0
!/////////
!          write(0,*) ' input0 read'
!//////////////
          endif
      enddo
      end
!     ***********************
      subroutine  mergehyb1(h1)
      use modHistogram1
      use modHistogram2
      use modHistogram3
      implicit none

      type(histogram1) h1
      type(histogram2) h2
      type(histogram3) h3

      integer ndepth
      integer nc
      parameter (ndepth= 50)
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

      integer klena
      integer j

      do while (h1%c%eventno .ne. EvNo0) 
         call get1hyb( h1%c%eventno .lt. EvNo0) 
      enddo
!
!       available variables
!     *      idx, ASdep(idx), Esize0(idx), age0(idx), 
!     *      SEloss0(idx), munit(idx), 
!     *      Ng0(idx), Ne0(idx), Nmu0(idx), 
!     *      ASdep(idx)/cog0, cog0
!
!        this part must be consistent with
!        FleshHist/interface.f output  for evid
!//////////////
!      write(0,*) ' id', h1%c%id
!/////////// 

      read(h1%c%id, '(i3)') j
!//////////////
!      write(0,*) ' j=',j
!///////////
      write(h1%c%id, 
     *   '(i3, i5,  f5.2, f5.2,
     *   i5,  i5)')
     *   j, int( ASDep(j) ),
     *   age0(j), ASDep(j)/cog0,
     *   int(munit(j)), int(cog0)
!//////////////
!      write(0,*) ' j=',j
!///////////
      return
!     *******************
      entry  mergehyb2(h2)
!     *******************
            

      do while (h2%c%eventno .ne. EvNo0) 
         call get1hyb( h2%c%eventno .lt. EvNo0) 
      enddo
      read(h2%c%id, '(i3)') j
      write(h2%c%id, 
     *   '(i3, i5,  f5.2, f5.2,
     *   i5,  i5)')
     *   j, int( ASDep(j) ),
     *   age0(j), ASDep(j)/cog0,
     *   int(munit(j)), int(cog0)

      return
!     *****************      
      entry mergehyb3(h3)
!     ****************
      do while (h3%c%eventno  .ne. EvNo0) 
         call get1hyb( h3%c%eventno .lt. EvNo0) 
      enddo

      read(h3%c%id, '(i3)') j
      write(h3%c%id, 
     *   '(i3, i5,  f5.2, f5.2,
     *   i5,  i5)')
     *   j, int( ASDep(j) ),
     *   age0(j), ASDep(j)/cog0,
     *   int(munit(j)), int(cog0)
      end


      subroutine openhyb(icon)
      implicit none
      integer icon  ! output.  1--> hybrid must be read
                    !          0--> hybrid need not be used
      integer leng
      integer ndepth
      parameter (ndepth= 50)
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

      character*120  hyb0
      integer kgetenv2
      
      fn1= 3
      leng = kgetenv2("HYBFILE0", hyb0)
      call copenfw2(fn1, hyb0, 1, icon)
      if(icon .ne. 1)  then
         write(0,*)
     *    "You haven't given env. var. HYBFILE0"
         write(0,*) 
     *    "or File specified by HYBFILE0"
         if( icon .eq. 0) then
            write(0,*) 'not exists'
         else
            write(0,*) ' cannot be opened '
         endif
         stop 9999
      else
         write(0,*)  hyb0(1:leng), ' opened'
         icon = 1
      endif

      EvNo0 =0 
      end
