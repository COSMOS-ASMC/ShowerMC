      module modZbin2ascii
      implicit none
      integer,parameter:: ndepth=50
      integer fn1
      real*8  ASdep(ndepth),  munit(ndepth)
      real*8  Esize0(ndepth),
     *     age0(ndepth),  cogdep0(ndepth),
     *     SEloss0(ndepth),
     *     Ng0(ndepth), Ne0(ndepth), Nmu0(ndepth),
     *     Nh0(ndepth), 
     *     cog, cog2
      integer EvNo0

      contains
      subroutine get1hyb( rew )
      implicit none
      logical rew
      character*128 input0
      integer i, klena

      integer code, subcode, charge

      real  E0, w1, w2, w3, firstz

!        data format
!        1   xxx 
!        1   xxx       
!        1   xxx 
!        1   xxx       
!      
!        2   xxx 
!        2   xxx       
!        2   xxx 
!        2   xxx       
!      
!        3   xxx 
!        3   xxx       
!        3   xxx 
!        3   xxx 
      if(rew) rewind fn1
      input0 = "x" 
      read(fn1,'(a)') input0
      read(input0(3:klena(input0)), *) 
     * EvNo0, code, subcode, charge,
     *    E0, w1, w2, w3, firstz, cog, cog2
      do while (input0(1:10) .ne. " ")
         input0=" "
         read( fn1 ,'(a)') input0
         if(input0(1:10) .ne. " ")  then
            read(input0(1:klena(input0)), *)
     *       i, ASDep(i),  munit(i),    age0(i),   cogdep0(i), 
     *       Ng0(i), Ne0(i), Nmu0(i),  Nh0(i),
     *       Esize0(i), SEloss0(i)
         endif
      enddo
      end subroutine get1hyb
!     ***********************
      subroutine  mergehyb1(h1)
      use modHistogram1
      use modHistogram2
      use modHistogram3
      implicit none

      type(histogram1) h1
      type(histogram2) h2
      type(histogram3) h3

      integer klena
      integer j

      do while (h1%c%eventno .ne. EvNo0) 
         call get1hyb( h1%c%eventno .lt. EvNo0) 
      enddo
!
!        this part must be consistent with
!        FleshHist/interface.f output  for evid
      read(h1%c%id, '(i3)') j
      write(h1%c%id, 
     *   '(i3, i5,  f5.2, f5.2,
     *   i5,  i4)')
     *   j, int( ASDep(j) ),
     *   age0(j), cogdep0(j),
     *   int(munit(j)), int(cog2)
                                          
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
     *   i5,  i4)')
     *   j, int( ASDep(j) ),
     *   age0(j), cogdep0(j),
     *   int(munit(j)), int(cog2)

      return
!     *****************      
      entry mergehyb3(h3)
!     ****************
            

      do while (h3%c%eventno  .ne. EvNo0) 
         call get1hyb( h3%c%eventno  .lt. EvNo0) 
      enddo

      read(h3%c%id, '(i3)') j
      write(h3%c%id, 
     *   '(i3, i5,  f5.2, f5.2,
     *   i5,  i4)')
     *   j, int( ASDep(j) ),
     *   age0(j), cogdep0(j),
     *   int(munit(j)), int(cog2)


      end subroutine mergehyb1


      subroutine openhyb(icon)
      implicit none
      integer icon  ! output.  1--> hybrid must be read
                    !          0--> hybrid need not be used
      integer leng

      character*120  hyb0
      integer kgetenv2
 
      
      fn1= 3
      leng = kgetenv2("HYBFILE0", hyb0)
      call copenfw2(fn1, hyb0, 1, icon)
      if(icon .ne. 1)  then
         write(0,*)
     *    '*************** caution ************'
         write(0,*)
     *    "You haven't given env. var. HYBFILE0"
         write(0,*) 
     *    "or File specified by HYBFILE0"
         if( icon .eq. 0) then
            write(0,*) 'not exists'
         else
            write(0,*) ' cannot be opened '
         endif
         write(0,*) 
     *   "It's ok if you don't merge hybrid data file"
         icon = 0
      else
         write(0,*)  hyb0(1:leng), ' opened'
         icon = 1
      endif

      EvNo0 =0 
      end subroutine openhyb
      end module modZbin2ascii
    
