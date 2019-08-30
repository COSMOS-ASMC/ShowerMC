      integer leng,  lengn, lengid
      integer lengdir
      character*120 dir
      integer kgetenv2
      character*8 numb
      character*64 execid
      character*128 msg
      logical takehist
      save

      do i= 1, nsites
         histdep(i) = 0
         indivdep(i) = 0
      enddo

      

      leng = kgetenv2("INDIVDEP", input)
      read(input, *) indivdep

      leng = kgetenv2("OUTPUT", input)
      read(input, *)   tkrtspec, tkarspec

      leng = kgetenv2("HISTDEP", input)
      read(input, *) histdep
      takehist = .false.
      do i= 1, nsites
         if(histdep(i) .gt. 0 ) then
            takehist=.true.
            exit
         endif
      enddo
      takehist = takehist .and. (tkrtspec .or. tkarspec)

      leng = kgetenv2("BINW", input)
      read(input, *)   binw



      input = ' '
#if defined (KEKB) || (KEKA)
      numb="rank-#"
      lengn=6
#else
      lengn =  kgetenv2("NUMB", numb)
#endif
      lengdir =  kgetenv2("OUTDIR", input)
      dir = input
      lengid = kgetenv2("EXECID", execid)


!      msg =input(1:lengdir)//"/"//execid(1:lengid)//
!     *     "-@."//numb(1:lengn)//".nrfai"
!      call copenfw2(fnonrfai, msg, 1, icon)
!      if(icon .gt. 1) then
!         write(0,*) ' icon=', icon
!         call cerrorMsg(msg, 1)
!         call cerrorMsg('could not be opened', 0)
!      endif
#if FNODATDEF > 0 && BUFSIZE < 200000
!      if(fnodat .gt. 0) then
         msg =input(1:lengdir)//"/"//execid(1:lengid)//
     *     "-@."//numb(1:lengn)//".dat"
         call copenfw2(fnodat, msg, 2, icon)
         if(icon .gt. 1) then
            write(0,*) ' icon=', icon
            call cerrorMsg(msg, 1)
            call cerrorMsg('could not be opened', 0)
         endif
!      endif
#endif
      if(takehist) then
         msg =input(1:lengdir)//
     *   "/"//execid(1:lengid)//"-@."//numb(1:lengn)//".hist"
         call copenfw2(fno, msg, binw, icon) ! binary write
         if(icon .gt. 1) then
            write(0,*) ' icon=', icon
            call cerrorMsg(msg, 1)
            call cerrorMsg('could not be opened', 0)
         endif
      endif
