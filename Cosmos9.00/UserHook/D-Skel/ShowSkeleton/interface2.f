      integer leng,  lengn, lengid
      integer kgetenv2
      character*5 numb
      character*64 execid
      character*128 msg
      data kd/' g', ' e', ' m'/
      save kd

      do i= 1, nsites
         histdep(i) = 0
         indivdep(i) = 0
      enddo

      leng = kgetenv2("HISTDEP", input)
      read(input, *) histdep

      leng = kgetenv2("INDIVDEP", input)
      read(input, *) indivdep

      leng = kgetenv2("OUTPUT", input)
      read(input, *) recxy, tklat, tkelosslat, tkrespec,
     *  tkrzspec, tkzfspec,  tkrfspec, tkefspec,
     *  tkretspec, tkrtspec,
     *  tkrezspec, tkrzfspec, tkrefspec, tkarspec

      binw = 2


      input = ' '
      lengn =  kgetenv2("NUMB", numb)
      leng =  kgetenv2("OUTDIR", input)
      lengid = kgetenv2("EXECID", execid)

      if(ObserveAS) then         
         msg =input(1:leng)//"/"//execid(1:lengid)//
     *        "-@."//numb(1:lengn)//".hyb"
         call copenfw2(fnoB, msg, 1, icon)
         if(icon .gt. 1) then
            write(0,*) ' icon=', icon
            call cerrorMsg(msg, 1)
            call cerrorMsg('could not be opened', 0)
         endif
      endif

      msg =input(1:leng)//"/"//execid(1:lengid)//
     *     "-@."//numb(1:lengn)//".nrfai"
      call copenfw2(fnonrfai, msg, 1, icon)
      if(icon .gt. 1) then
         write(0,*) ' icon=', icon
         call cerrorMsg(msg, 1)
         call cerrorMsg('could not be opened', 0)
      endif
#if FNODATDEF > 0
!      if(fnodat .gt. 0) then
         msg =input(1:leng)//"/"//execid(1:lengid)//
     *     "-@."//numb(1:lengn)//".dat"
         call copenfw2(fnodat, msg, 2, icon)
         if(icon .gt. 1) then
            write(0,*) ' icon=', icon
            call cerrorMsg(msg, 1)
            call cerrorMsg('could not be opened', 0)
         endif
!      endif
#endif

      msg =input(1:leng)//
     *   "/"//execid(1:lengid)//"-@."//numb(1:lengn)//".hist"
      call copenfw2(fno, msg, binw, icon)  ! binary write
      if(icon .gt. 1) then
         write(0,*) ' icon=', icon
         call cerrorMsg(msg, 1)
         call cerrorMsg('could not be opened', 0)
      endif
