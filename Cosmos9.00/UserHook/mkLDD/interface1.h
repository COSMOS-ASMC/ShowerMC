      do i= 1, nsites
         histdep(i) = 0
      enddo
                                             
      call cqUhookc(1, input)
      read(input, *) histdep


      call cqUhookc(2, input)
      read(input, *)  tkrtspec, tkarspec, tkweb

      call cqUhookc(3, input)
      basefilename = input

      call cqUhooki(1, binw)  ! by binw,  specify binary or ascii write of histogaram
                          
      call cqUhooki(2, reducedTime)
      if(reducedTime .ne. 1 .and. reducedTime .ne. 2) then
         write(0,*) ' UserHooi(2)=',UserHooki(2), ' invalid'
         stop
      endif

      
