      subroutine creadSoftenPara(io)
      implicit none
      integer,intent(in):: io   ! logical dev. #
      character*24 vname
      character*100 vvalue


       call cskipsep(io)
       do while( cgetParmN(io, vname, vvalue ) )
          select case(vname)
          case('mode')
             call creadParaI(vvalue, mode)
          case('xth') 
             call creadParaR(vvalue, xth)
          case('E0th') 
             call creadParaR(vvalue, E0th)
          case('fwbw')
             call creadParaI(vvalue, fwbw)
          case('pw')
             call creadParaR(vvalue, pw)
          case('repeat')
             call creadParaR(vvalue, repeat)
          case('special')
             call creadParaL(vvalue, special)
          case('leadingPiK')
             call creadParaL(vvalue, leadingPiK)
          case('useXinCMS')
             call creadParaL(vvalue, useXinCMS)
          end select
       enddo
       end       subroutine creadSoftenPara
!      *************
       subroutine  cwriteSoftenPara(io)
       implicit none
       integer,intent(in):: io

       write(io,*)'----------------------'
       call cwriteParaI(io,'mode', mode)
       call cwriteParaR(io,'xth', xth)
       call cwriteParaR(io,'E0th', E0th)
       call cwriteParaI(io,'fwbw', fwbw)
       call cwriteParaR(io,'pw', pw)
       call cwriteParaR(io,'repeat',repeat)
       call cwriteParaL(io,'special',special)
       call cwriteParaL(io,'leadingPiK',leadingPiK)
       call cwriteParaL(io,'useXinCMS', useXinCMS)
       end       subroutine  cwriteSoftenPara


       subroutine cskipsep(io)
       implicit none
       integer io
       character(10)  sep
       do while (.true.)
          read(io, '(a)') sep
          if(sep(2:10) == '---------') exit
       enddo
       end  subroutine cskipsep
!        ************************* real*8 data
       subroutine creadParaR(vvalue, x)
        implicit none
        integer io
        character*(*) vvalue
        real*8 x
!        read(vvalue, *)   x, x
        read(vvalue, *)   x
        end       subroutine creadParaR
       subroutine creadParaR2(vvalue, x, n)
        implicit none
        integer io
        character*(*) vvalue
        integer n
        real*8 x(n)
        read(vvalue, *)   x
        end       subroutine creadParaR2

!     ************************* complex data
      subroutine creadParaCx(vvalue, c)
      implicit none
      character*(*) vvalue
      complex*8 c
      read( vvalue, *)   c
      end      subroutine creadParaCx
!     ************************ integer data
      subroutine creadParaI(vvalue, i)
      implicit none
      character*(*) vvalue
      integer i
      read(vvalue, *)   i
      end      subroutine creadParaI
!        ************************* character data
      subroutine creadParaC(vvalue, cha)
      implicit none
      character*(*) vvalue
      character*(*) cha
      read(vvalue, *)  cha
      end      subroutine creadParaC
!        ***************************** logical data
      subroutine creadParaL(vvalue, logi)
      implicit none
      character*(*) vvalue
      logical logi
      read(vvalue, *)  logi
      end           subroutine creadParaL
!        ---------------------------------------------
      subroutine cwriteParaR(io, vname, x)
      implicit none
      integer io
      character*(*) vname
      real*8  x
      
      write(io, *) ' ', vname,' ', x,' /'
      end      subroutine cwriteParaR
      subroutine cwriteParaR2(io, vname, x, n)
      implicit none
      integer io
      integer n  ! arra size of x
      character*(*) vname
      real*8  x(n)
      
      write(io,*) ' ', vname,' ', x,' /'
      end      subroutine cwriteParaR2

      subroutine cwriteParaCx(io, vname, c)
      implicit none
      integer io
      character*(*) vname
      complex*8  c
      write(io,  *) ' ', vname,' ', c,' /'
      end subroutine cwriteParaCx

      subroutine cwriteParaI(io, vname, i)
      implicit none
      integer io
      character*(*) vname
      integer i
      
      write(io,  *) ' ', vname,' ', i,' /'
      end      subroutine cwriteParaI

      subroutine cwriteParaC(io, vname, cha)
      implicit none
      integer io
      character*(*) vname
      character*(*) cha
      integer klena
      character*2 qmk/" '"/             ! ' 
      if(klena(cha) .gt. 0) then
         write(io,  *) ' ', vname, qmk, cha(1:klena(cha)),
     *        qmk,' /'
      else
         write(io, *) ' ', vname, qmk, ' ', qmk, ' /'
      endif
      end      subroutine cwriteParaC
      subroutine cwriteParaL(io, vname, logi)
      implicit none
      integer io
      character*(*) vname
      logical  logi

      write(io,  *) ' ', vname,' ', logi,' /'
      end      subroutine cwriteParaL

       function cgetParmN( io,  vname, vv ) result(ans)
!          get parameter variable name and given value(s)  from io
       implicit none
       integer io
       character*(*)  vname, vv  ! output
       logical ans

       integer linel
       parameter( linel = 100)
       character*(linel)  line
       integer loc, loc2
       vname = " "
       do while(.true.)
          read(io, '(a)', end=100 ) line
          if( line(1:1) .eq. " " .and. line(2:2) .ne. " ") then
             loc = index( line(3:linel), " ")  + 2
             vname = line(2:loc-1)
             loc2 = index( line, "/")
             if(loc2 .eq. 0 ) then
                write(0,* ) ' "/"  is missing in the input data file '
                write(0,*)  ' The line is: ', line
                stop 1234
             endif
             vv = line(loc+1:linel)  !  some data containes '/' such as '../../Media' so put all
                                     ! data.
             goto 50
          endif
       enddo
 50    continue
       ans = .true.
       return
 100   continue
       ans =.false.
       end function cgetParmN

