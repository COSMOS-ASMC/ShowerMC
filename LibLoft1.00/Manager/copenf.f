!  ***********************************************************
!  *     open a sequential disk file.
!  *     This is intended to open a file that exists already
!  *     for formatted mode. If it doen not
!  *     exist or cannot be opened, return cond =1 
!  ***********************************************************
         subroutine copenf(io, fnin, icon)
!           io: integer. input.  Fortran logical device number
!          fnin: character(*). input. Disk file name to be openend.
!               top blanks, if any, are removed.
!               This is needed for Gfortran(see note below)
!                     All %, #, #1, #2, @  $  treated as follows:
!                     $(XXX) type.  XXX is regarded as enevrionmental
!                       variable and replaced actual valu
!                     $XXX/ or $XXX   type.  same as above
!
!                     @ is replaced by the hostname if AtEnv is ' '
!                          If hostname contianes domainname, only hostname
!                          is extracted.   
!                     #1 is replaced by the initial Seed of the random number
!                        (1st of 2 32-bit integers)
!                     #2 is //
!                        (2nd of 2 32 bit integers)
!                     If # is not followed by 1 or 2, it
!                     is replaced by unix process number if SharpEnv is ' '.
!
!                     % is replaoced  by YYMMDDHHMMSS if PercentEnv is ' '
!                  
!                     In all cases above, if the corresponding variable
!                     (AtEnv etc is non blank (say, 'XYZ'),
!                     the envrionmental variable XYZ is assumed to exist
!                     and its value is used instead of hostname etc.
!     NOTE:
!    Suppose a progam
!     integer:: i
!     character(10):: fname
!     read(*, '(i2, a)')  i, fname
!     call copenf(11, fname, icon)
!     and this executable is a.out:     echo 10 "input" | ./a.out
!     will get fname as " input". and the open statement
!     result in error because file " input" is non existent. (Gfortran on mac).
!     (don't know  Linux Gfortran).   If echo "input" | ...  top blank is not
!     added.         
!         icon: integer. output. 0--> ok
!                                1--> cannot be opened.
         implicit none
#include "Zreadonly.h"

         character(*),intent(in):: fnin
         logical opn, ex
         integer io, ios, icon, klena, fornamelist
         character*300 msg
         character*256 fn

         fornamelist = 0
         goto 10
!        ***************
         entry  copenNLf(io, fnin, icon)
!        ***************
         fornamelist = 1
 10      continue
         call cgetfname(fnin, fn)
!                  see if already opened.
             inquire(file=fn(1:klena(fn)), opened=opn, exist=ex)
             if(opn) then
                icon = 0
             elseif(ex) then
#if defined (PCLinux) || defined (PCLinuxIFC) || defined (MACOSX) || defined (PCLinuxIFC8) || defined (PCLinuxIFC64) || defined (MacIFC)
#define SPECIAL 1
#define DELIM ,delim='apostrophe'
#else
#define SPECIAL 0
#define DELIM 
#endif

#ifdef ACTION_READ
!                     for non-writable file action ='read'
!                      is needed.
                if(fornamelist .eq. SPECIAL) then
                   open(io, file=fn(1:klena(fn)), 
     *            iostat=ios, access='sequential',
     *            form='formatted', action='read' DELIM)
                else
                   open(io, file=fn(1:klena(fn)), 
     *                  iostat=ios, access='sequential',
     *                  form='formatted', action='read')
                endif
#else
                if(fornamelist .eq. SPECIAL) then
                   open(io, file=fn(1:klena(fn)), 
     *              iostat=ios, access='sequential',
     *              form='formatted' DELIM )
                else
                   open(io, file=fn(1:klena(fn)), 
     *                  iostat=ios, access='sequential',
     *                  form='formatted')
                endif
#endif
                 if(ios .eq. 0) then
                    icon = 0
                 else
                     write(msg, *)' file=',fn(1:klena(fn)),
     *               ' exists but cannot be opened'
                     call cerrorMsg(msg, 1)
                     write(msg,*) ' see copnef.f in Manager dir'
                     call cerrorMsg(msg, 1)
                     icon =1
                 endif    
             else
                 write(msg, *) ' file=', fn(1:klena(fn)),' not exist'
                 call cerrorMsg(msg, 1)
                 icon = 1
             endif    
         end
!  ***********************************************************
!  *     open a sequential disk file.
!  *     This is intended to open a file 
!  *     for formatted i/o mode. If it
!  *     cannot be opened, return cond =1 
!  ***********************************************************
         subroutine copenfw(io, fnin, icon)
!      
!           io: integer. input.  Fortran logical device number
!           fnin: character(*). input. Disk file name to be openend.
!         icon: integer. output. 0--> ok
!                                1--> cannot be opened.
         implicit none
#include "Zreadonly.h"

         character*(*) fnin
         logical opn, ex
         integer io, ios, icon, klena, fornamelist
         character*100 msg, fn


         fornamelist = 0
         goto 10

!        *******************
         entry copenNLfw(io, fnin, icon)
!        *******************
         fornamelist = 1
 10      continue

         call cgetfname(fnin, fn)
!                  see if already opened.
             inquire(file=fn(1:klena(fn)), opened=opn, exist=ex)
             if(opn) then
                icon = 0
             elseif(ex) then
                if(fornamelist .eq. SPECIAL) then
                   open(io, file=fn(1:klena(fn)), 
     *             iostat=ios, access='sequential',
     *             form='formatted' DELIM )
                else
                   open(io, file=fn(1:klena(fn)), 
     *               iostat=ios, access='sequential',
     *               form='formatted')
                endif

                 if(ios .eq. 0) then
                    icon = 0
                 else
                     write(msg, *)' file=',fn(1:klena(fn)),
     *               ' exists but cannot be opened'
                     call cerrorMsg(msg, 1)
                     write(msg,*) ' see copnef.f in Manager dir'
                     call cerrorMsg(msg, 1)
                     icon =1
                 endif    
             else
                if(fornamelist .eq. SPECIAL) then
                   open(io, file=fn(1:klena(fn)), 
     *              iostat=ios, access='sequential',
     *              form='formatted' DELIM )
                else
                   open(io, file=fn(1:klena(fn)), 
     *               iostat=ios, access='sequential',
     *               form='formatted')
                endif
                if(ios .eq. 0) then
                   icon = 0
                else
                   icon = 3
                endif
             endif    
         end
      subroutine cskiptoEOF(iodev)
      implicit none
      integer iodev

!          skip to the end of previous write
       do while(.true.)
          read(iodev, *, end=100)
       enddo
 100   continue
       end
!  ***********************************************************
!  *     open a sequential disk file.( upgraded verson of
!  *     copenfw:
!  *     This is intended to open a file
!  *     for formatted or unformatted i/o mode.
!  ***********************************************************
         subroutine copenfw2(io, fnin,  form, icon)
         implicit none
!
         integer  io ! input.  Fortran logical device number
         character*(*)  fnin !  input. Disk file name to be openend.
         integer  form !  input. if 1--> formatted file
                       !            2--> binary file
         integer  icon !. output. 0  file is newly created and  opened
                       !          1  file exists and  opened
                       !          2  file has been already opened
                       !          3  file cannot be opened.
         logical opn, ex
         integer  ios,  klena
         character*11 format
         character*256 fn

         if(form .eq. 1) then
            format='formatted'
         elseif(form .eq. 2) then
            format='unformatted'
         else
            call cerrorMsg(
     *      'form input to chookopenfw is  invalid',0)
         endif
         call cgetfname(fnin, fn)  !  replace @ # etc to hostname etc
!                  see if already opened.
         inquire(file=fn(1:klena(fn)), opened=opn, exist=ex)
         if(opn) then
            icon = 2
         elseif(ex) then
            open(io, file=fn(1:klena(fn)),
     *           iostat=ios, access='sequential',
     *           form=format)
            if(ios .eq. 0) then
               icon = 1
            else
               call cerrorMsg(fn, 1)
               call cerrorMsg(
     *         'exists but cannot be opened', 1) 
               icon =3
            endif
         else
            open(io, file=fn(1:klena(fn)),
     *           iostat=ios, access='sequential',
     *           form=format, status='new' )
            if(ios .eq. 0) then
               icon = 0
            else
               icon = 3
            endif
         endif
         end
!           upgraded version of cskiptoEOF
      subroutine cskiptoEOF2(iodev, form)
      implicit none
      integer iodev  ! input  dev. no
      integer form   ! input  1--> ascii file
                     !        2--> binary file

!          skip to the end of previous write
       do while(.true.)
          if(form .eq. 1) then
             read(iodev, *, end=100)
          elseif(form .eq. 2) then
             read(iodev, end=100)
          endif
       enddo
 100      continue
       end

      subroutine cgetfname(fnin,  fn)
      implicit none

      character*(*) fnin  ! input. for  %,  #,  @.
                          !  see the top of file. 

      character*(*) fn    ! output. 


      integer j
      fn = ' '
      if( fnin(1:1) == "~" ) then
         fn = "$HOME/"//fnin(2:)
      else
         fn = fnin
      endif
      j = index( fn, '$')
      do while ( j > 0 ) 
         call creplst( fn, j, '$')
         j = index( fn, '$')
      enddo

      j = index (fn, '%') 
      do while ( j .gt. 0 )
         call creplst( fn, j,  '%')
         j = index (fn, '%') 
      enddo   
      
      j = index (fn, '#1') 
      do while ( j .gt. 0 )
         call creplst( fn, j, 'R')
         j = index (fn, '#1') 
      enddo   
      j = index (fn, '#2') 
      do while ( j .gt. 0 )
         call creplst( fn, j, 'r')
         j = index (fn, '#2') 
      enddo   

      j = index (fn, '#') 
      do while ( j .gt. 0 )
         call creplst( fn, j, '#')
         j = index (fn, '#') 
      enddo   

      j = index (fn, '@') 
      do while ( j .gt. 0 )
         call creplst( fn, j, '@')
         j = index (fn, '@') 
      enddo

      do 
         if(fn(1:1) /= " ") then
            exit
         endif
         fn = fn(2:)            
      enddo
      
      end
      subroutine creplst( fn, j, ch )
      implicit none
#include  "Zmanagerp.h"

      character*(*) fn   ! input.  must be < 256
                         ! output
      integer       j    ! input.  j-th chr pos. has %, # or @  $ 
      integer       jj
      character*1   ch   ! input. one of %, #, @. R, r  or $

      integer  klena, leng, kgetpid, dummy, kgetenv2, kgetnow
 
      character*64 replst   ! to contain hostname, etc to replace %, # of or @
      integer ir(2)
      character*16 envname
      character*256 fntemp
      integer jp

      jj = 1
      if( ch  == '$' ) then
         if( fn(j+1:j+1) == '(' ) then
!                 $(USER) style
            jp =index( fn(j+1:), ')') 
            if( jp > 0 ) then
               jj = jp + j -1 
               envname = fn(j+2:jj)
               leng = kgetenv2(envname, replst)
               if(leng == 0 ) then
                  write(0,*) ' envrironmetal variable =', envname
                  write(0,*) ' not found for fn=',trim(fn)
                  goto 100
               endif
               jj = jj + 1
               jj = jj-j+1
            else
               write(0,*) ' $( is used in file name  fn=',fn
               write(0,*) ' but  no counterpart ) '
               stop
            endif
         else
!            $USER/  type or $USER only
            jp =index( fn(j+1:), '/') 
            if( jp > 0 ) then
               jj = jp + j-1 
!             1       jj           
!             /tmp/$USER/ 
               envname = fn(j+1:jj)
            else
               envname = trim(fn(j+1:))
!                  USER
               jj = len(trim(fn(j+1:)))+j
            endif
            leng = kgetenv2(envname, replst)
            if(leng == 0 ) then
               write(0,*) ' Env. =', envname, ' not exist'
               goto 100
            endif
            jj = jj - j + 1
         endif
      endif

      if( ch  .eq. '@'  ) then
         if(AtEnv .ne. ' ' ) then
            leng = kgetenv2(AtEnv, replst)
            if(leng .eq. 0) then
               call cerrorMsg(
     *         'Environmental variable specified by AtEnv=', 1)
               call cerrorMsg(AtEnv, 1)
               call cerrorMsg(' Not exist ', 1)
               goto 100
            endif
         else   
            call cgetHost(leng, replst)
         endif
      endif
      if( ch .eq. 'R') then
         call cqIniRn(ir)
         replst = ' '
         write(replst,'(i11)') ir(1)
         call kseblk(replst, '{', leng)
         jj=2
      endif

      if( ch .eq. 'r') then
         call cqIniRn(ir)
         replst = ' '
         write(replst,'(i11)') ir(2)
         call kseblk(replst, '{', leng)
         jj=2
      endif

      if(  ch  .eq.  '#' )  then
         if(SharpEnv .ne. ' ' ) then
            leng = kgetenv2(SharpEnv, replst)
            if(leng .eq. 0) then 
               call cerrorMsg(
     *         'Environmental variable specified by SharpEnv=', 1)
               call cerrorMsg(SharpEnv, 1)
               call cerrorMsg(' Not exist ', 1)
               goto 100
            endif
         else
            replst = ' '
            write(replst, '(i10)') kgetpid(dummy)
            call kseblk(replst, '{', leng)
         endif
      endif

      if( ch  .eq. '%' ) then
         if(PercentEnv .ne. ' ' ) then
            leng = kgetenv2(PercentEnv, replst)
            if(leng .eq. 0) then 
               call cerrorMsg(
     *         'Environmental variable specified by PercentEnv=', 1)
               call cerrorMsg(PercentEnv, 1)
               call cerrorMsg(' Not exist ', 1)
               goto 100
            endif
         else
            replst = ' '
            leng =  kgetnow(replst)  ! get YYMMDDHHMMSS
            leng =  klena(replst)
         endif
      endif
      if( j .eq. 1 ) then
         fntemp = replst(1:leng)//fn(j+jj:klena(fn))
         fn = fntemp
      else
         fntemp =
     *        fn(1:j-1)//replst(1:leng)//fn(jj+j:)
         fn = fntemp
      endif
      return
 100  continue
      write(0,*) ' input fn=',fn
      stop
      end
!      **************
      subroutine cgetHost(leng, hostn)
      integer  leng       ! output
      character*(*) hostn ! output

      character*1 NULL
      integer  kgetenv, j

      NULL = char(0)
      leng = kgetenv("HOSTNAME"//NULL, hostn)
      if(leng .eq. 0) then
         leng = kgetenv("HOST"//NULL, hostn)
         if(leng .eq. 0) then
            call
     *      cerrorMsg('Env. var. HOST or HOSTNAME not found',0)
         endif
      endif

      j =index(hostn, '.') 

      if(j .gt. 0) then
         leng = j-1
         hostn = hostn(1:leng)
      endif
      end
