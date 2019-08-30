      subroutine epReadXXsec(mediain, icon)
      implicit none
#include "ZepManager.h"
#include "ZmediaLoft.h"
       type(epmedia)::  mediain

      character*8 name   !   media name.
      integer n            !-->size  output.  number of rows of xcomtab
      integer icon          ! output.  icon = 0.  xcombtab obtained
                            !  -1.  no. tab found.

      integer klena
      integer i, j, result, k, ns, is
      character*100 mdpath
      
      icon = -1
      name = mediain%name

      do i = 1, MaxMediaDir
         if( klena(MediaDir(i)) .gt. 0 ) then
            mdpath = ' '
            mdpath = MediaDir(i)(1:klena(MediaDir(i)))
     *            // '/'  //name(1:klena(name))//".xcom"
            call copenf(iowk, mdpath, result)
            if( result .ne. 0) then
               write(0,*) MediaDir(i)(1:klena(MediaDir(i)))
               write(0,*) name(1:klena(name))
               write(0,*) mdpath
               write(0,*) ' cannot be opended'
               stop 999
            endif
            call cskipComment(iowk, icon)
            if( icon .ne. 0) then
               write(0,*)
     *              "#------------- line not found in .xcom"
               stop 987
            endif
            n = 0
            do while(.true.)
               read(iowk, *, end=100)
     *            ( mediain%xcom%tab(k,n+1), k = 1, 8)
               n = n + 1
            enddo
 100        continue
            close(iowk)

            mediain%xcom%size = n
!              count shell energy levels
            ns = 0
            do is = 1, n-1
               if( mediain%xcom%tab(1,is) .eq.
     *             mediain%xcom%tab(1,is+1) ) then 
                  ns = ns +1
                  if(ns .gt. maxNoOfShells) then
                     write(0,*) ' media=',name, ' has ',
     *                 ' shells > ', maxNoOfShells
                     stop 1235
                  endif
                  mediain%pe%shellE(ns)=mediain%xcom%tab(1,is)
               endif
            enddo

            mediain%pe%noOfShells=ns

!             we use log table.
            icon = 0
            do k = 1, n
               do j= 1, 8
                  if( mediain%xcom%tab(j, k) .eq. 0.) then
                     mediain%xcom%tab(j,k) = -100
                  else
                     mediain%xcom%tab(j,k) = log( mediain%xcom%tab(j,k))
                  endif
               enddo
            enddo
!           *****
            write(0,*) mdpath(1:klena(mdpath)),
     *       ' has been read, no of shells=',ns
            exit
!           ****
         endif
      enddo

      end
