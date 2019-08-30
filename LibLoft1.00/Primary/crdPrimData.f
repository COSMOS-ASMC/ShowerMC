!     include '../Sysdep/copenf.f'
!     include '../Sysdep/cskipComment.f'
!c            test crdPrimData etc
!c
!     include '../Particle/Zptcl.h'
!     include 'Zprimary.h'
!     type(primaries):: Prm
!     call copenf(TempDev, 
!    * '../Data/Primary/sample.d')
!     call cskipComment(TempDev, icon)
!     call crdPrimData(Prm)
!     write(*, *) ' no. of comps=', Prm.no_of_comps
!     end
!       *************************************************
!       *   crdPrimData:  read primary information
!
        subroutine crdPrimData(prm)
        implicit none
#include  "Zmanagerp.h"
#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryc.h"
#include  "Zprimaryv.h"

        type(primaries):: prm
        type(component):: onecomp
!
        integer icon, n
!
        n = 0
        icon = 0
        do while(icon .eq. 0)
!            call crdCompPrim(prm.each(n+1), 
            call crdCompPrim(onecomp,
     *      n+1, icon)        ! each component of primaries.
            if(icon .eq. 0) then
               n = n + 1
               if(n .gt. maxNoOfComps) then
                 write(0,*) ' # of 1ry comp. > maxNoOfComps=',
     *            maxNoOfComps
                 write(0,*) ' change the value in Zmaxdef%h '
                 stop
              endif
              prm%each(n)= onecomp
            endif
        enddo

        close(TempDev)

        if(n .le. 0) then
           write(*,*) ' no primary is given'
           stop 9999
        endif
        prm%no_of_comps = n
      end
!     ***********************************************
      subroutine crdCompPrim(each, n, icon)
!     ***********************************************
!          each: /component/. Output.  one composition of primary data is read
!            n: integer. Input. n+1-th block in the table is next data
!         icon: integer. output. 0--> more data
!                                1--> eof reached.
        implicit none

#include  "Zmanagerp.h"
#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryc.h"
        integer n   !n+1-th primary data in the table is being read
        integer  icon
        type(component):: each
        integer ios, np    ! np is the data points counter
!                 read 1ry type, energy unit, energy type, etc
         each%cut = 0.
         each%cut2 = 0.
         read(TempDev, *, iostat=ios) each%symb, each%eunit,
     *   each%etype, each%diff_or_inte,  each%flatterer, each%cut,
     *   each%cut2
!
         if(ios .eq. 0) then
!                 ! not yet eof or no error
              np = 0       ! segment counter
              do while(.true.)  ! read until end of 1 block
                  read(TempDev, *, iostat =ios)
     *            each%energy(np+1),
     *            each%flux(np+1)
                  if(ios .ne. 0) then
                      write(*, *) ' error in the primary data'
                      write(*, *) ' data position: Primary block=', 
     *                n+1,' segment # =', np + 1
                      close(TempDev)
                      stop 9999
                  endif
!                     see if end of the composition
                  if(each%energy(np+1) .le. 0.) goto 100
                  np = np + 1
                  if(np .gt. maxSegments+1) then
                     write(*, *) 'too many segments in ',n+1,
     *              '-th composition. must be <',maxSegments+1
                     close(TempDev)
                     stop 9999
                  endif   
             enddo
 100         continue
             each%no_of_seg =max(0, np -1) ! segments
!                 note that each.no_of_seg is data points -1
!                (2 points makes a segment !)
             icon = 0
        else
             icon = 1    ! eof
        endif
      end
 









