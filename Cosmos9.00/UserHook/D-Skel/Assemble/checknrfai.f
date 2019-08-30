#include "../FleshHist/asinfo.f"
#include "../FleshHist/asdensity.f"
#include "Zmaxdef.h"
#include "../FleshHist/crecprob.f"

!       read final .nrfai data and
!       see if the probabilty of recording partilces;
!       is it not too small or not too large ?
      implicit none
#include "Zobs.h"
#include "../FleshHist/Zprivate0.h"
      integer ndepth
      parameter (ndepth= nsites)
      real nrfaiRec0(nrbin, nfai, 4, ndepth)
      real nrfaiAll0(nrbin, nfai, 4, ndepth)
      real recprob(nrbin, 4, ndepth)

      integer EvNo0

      real limit(4), E0, Emin
      real rat, all, rec, prob
      integer klena
      real intdep(ndepth)
      integer nrbina, nfaia, ansites0
      integer leng, i, j, k, l
      integer i0, j0, k0, l0
      integer NN
      integer newfmt, maxsites, nr
      integer icon0, iconx, icont
      character*128 input0
      real  cosz, age, sum, Nx, depth
      character*20 field(15)
      real nptcls(nrbin, 4,  ndepth)
!      limit= 5000.
      E0=1.e8
      cosz= 1.0
      do while(.true.)
         input0 = ' '
         read( *, '(a)',  end=1000 ) input0
         if(input0 .ne. " ") then
            call  ksplit(input0,  20, 15, field,  nr)
            if(nr .eq. 8) then
               read(input0(1:klena(input0)), *)
     *           EvNo0, E0,NN, cosz, limit(1), nrbina, nfaia, ansites0
               maxsites = ansites0
               Emin=500.d-6
               newfmt = 0
            elseif(nr .eq. 9 ) then
               read(input0(1:klena(input0)), *)
     *           EvNo0, E0,NN, cosz, limit(1), nrbina, nfaia, ansites0,
     *           maxsites
               Emin=500.d-6
               newfmt = 1
            elseif( nr .eq. 13) then
               read(input0(1:klena(input0)), *)
     *           EvNo0, E0,NN, cosz, Emin, nrbina, nfaia, ansites0,
     *           maxsites, limit
               newfmt = 2
            endif

            if(nrbina .ne. nrbin .or. nfaia .ne. nfai) then
               write(0,*)' nrbina=',nrbina, 'or  nfaia=',nfaia,
     *              ' differ from the def. in this prog'
               stop 5555
            endif
!           ********
            do i = 1, ansites0
               do j = 1, 4
                  do k= 1, nfai
                     read(*, '(3x, f7.1,  4i4)' )
     *                    intdep(i), l0, i0, j0, k0
                     if(i0 .ne. i .or. j .ne. j0 .or. k .ne. k0) then
                        write(0,*) ' intdep, i0,j0,k0=',
     *                      intdep(i), i0, j0, k0, ' strange'
                        stop 8888
                     endif
!
                     read(*, *)
     *                    ( nrfaiRec0(l,k,j,i), l=1,nrbin )
                  enddo
               enddo
            enddo
!    ************
            do i = 1, ansites0
               do j = 1, 4
                  do k = 1, nfai
                     read(*, '(3x,f7.1, 4i4)' )
     *                 intdep(i), l0, i0, j0, k0
                     if(i0 .ne. i .or. j .ne. j0 .or. k .ne. k0) then
                        write(0,*) ' intdep, i0,j0,k0=',
     *                      intdep(i), i0, j0, k0, ' strange'
                        stop 9876
                     endif
                     read(*, *)
     *                    ( nrfaiAll0(l,k,j,i), l=1,nrbin )
                  enddo
               enddo
            enddo
         else
!            1 event read
            if(newfmt .eq. 0) then
               write(0, '(i2, 1pE11.3, 0pf7.1, 1pE11.3, 3i4)')
     *       EvNo0, E0,NN, cosz, limit(1),  nrbina, nfaia, ansites0
            elseif( newfmt .eq. 1) then
               write(0, '(i2, 1pE11.3, 0pf7.1, 1pE11.3, 4i4)')
     *       EvNo0, E0,NN, cosz, limit(1),  nrbina, nfaia, ansites0,
     *           maxsites
            elseif( newfmt .eq. 2) then
               write(0,
     *            '(i2, 1pE11.3, 0pf7.1, 1pE11.3, 4i4, 4f10.0)')
     *           EvNo0, E0,NN, cosz, Emin,  nrbina, nfaia, ansites0,
     *           maxsites, limit 
            endif
            write(*,'(12a)') 'prob ', ' accpt ', ' accuracy ', 
     *      ' rec',  ' all', ' exp-all', 
     *       ' dep ', ' code', ' fai',' r ','  age', ' dep' 
            do i = 1, ansites0
               depth =intdep(i)
               do j = 1, 4
                  do k = 1, nfaia
                     if(k .eq. 1) then
                        call crecprob(depth, j, limit, 
     *                  360.0/nfai, E0,
     *                  cosz, nrbin, recprob(1, j,  i),
     *                  nptcls(1, j, i),  age,  sum, Nx)
                     endif

                     do l = 1, nrbin
                        if( nrfaiAll0(l, k, j, i) .gt. 0.) then
                           prob=recprob(l, j, i)
                           rec = nrfaiRec0(l,k,j,i)
                           all = nrfaiAll0(l,k,j,i)
                           rat=rec/limit(j)
                           write(*,'(1p6E11.3,4i4, 0p,f7.3,f7.1)')
     *                      prob, rec/all, rat, rec,
     *                      all, nptcls(l, j, i), i,j,k,l, age, depth
                        endif
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
 1000 continue
      write(0,*) 'all events processed '
      end
