      subroutine creadAtmosD
      use modAtmosDef
      implicit none
#include  "Zmanagerp.h"
#include  "ZmediaLoft.h"
      
      integer  ios, nodes, icon, nf

      character*150 msg
      character(100):: line
      character*60 adata(maxnodes)

      integer i
!          adata is the same one as
      nodes = 0
      atmos%matter(:) = "  "
      NoOfMedia = 0      
      do while(.true.)
         if(AtmosFile .ne. '     ') then
            read(TempDev, '(a)',end=10, iostat=ios) line
         else
            line = adata(nodes+1)
         endif
         call kcountFields(line, nf)
         if( nf > 0 ) then
            nodes = nodes + 1
            if( nf == 4 ) then
               read(line, *)
     *           atmos%z(nodes), atmos%T(nodes),
     *           atmos%P(nodes), atmos%rho(nodes)
            elseif( nf == 5 ) then
               read(line, *)
     *           atmos%z(nodes), atmos%T(nodes),
     *           atmos%P(nodes), atmos%rho(nodes), atmos%matter(nodes)
            else
               write(0,*) ' error in adata or input data creadAtmosD.fh'
               write(0,*) ' data: ', line
               write(0,*) ' nf=', nf
               stop
            endif
         else
            goto 10
         endif
!         if(ios .ne. 0) then
!            call cerrorMsg(
!     *       'something wrong about atmosphere data file:',1)
!            call cerrorMsg(AtmosFile, 0)
!        endif
      enddo
 10   continue

      atmos%nodes = nodes
         !         assume default is Air
      if( atmos%matter(1) == "  " )  atmos%matter(1) = "Air"
      do i = 2, nodes
         if( atmos%matter(i) == " " )
     *            atmos%matter(i) = atmos%matter(i-1)
      enddo
!     Count # of uniq media (normally 1; only Air)
!     and read media data (by eprdMFfile) to be ready for EM
!      process sampling
      call cgetMedia

      
      if(AtmosFile .ne. '  ') then
         close(TempDev)
         write(msg, *) 
     *   "Atmosphere data has been read: # of nodal points =",
     *    nodes
      else
         write(msg, *) 
     *   "Atmosphere data stored internally is read nodes=",
     *    nodes
      endif
      call cerrorMsg(msg, 1)
      end   subroutine creadAtmosD


