      subroutine epBPgeneini(file, media, how )
      use BPLPM
      implicit none
#include "Zmedia.h"

      character(*),intent(in)::file ! path to the basic media file
              ! may be Air*0.01 type
      type(epmedia),intent(out)::  media  !  media data is read

!          initialization before using epBrgenex, epBrgene
!                                      epPrgenex, epPrgene
!         Besides this,
!         you have to give Eeme for epBrgenex, epBrgene
!                          Egme for epPrgneex, epPrgene
!         in ZBPgene.h
!
      integer,intent(in):: how  !  =-1 is normal usage.  see epNormBPxs below.
!
      integer i, icon
      integer,parameter::io=5
      real(4):: rhoc
      character(:),allocatable:: tfile

      tfile=file  
      call epgetRhoc(file, tfile, rhoc)  !../Air*0.01==>../Air and rhoc  
      call copenf(io, tfile, icon)
      if(icon /= 0 ) then
         write(0,*) 'tfile=',trim(tfile), '---- could not  be opend---'
         stop
      endif

      write(0,*) ' entering epReadMTbl'
      call epReadMTbl(io, media)
      close(io)
!
!     media%rho = media%rho*rhoc  ! this is easy way but
! in actual M.C where the same media may be
!           used with diff. rhoc, N.G      
      media%rhoc =rhoc


      write(0,*) ' entering epGetEffZA'
      call epGetEffZA(media)
      write(0,*) ' entering epSetSTblCns'
      call epSetSTblCns(media, media%cnst)
!          init for Seltzer function
      write(0,*) ' entering epZforSelt'
      do i = 1, media%noOfElem
         call epZforSeltz(i, media%elem(i)%Z)
      enddo

      write(0,*) ' bpsetmedi '
      call epBPSetMedia(media)  ! copy media into BPgene.h

      call epBPSetForce('?')
      call epBPnormXs(  media, how )
      write(0,*) 'exiting  geneini'
      end       subroutine epBPgeneini
