!      read matter table

      subroutine epReadMTbl(io,media)
      implicit none
#include "Zmedia.h"
      integer io  ! input.  logical  dev. number of input data
       type(epmedia):: media  !  output. obtained mediea characteristics
!    if media.format = 1,  the following is obtained.
!      media.noOfElem,
!      media.rho
!                                                            H2O
!      media.elem(i).A, i = 1,media.noOfElem  mass number   1, 16  
!      media.elem(i).Z, //         atomic number            1, 8
!      media.elem(i).N, //         nucleon number           1  16
!      media.elem(i).No, //   number of elem.       2  1
!      media.gasF        gasF =0 for solid, 1 for gas. 
!      media.n           reflactive index. if unknown, 0
!      media.BirksC1, C2, CC Birks correction coeff. for organic sciniti.
!      media.name
!   if media.format = 2  following additional data is assumed to follow.
!#      (if some are  not known,  -100 is given)
!#    œò§èZ/Aœò§é   I[eV]   a       k      x0    x1     Cbar  delta0
!    .42065  534.1 0.09569 3.0781 0.0456 3.7816  5.7409 0.00   
!     example above   is for BGO.
      character*80  error
      character*120 line

      character*24 field(20)
      integer nf, icon
      integer klena,  i
      data error/'not enough data for basic media'/
      

      call epget1line(error, io, line, icon)
!           decompose into field
      call kgetField(line, field, 20, nf)         
      if(nf .ne. 2) then
         call cerrorMsg('First data in media def. file invalid',1)
         call cerrorMsg(line, 0)

      endif
      media%name = field(1)(1:min(klena(field(1)),8))

      read(field(2), *)  media%format
      if(media%format .ne. 1 .and. media%format .ne. 2 ) then
         write(0,*) 'media foramt =', media%format, ' not acceptable'
         stop
      endif
      call epget1line(error, io, line, icon)
!         this should be  
!         Elem  rho  Gas/Solid line
      call kgetField(line, field, 20, nf)
      if(nf .lt. 4) then
         write(*,*) ' nf=', nf,' :  ', field(1), ':',
     *           field(2),' : ', field(3), ' : ', field(4)
         call cerrorMsg(line, 1)
         call cerrorMsg('Elem, rho, gas/solid  line  is short', 0)
      endif
      if( nf == 8 ) then
         if( field(8) == "B" ) then
!              Birks formula
            read(line, *) media%noOfElem, media%rho, media%gasF,
     *               media%n, 
     *              media%BirksC1, media%BirksC2, media%BirksCC
            media%Birks="B"

         elseif( field(8) == "L" ) then
!                  log formula:  
!           dE/dx)q = X* (X*cc)  **( - c2*log(c1*X) )
!           where  X = dE/dx + 1./cc and dE/dx is unquenched one in GeV/(g/cm2)
!           c1 and cc in g/cm2/GeV and c2 is simple number
!           If fitting is done, say, for dE/dx in MeV/(g/cm2)
!           and C1, C2, CC are obtained,  c1=C1*1000, c2=C2, cc=CC*1000
!           must be used 
!
!              typical possilbe values: c1=23.5 c2= 0.086 cc = 4.61
!
            read(line, *) media%noOfElem, media%rho, media%gasF,
     *               media%n, 
     *              media%BirksC1, media%BirksC2, media%BirksCC
            media%Birks="L"
         else
            write(0,*) " media name =", media%name
            write(0,*) " input line=",line
            write(0,*) " is strange: # of fields =8, then the last one"
            write(0,*) " should be 'B' or 'L'  which stands for Birks"
            write(0,*) " or log formula "
            stop
         endif
      elseif( nf == 7 ) then
         if( field(7) == "T" ) then
!                  Tarle  dE/dx)q =  (1-c2)X/(1 + c1(1-c2)X) + c2X
!            where X is unquenched dE/dx. in GeV/(g/cm2).
!            c1 is in g/cm2/GeV C2 and simple number
            read(line, *) media%noOfElem, media%rho, media%gasF, 
     *              media%n, 
     *              media%BirksC1, media%BirksC2
            media%Birks = "T"

         else
!              old Birks format ( without "B") is assumed
            read(line, *) media%noOfElem, media%rho, media%gasF,
     *              media%n, 
     *              media%BirksC1, media%BirksC2, media%BirksCC
            if( media%BirksC1 > 0.) then
               media%Birks ="B"
            else
               media%Birks =" "
            endif

         endif
      elseif( nf == 4 ) then
         read(line, *) media%noOfElem, media%rho, media%gasF, 
     *              media%n 
         media%Birks = " "
      else
         write(0,*) " media name =", media%name
         write(0,*) " input line=",line
         write(0,*) " strange "
         stop
      endif
!         for gas n is   (n-1)*10^6. 
      if(media%gasF .eq. 1) then
         media%n = media%n/1.d6 + 1.d0
      endif

      do  i = 1,  media%noOfElem
         call  epget1line(error, io, line, icon)

         read(line, *)
     *           media%elem(i)%Z, media%elem(i)%A, media%No(i)
!            Nucleon number
         media%elem(i)%N = media%elem(i)%A + 0.4999999 +
     *           media%elem(i)%Z
      enddo
      if(media%format .eq. 2) then
         call epget1line(error, io, line, icon)
!#    œò§èZ/Aœò§é   I[eV]   a       k      x0    x1     Cbar  delta0
         read(line, *) media%ZbyAeff, media%I, media%sh%sa, media%sh%k,
     *           media%sh%x0, media%sh%x1, media%sh%c, media%sh%delta0
!          we convert media.I to GeV here to avoid later confusion
         if( media%I /= -100.0 ) then
            media%I = media%I*1.d-9
         endif
!       NOTE:    our sh.c(<0) is -Cbar, then
         if( media%sh%c /= -100.0 ) then
            media%sh%c = -media%sh%c
         endif
      endif
!        if this is used not via config file, rhoc must be give
!       (e.g when making sampling table).
      media%rhoc = 1.

      end
!     **********************************************
      subroutine epget1line(error,  io, line, icon)   
      implicit none
!        get 1 line (skip # lines or blank line)
      character*(*)  error ! input. if E.O.F is encountered
                           !        this message is output.  and stop is made
                           !     if error=' ', no output, no stop
      integer io          ! input. dev number
      character*(*) line  ! output. 1  line read from io

      integer icon        ! output. 0 if line is obtained

                          !         1 if E.O.F 
      icon = 0
      line = ' '
      do while(.true.)
         read(io, '(a)', END=110)  line
         if(line(1:1) .ne.  '#' .and. line .ne. ' ') goto 10
      enddo

 110  continue
      if(error .ne. ' ') then
         call  cerrorMsg(error, 0)
      endif
      icon = 1
 10   continue
      end

      
         
