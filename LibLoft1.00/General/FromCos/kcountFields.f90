      subroutine kcountFields(cdata, nf)
!     count # of 'blank' or 'Tab' separated  fields in cdata
      implicit none
      character(len=*),intent(in)::cdata  ! given data. one or more blanks
                                     ! separate field.  (and/or Tab(s))
      integer,intent(out):: nf   !  # of fields 

      integer::klena    ! acutal length counting function.
      integer nc, m1, m2
      nc = klena(cdata)
      m2 = 1
      nf = 0
!     buf = " abc  xxxx  yyy  zzz "
      do while (m2  <=  nc)
         call kcgetCpos(cdata(m2:nc), m1)  ! same as kgetCpos in Epics
         m1 = m1 + m2 -1
         call kcgetBpos(cdata(m1:nc), m2)  ! same as kgetBpos in Epics
         m2 = m1 + m2 -1
         nf = nf + 1
      enddo
      end
      subroutine kcgetCpos(ichrs,  idx)
      implicit none
!        see ichrs from 1st char pos to 2nd ...
!        until non blank or Tab character appears
!        and set idx to be the such char. pos (1st char 
!        is non blank or non Tab. idx = 1)
      character(len=*),intent(in):: ichrs
      integer,intent(out)::  idx   ! idx-th pos in ichrs is first non blank
                                   ! or non Tab. 
                                   ! 0 if all are blanks

      integer k, i, klena
      character*1 tab
      tab = char(9)
      k = klena(ichrs)          ! effective length
      idx = 0 
      do i = 1, k
         if(ichrs(i:i) /= ' ' .and. ichrs(i:i) /= tab ) then
            idx = i
            exit
         endif
      enddo
      end
      subroutine kcgetBpos(ichrs, idx)
      implicit none
!            similar one as cgetCpos.   for Blank / Tab
!           
      character(len=*),intent(in):: ichrs
      integer,intent(out)::  idx ! idx-th pos in ichrs is first blank
                                   ! or Tab. 
                                   ! 0 if all are non blank or non Tab

      integer k, i, klena
      
      character*1 tab

      tab = char(9)
      k = klena(ichrs)          ! effective length
      idx = k+1
      do i = 1, k
         if(ichrs(i:i) == ' ' .or. ichrs(i:i)  == tab) then
            idx = i
            exit
         endif
      enddo
      end
      
      subroutine  kcgetaField(buf, nf, n1, n2)
      implicit none
!        find n1 and n2 such that buf(n1:n2) is the
!        nf-th field in buf.
!        If there is no nf-th field in buf, n1=0 
!        will result.
      character(len=*),intent(in):: buf
      integer,intent(in):: nf    
      integer,intent(out):: n1, n2

      integer:: klena
      integer:: nfc
      integer nc, m1, m2

      nc = klena(buf) 
      m2 = 1
      nfc = 0
      n1 = 0 
      do while (m2  .le.  nc)
         call kcgetCpos(buf(m2:nc), m1)  ! same as kgetCpos in Epics
         m1 = m1 + m2 -1
         call kcgetBpos(buf(m1:nc), m2)  ! same as kgetBpos in Epics
         m2 = m1 + m2 -1
         nfc = nfc + 1
         if(nfc == nf ) then
            n1 = m1
            n2 = m2 - 1
            exit
         endif
      enddo
      end
!!            test kcountFields etc
!      implicit none   
!      integer nf, n1, n2
!      character*45 buf
!!      buf = "abc xxxx   yyy  zzz "
!!      buf = " abc xxxx   yyy  zzz"
!      buf = "   "
!      call kcountFields(buf, nf)
!      write(*,*) ' nf =',nf
!      call kcgetaField(buf, 1, n1, n2)
!      if(n1 > 0 ) then
!         write(*,*) buf(n1:n2)
!      endif
!
!      call kcgetaField(buf, 2, n1, n2)
!      if(n1 > 0 ) then
!         write(*,*) buf(n1:n2)
!      endif
!      call kcgetaField(buf, 3, n1, n2)
!      if(n1 > 0 ) then
!         write(*,*) buf(n1:n2)
!      endif
!      call kcgetaField(buf, 4, n1, n2)
!      if(n1 > 0 ) then
!         write(*,*) buf(n1:n2)
!      endif
!
!      call kcgetaField(buf, 5, n1, n2)
!      if(n1 > 0 ) then
!         write(*,*) buf(n1:n2)
!      endif
!
!      call kcgetaField(buf, 0, n1, n2)
!      if(n1 > 0 ) then
!         write(*,*) buf(n1:n2)
!      endif
!
!      call kcgetaField(buf, -1, n1, n2)
!      if(n1 > 0 ) then
!         write(*,*) buf(n1:n2)
!      endif
!      end
