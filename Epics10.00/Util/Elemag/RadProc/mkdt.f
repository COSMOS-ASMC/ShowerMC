!     ****************************************************************
!     *                                                              *
!     *   mkdt:  makes 'data' statement usable in fortran program    *
!     *                                                              *
!     ****************************************************************
!
!   usage:
!           call mkdt( vname, d, itv, n, bfmti, lengi, jpunch )
!
!    vname: a max of 6 charcter data which is to appear in the data
!           statement like   data vname/1.0, 2.0,.../
!    d:     data to be put in data statement is given in d as d(1),
!           d(1+itv),...d(1+(n-1)*itv).
!    bfmti: a max of 8 charcter data indicating the format for one data
!           e.g., 'f9.3,   ',  'e12.3,  ', etc.  blank should follow
!           the comma to make up 8 characters.
!    lengi: width of the format above.  if this is 0, bfmti will not
!           be referenced, but standard one ('f7.4,   ') is used.
!   jpunch: if 0, data statement appears only on main sysout.
!           if non 0, the data statment will appear on device
!           specified by the device code 7
!
!
!
!
       subroutine mkdt( vname, d, itv, n, bfmti, lengi, jpunch )
!
!
!
!      real*8 vname
!       character*(*) vname, bfmti
       real*8 vname, bfmti
      dimension d(1),fmt1(6),fmt(13)
       real*8 fmt1, fmt, bfmts, cls1, cls2,  bfmt
      data fmt1/'(i6,    ','replace','replace','1h,),   ','replace','re
     1lace'/
!     'replace' means that part will be replaced on execution.
      data bfmts/'f7.4,   '/,lengs/7/,cls1/'1h,)    '/,cls2/'1h/)    '/
      data fmt/'replace','      1(','      2(','      3(','      4(','
     1    5(','      6(','      7(','      8(','      9(','     10(','
     2   11(','     12('/
!
!
      bfmt=bfmts
      if(lengi .ne. 0) bfmt=bfmti
      fmt1(5)=bfmt
      fmt(1)=bfmt
      leng=lengs
      if(lengi .ne. 0) leng=lengi
      npcd=66/(leng+1)
      npcdm=npcd-1
      npcd9=npcd*9-1
      iunit=6
      if(jpunch .ne. 0) iunit=7
      i1=1
      nc=0
      n1=1
    5 i2=min0(i1+npcdm,n)
      if(nc .eq. 9) nc=0
      nc=nc+1
      if(nc .ne. 1) go to 30
      n2=min0(n1+npcd9,n)
      do 7 iuni=iunit,6,-1
    7 write(iuni,10) vname, n1,n2
   10 format(6x,6hdata (,a6,6h(i),i=, i4,1h,i4,2h)/)
      n1=n2+1
   30 j1=(i1-1)*itv+1
      j2=(i2-1)*itv+1
      i=i2-i1+1
      fmt1(2)=fmt(i)
      fmt1(3)=bfmt
      fmt1(6)=cls2
      if(i .eq. 1) fmt1(3)=cls2
      if(nc .ne. 9 .and. i2 .ne. n) fmt1(6)=cls1
   33 do 35 iuni=iunit,6,-1
   35 write(iuni,fmt1) nc,(d(j),j=j1,j2,itv)
      if(i2 .eq. n) return
   40 i1=i2+1
      go to 5
      end
=j1,j2,itv)
      if(i2 .eq. n) return
   40 i1=i2+1
      go to 5
      end
end
