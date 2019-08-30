      subroutine callfort
      structure /xxxx/
         real*4 z
         integer ar(5)
         integer j
         logical  tf
         character*10  string(2)
      end structure
      structure /yyyy/
         record /xxxx/  a1, a2
         real*4  a3
      end structure

      record  /yyyy/ ppp
      common /abc/ ppp


      real*4  x, y
      integer i
      x= -11
      y = -22
      i = -33
      ppp.a1.z = 123
      ppp.a2.z= 345
      ppp.a2.j = -123
      ppp.a1.ar(1)=-1
      ppp.a1.ar(2)=-12
      ppp.a2.ar(3)=-123
      ppp.a2.ar(4)=-1234
      ppp.a2.ar(5)=-12345

      ppp.a2.tf = .true.
      ppp.a3 = -200
      ppp.a1.string(1)="123456789a"
      ppp.a1.string(2)="ABCDEFGHIJ"
      write(*,*) ppp.a3, ppp.a1.ar(1),  ppp.a2.z
      call xyz(ppp,  ppp.a2, 3)
      end
      integer function cleng(string)
      character*(*)  string
      write(*,*) ' in cleng ', string
      cleng = len(string)
      end
      
