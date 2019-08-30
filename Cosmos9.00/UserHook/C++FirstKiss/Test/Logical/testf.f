      subroutine ddd
      structure /abc/
         real*8  dv, dv2
         logical torf, torf2
         integer  iii
      end structure
      
      record /abc/ test1(2)
      common /yyy/ test1

      test1(1).dv = 2.7182d0
      test1(1).dv2 = 3.14d0
      test1(1).torf = .true.
      test1(1).torf2 = .false.
      test1(1).iii = -100

      test1(2).dv = 2.7182d0
      test1(2).dv2 =-3.14d0
      test1(2).torf = .true.
      test1(2).torf2 = .false.
      test1(2).iii = -200
      end

