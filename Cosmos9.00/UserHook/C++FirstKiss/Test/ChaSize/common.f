      subroutine fort
      structure /abc/
         real*8   x
         character*10 ch10
         integer  ii
         character*8  ch8
         character*4  ch4
c         integer  kk
      end  structure

      structure /coord/
         real*8 r(3)
         character*4 sys 
         character*1024 ca(5)
      end structure


      record /abc/  mix  
      record /coord/ mycoord
      common /ZZ/  mycoord, mix 

         mycoord.ca(2) = 'mycoord'
         mix.x = 120.
         mix.ch10 =  'mix10'
         mix.ii = -10
c         mix.kk = -20
         mix.ch8 = 'ch8'
         mix.ch4 = '4ch'
         mycoord.r(1) = 124
         mycoord.sys = 'abc'
       end

