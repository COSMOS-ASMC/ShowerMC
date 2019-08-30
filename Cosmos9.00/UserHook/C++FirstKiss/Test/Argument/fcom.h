         structure /abc/
	 real*8   x8
	 integer  i8
	 end structure

	 record /abc/ instabc(3)
         integer xxx
         common /com/ instabc,  xxx
       
