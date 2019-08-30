module modepLightMaxDef
  integer, parameter::maxLightComp=500  ! max # of comp. wiht mn>0
  integer,parameter::flen=64  ! filename length including path
  integer,parameter::maxSurfaces=10 ! max # of surfaces of a component
                                   ! where light can travel.  
                                   ! as of May. 2011. octagon
  integer,parameter::maxProperties=24 ! max # of unique propterties of
                               !  component (irres. of size)
  integer,parameter::SensorNoMin=50 ! mn of sensor must >= this
  integer,parameter::SensorNoMax=99 ! mn of sensor must <= this
  integer,parameter::LightNoMin=10  ! mn of property file for light >=this
                                     ! and < SsensorNoMin 
!  real(8),parameter::EpsLength=1.d-10
  real(8),parameter::EpsLength=1.d-6
             ! (cm) same purpose as EpsLeng in epicsfile
!  real(8),parameter::EpsLength2=EpsLength*1.d5
  real(8),parameter::EpsLength2=EpsLength
             !  used to judge if EpsLength adjustment is within this range
  integer,parameter::DdiskNoDefault=13  ! direct access disk file # for stacking
  integer,parameter::MaxStackSize=3000 ! memory stack size
  integer,parameter::MaxVolAttr = 11    ! max #  vol attributes acceptable
                               !  for box  3(a,b,c), for cyl 2(r,h) 
                               !  for ecyl 3(r1,r2,h), for pipe 3(r1,r2,h)    
end module modepLightMaxDef
