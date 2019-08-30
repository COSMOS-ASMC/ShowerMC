#  give data before : in the lines below #-----------
#  :fromfile 
#        is to specify that basic atmomospheric data
#        is in a file or not. If t, filepath must be given in the
#        next filepath line.  The file must contian basic atmospheric data
#        which is a table crated by this program when fromfile = f. 
#  :filepath 
#        as explained above.
#  :lat & long (real)
#        If fromfile is f, lat, long,  period must be given
#        lat and long is latitude and logintude in deg. of the place
#        period specfies the period over which the average of the
#        atmospheric data is taken.
#  :period (integer)
#    period(1) is starting day (1~365). 1 is Jan.1 365 is Dec.31
#    period(2) is ending day (1~365). 
#    period(3) is starting hour (0~24). 0,24 is midnight 12 is noon
#    period(4) is ending hour. (0~24)
#    Let dayi=period(i) (i=1,2), hourj-2=period(j) (j=3,4)
#    If day1 > day2, the period is understood as straddling Dec. to Jan.
#    If hour1 > hour2, the period is understood as straddling midnight
#           
#    The data is generated by taking the average of data for sampled
#      times in the period:
#   day1                                                day2
#   *------*------*------*------*-----*------*----------*
#
#   hour1               hour2
#   *---*---*---*---*--*
#      a)  Samples days from day1 to day2 (7 day step; day1 and day2 is always
!          included. If the last but one day is close to day2 (< nearday) 
!          such day is skipped. 
!      b)  For each sampled day, get data for sampled hours between 
!        　hour1 and hour2. (4 hour step; hour1 and hour2 are always included.
!          If the last but one hour is close to hour2 (< nearhour) such  hour
!            is skipped.
!      The user can generate data, e.g, 
!      the average during day time of  the winter season (say,
#      from Dec to Feb),  or night time of the same period.
#      Also it is possible to take the   average of whole days during 
#      the summer season etc.
!      As to the height, data is generated  from lowh(-400m) to
!      heightx(500km),  above heightx extrapolation by exp formula  is used.
!      Also below lowh, extrapolation is made but it must not be large.

#---------------------------------------------
f  :fromfile  If t,  next dummy must be given 
dummy : COSMOSTOP/Data/Atmos/NRLAtmos.d  :filepath for basic atmospheric data
0   45  :lat &  long in deg. needed if fromfile is f. 
1   30   0  4  :period   Jan.  0h to 4h
