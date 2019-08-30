!                              **************
                               common /Zmfci/
!                              **************
!

     1                 suppos,  mupol,
!
     3 hgtoa(maxhg), ictohg(maxhc), ihgtoc(maxhg),
!                        - physical consts -
     4 amfncl, amfpch, amfkch,  deltn, deltpi, deltk,
!
     5 ptavni,                  ptavfg
!
      logical suppos, mupol
!
!
! --- end of Zmfci ---
!
!    ictohg:  charge index to heavy group index
!             1, 2, 3,3, 4,4,4, 5*5, 5*6, 9*7
!    ihgtoc:  heavy group index to charge
!             1, 2, 4, 7, 12, 17, 26
!     hgtoa:  heavy group to mass #.  1, 4, 8, 14, 25, 35, 56
!
       -inc Zscl
