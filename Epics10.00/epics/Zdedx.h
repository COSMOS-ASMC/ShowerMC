!           common used in dedxnew dedx2

	integer nm
 
      parameter (nm=13)
 
      dimension stha(nm),  sthb(nm), sthc(nm), sthx0(nm),
     *          sthx1(nm), sthsa(nm), peak(nm), peake(nm),
     *          peak2(nm), peak2e(nm), jdef(nm), jdef2(nm)

      common /Zdedx/ stha, sthb, sthc, sthx0, sthx1, sthsa,
     *  peake, wlg0, w0, peak2, peak2e, peak, 
     *   midx, Knckon, jdef, jdef2

	real*8  stha, sthb, sthc, sthx0, sthx1, sthsa, peak
	real*8  peake, wlg0, w0,  peak2, peak2e
	integer midx, jdef, jdef2
        logical Knckon

       common /Zdedx2/ med(nm)

       character*4  med 
