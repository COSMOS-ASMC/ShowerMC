#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Z90histfuncdef.h"

/*
c     
c            weighted histograming  fortan 90 version
c            (Not work under Absoft fortran 90) 
c      Usage:  kwhisti:   insatnciate one histogram
c              kwhistc:   clear histogram area
c              kwhist:    take histogram
c              kwhists0:  sepcify the integaral from -inf or +inf
c              kwhists:   compute statistical result.
c                        This can be used more than once 
c                        with a differennt normalization factor
c                        for the same histogram.   
c                 Therefore, in bin2ascii you may use inorm=-1 to 
c                    keep the previous value.
c              kwhistIxy: get x,y for integral distribution
c              kwhistxy:  get x,y for x and dn/dx
c              kwhistp:   print statistical result or binary write histogram result
c                         by calling kwhistpr or kwhistw
c              kwhistw:  write histogram with binary format
c                        for later use.
c              kwhistr:  read histogram written by kwhistw
c              kwhista:  add two histograms with identical
c                        structure.
c              kwhistpr: print histogram 
c              print format
c         #hist1 

c         ...
c         ...
c         0 0 0 0 0 0
c
c      MinIndex: min. bin index where non zero data is stored
c      MaxIndex: max. bin index where non zero data is stored
c
*/

//      subroutine kwhisti(h, ixmin, ibinORxmax, inbin, itklg )
#include "Z90histc.h"
#include "Z90histo.h"
#include "Z90hist1.h"
#define TRUE 1
#define FALSE 0


static float normf;
static int BinWrite;
static int fromwhich=0;

static  char dirstr[96];



void kwhisti(struct histogram1 *h, float ixmin, float ibinORxmax,int inbin, int itklg) {
  int ndiv;
  int asmax;
  //      implicit none

  //c         instanciate
  //    integer inbin  ! input. request inbin histogram area
  // int inbin;
  //      real ixmin     ! input. xmin. not in log even if log10(variable) is taken
  //  float ixmin;
  //                     !         see itklg 
  //      real ibinORxmax  ! input. bin or ixmax. depends on itklg.
  //  float ibinORxmax;
  /*
                     !  If bin and log10 is taken, bin is for log10 

      integer itklg  ! input.  bit pattern. give it like b'10001'
  */
  //  int itklg;
  /*
                     !         bit 1 is LSB.
                     !         bit 1: 0--> not take log10 of variable
                     !                1--> take log10   //
                     !             2: 0--> ixmin is the min of lowest bin
                     !                    |---|---|---|....     |...|
                     !                    |                         |
                     !                    ixmin                     ixmax
                     !                1--> ixmin is the center of the lowest bin
                     !                  |--*--|-----|-----|....    |--*--|
                     !                     |                          |
                     !                     ixmin                      ixmax
                     !            max follows the same rule.
   
                     !             3: 0--> neglect underflow
                     !                1--> underflow is put in lowest bin
                     !                     mean bin value is affected by
                     !                     those with underflowed values
                     !             4: 0--> neglect overflow
                     !                1--> overflow is put in the highest bin
                     !                     mean bin value is affected by
                     !                     overflowed ones    
                     !             5: 0-->ibinORxmax  is the bin
                     !                        xmax is determined by bin,
                     !                        xmin and inbin
                     !                1-->ibinORxmax  is ixmax. 
                     !                        bin is determined by xmax xmin
                     !                        and inbin.

c     ******************
  */
  //      type(histogram1) h, h1, h2  ! no need


//c     ====================      
//      integer fno  !  if < 0, standard output is used else fno is used for histogram output

   //      character*96  dirstr

   //      integer klena  ;   not used in C
   // save normf
   /*
      if( h%c%init .eq. 'initend') then
         write(0, *) '1D hist already instanciated'
         write(0, *) ' title=',h%c%title
         write(0, *) ' categorye=',h%c%categ
         write(0, *) ' id=',h%c%id
         stop 9999
      else
         h%c%init = 'initend'
      endif
   */
  if( strcmp(h->c.init, "initend" ) == 0 ){
    fprintf(stderr, "1D hist alrady instanciated\n");
    exit(1);
  }
  else {
    strcpy(h->c.init, "initend");
  }
 //      h%x%nhist = inbin
  h->x.nhist = inbin;
  h->xw = (float *) malloc (inbin*sizeof (*h->xw));
 //      allocate( h%xw(inbin) )
  h->dnw = (float *) malloc (inbin*sizeof (*h->dnw)); 
 //      allocate( h%dnw(inbin) )
  h->mean = (float *) malloc (inbin*sizeof (*h->mean));

 //      allocate( h%mean(inbin) )
  h->dndx = (float *) malloc (inbin*sizeof (*h->dndx) );
  h->xforI = (float *) malloc (inbin*sizeof (*h->xforI) );
  h->integ = (float *) malloc (inbin*sizeof (*h->integ) );
 //      allocate( h%dndx(inbin) )


  h->x.tklg= ( itklg - (itklg/2)*2  ) != 0;
 //          ( itklg & 01 ) != 0
  h->x.cent = ( (itklg/2)*2 - (itklg/4)*4 ) /2;
  //           (itklg & 02)  != 0
  h->x.ufl =  ( (itklg/4)*4 - (itklg/8)*8 ) !=0 ;
 //            ( itklg & 04) !=0
  h->x.ofl = ( (itklg/8)*8 - (itklg/16)*16 )  != 0;
 //            (itklg & 05) !=0
  asmax = ( (itklg/16)*16 - (itklg/32)*32 ) != 0;
 //           ( itklg & 06)  !=0
  h->x.xmin = ixmin;      
 //      h%x%xmin = ixmin    !  not used at present
  if( asmax ) {
    if( ixmin >= ibinORxmax ) {
      fprintf(stderr,"ixmin=%f ibinORxmax=%f\n", ixmin, ibinORxmax);
      fprintf(stderr,"ibinORxmax is regarded as ixmax but <= ixmin\n");
      exit (1);
    }
    else {
      if( h->x.cent )  ndiv = inbin - 2;
      else ndiv = inbin-1;
      if (h->x.tklg ) h->x.bin = log10(ibinORxmax/ixmin)/ndiv;
      else       h->x.bin = (ibinORxmax - ixmin )/ndiv;
    }
  }
  else {
    h->x.bin = ibinORxmax;     
  }
   /*
      if(asmax) then
         if(ixmin .ge. ibinORxmax ) then
            write(0,*) ' ibinORxmax is regarded as ixmax but <= ixmin'
            stop 99999
         else
            if( h%x%cent .eq. 1 ) then
               ndiv= inbin - 1
            else
               ndiv = inbin
            endif
            if(h%x%tklg) then
	       h%x%bin = log10(ibinORxmax/ixmin)/ndiv
            else
               h%x%bin = (ibinORxmax - ixmin )/ndiv
            endif
         endif
      else
         h%x%bin = ibinORxmax
      endif
   */
  if( h->x.tklg ) {
    if( h->x.xmin <= 0.0 ) {
      fprintf(stderr, "min must be >0 for log option\n");
      exit(1);
    }
    h->x.xm = log10(h->x.xmin) -
      h->x.cent * h->x.bin/2;
    
    h->x.inc = pow(10., h->x.bin);
  }
  else {
    h->x.xm = h->x.xmin  -  h->x.cent * h->x.bin/2;
    h->x.inc = h->x.bin;
  }
   /* 
   if( h%x%tklg  ) then
         if( h%x%xmin <= 0.0 )  then
            write(0,
     *       '("min must be > 0 for log option")')
            stop
         endif
         h%x%xm = log10(h%x%xmin) - h%x%cent * h%x%bin/2
         h%x%inc = 10.**h%x%bin
      else 
         h%x%xm = h%x%xmin  -  h%x%cent * h%x%bin/2
         h%x%inc  = h%x%bin
      endif
   */

  strcpy(h->c.id, " ");
  h->c.eventno = 1;
    //c      h%c%eventno = 1
  strcpy(h->x.label, "x");

   //      h%x%label =' '

  strcpy(h->x.unit, " ");
   //      h%x%unit = ' '

  strcpy(h->c.title, " ");
   //      h%c%title = ' '
  strcpy(h->c.dNunit, " ");

   //      h%c%dNunit=' '
  strcpy(h->c.categ, " ");
   //      h%c%categ = ' '
  strcpy(h->c.dir, " ");
   //      h%c%dir = ' '
  h->c.pw = 0.;
    //      h%c%pw = 0
  h->c.logv = TRUE;
   //      h%c%logv = .true.
  h->c.norm = - 1.0;
   //      h%c%norm = -1.0
  return;
}

 /*
c    ************************
      entry kwhistc(h)
c    ************************
 */
void  kwhistc(struct histogram1 *h) {
  int i;
  for(i= 0; i< h->x.nhist; i++) {
    *(h->xw+i) = 0.;
    *(h->dnw+i) = 0.;
  }
  return ;
}
/*
  do i = 1, h%x%nhist
         h%xw(i) = 0.
         h%dnw(i) = 0.
      enddo
      return
*/
/*
c    *************************
      entry kwhist( h, x, w )
c    *************************
*/
void kwhist(struct histogram1 *h, float x, float w){
  float xx;
  int i;
  if( h->x.tklg && x <= 0. ){
     // neglect this case
  }
  else {
    if( h->x.tklg  ) xx = log10(x);
    else  xx = x;

    i = ( xx -  h->x.xm ) / h->x.bin  + 1;

    if(i <= 0 && h->x.ufl ) i = 1;
    else if(i >= h->x.nhist && h->x.ofl )  i = h->x.nhist - 1;

     /*
      if( h%x%tklg  .and. x .le. 0.) then
c         neglect this data
      else
         if( h%x%tklg  ) then
            xx = log10(x)
         else
            xx = x
         endif
         i = ( xx-h%x%xm ) / h%x%bin  + 1

         if(i .le. 0 .and. h%x%ufl ) then
            i = 1
         elseif(i .gt. h%x%nhist .and. h%x%ofl ) then
            i = h%x%nhist
	    endif
     */

    if(i >= 1 &&  i < h->x.nhist ) {
      *(h->xw+i-1) += x*w;
      *(h->dnw+i-1) +=  w;
    }
  }

     /*                                                        
         if(i .ge. 1 .and.  i  .le. h%x%nhist )  then
            h%xw(i) = h%xw(i)  +  x*w
            h%dnw(i) = h%dnw(i) + w
         endif
     endif
     */
  return;
}
/*
c     ***********************
      entry kwhists( h, inorm )
c     ************* take statistics
c         if inorm = -1.0, use alredy fixed one.
 */
void kwhists0(int from) {
  fromwhich = from;
  return;
}
void kwhists(struct histogram1 *h, float inorm){
  float xx;
  float dx;
  double integral, xl;
  double isumw;

  if( inorm != -1.0)  h->c.norm = inorm;
  h->x.imin = 1;
  while ( h->x.imin < h->x.nhist &&
	   *(h->dnw + h->x.imin-1) == 0.) {
     h->x.imin++;
   }
   h->x.imax = h->x.nhist;
   while (h->x.imax > 1 &&  *(h->dnw + h->x.imax-1) ==  0.) {
     h->x.imax--;
   }
   
   /*
      if( inorm .ne. -1.0) then
         h%c%norm = inorm
      endif

      h%x%imin = 1
      do while( h%x%imin .lt. h%x%nhist .and.  h%dnw(h%x%imin) .eq. 0.) 
         h%x%imin = h%x%imin + 1
      enddo

      h%x%imax = h%x%nhist
      do while (h%x%nhist .gt. 1 .and.  h%dnw(h%x%imax) .eq.  0.)  
         h%x%imax = h%x%imax -1
      enddo
   */
   h->x.sumw = 0;
   int i;
   isumw = 0.;
   for(i = h->x.imin; i<= h->x.imax; i++){
     isumw += *(h->dnw+i-1);
   }
   h->x.sumw = isumw;

   if(h->c.norm == 0. && h->x.sumw > 0.)  normf = h->x.sumw; 
   else if(h->c.norm <= 0. ) normf = 1.0;
   else  normf = h->c.norm;

   /*
      h%x%sumw = 0
      do i = h%x%imin, h%x%imax
         h%x%sumw = h%x%sumw +  h%dnw(i)
      enddo

      if(h%c%norm .eq. 0. .and.  h%x%sumw  .gt. 0.) then
        normf = h%x%sumw 
      elseif(h%c%norm .le. 0. ) then
         normf = 1.0
      else
         normf = h%c%norm
      endif
   */
   //        bin center value
   if( h->x.tklg ){
     xx = pow(10.0,(h->x.xm + h->x.bin/2.) )
       * pow(h->x.inc,(h->x.imin-1));
     xl = pow(10.0,h->x.xm)  * pow(h->x.inc,(h->x.imin-1));
   }
   else {
     xx = h->x.xm + h->x.bin/2 + h->x.inc*(h->x.imin-1);
     xl = h->x.xm + h->x.inc*(h->x.imin-1);
   }

   /*
c        bin center value
      if( h%x%tklg ) then
         xx =10**(h%x%xm + h%x%bin/2.) * h%x%inc**(h%x%imin-1)
      else
         xx = h%x%xm +   h%x%bin/2 + h%x%inc*(h%x%imin-1)
      endif
   */
   if(fromwhich == 0 )   integral = 0.;
   else integral = isumw;

   dx = h->x.bin;
   for(i = h->x.imin; i<=h->x.imax; i++) {
     if( h->x.tklg ) 
            dx  = pow(10.0, (h->x.xm + i * h->x.bin)) -
	      pow(10.0,(h->x.xm + (i-1)*h->x.bin));

     if( inorm == -1.0 ) {
       // data is probably  read from a file so that mean has been given
       // while xw is not computed we get it here
       if( *(h->dnw +i -1) != 0.)  *(h->xw+i-1) = *(h->mean + i-1) * *(h->dnw+i-1);
       else *(h->xw+i-1) = xx;
     }
     else {
       if(*(h->dnw +i-1) == 0) *(h->mean +i-1) = xx;
       else *(h->mean +i-1) = *(h->xw+i-1)/ *(h->dnw+i-1);
     }
     *(h->xforI +i -1) = xl;
     *(h->integ +i -1) = integral/normf;

     *(h->dndx +i-1)=*(h->dnw+i-1)/dx/normf;
     if( h->x.tklg ) {
       xx = xx * h->x.inc;
       xl = xl * h->x.inc;
     }
     else {
       xx = xx + h->x.inc;
       xl = xl + h->x.inc;
     }
     if(fromwhich == 0)  integral += *(h->dnw+i-1);
     else  integral -= *(h->dnw+i-1);
   }
   int ii=h->x.imax;
   *(h->xforI + ii) =xl;
   *(h->integ + ii) =integral/normf;
   *(h->dndx + ii) =0.;
   *(h->dnw + ii) =0.;
   *(h->mean + ii) =xx;

   /*
      dx = h%x%bin      
      do i = h%x%imin, h%x%imax
         if( h%x%tklg ) then
            dx  = 10.0**(h%x%xm + i * h%x%bin) -
     *            10.0**(h%x%xm + (i-1)*h%x%bin)
         endif
         if(h%dnw(i) .eq. 0) then
            h%mean(i) = xx
         else
            h%mean(i) = h%xw(i)/h%dnw(i)
         endif

         h%dndx(i) = h%dnw(i)/dx/normf
         if( h%x%tklg ) then
            xx = xx * h%x%inc
         else
            xx = xx + h%x%inc
         endif
      enddo
   */
   return;
 }
 void kwhistev(struct histogram1 *h, int evno){
   /*
c     ********************
      entry kwhistev(h, evno)
c     ********************
   */
   h->c.eventno = evno;
   /*
      h%c%eventno = evno
      return
   */
   return;
 }
 void kwhistid(struct histogram1 *h, char *id) {
   /*
c     ********************
      entry kwhistid(h,  id )
c     *******************
   */
   strcpy(h->c.id, id);
   return;
   /*
      h%c%id = id
      return
   */
 }
 void kwhistai(struct  histogram1 *h, char *title, char *categ, char * dNunit,
		int logv, float pw, char * label, char *unit) {
   /*/
c     ********************
      entry kwhistai(h,  title, categ, dNunit, logv, pw, label, unit)
c     *******************
c       additional info.
   */
   strcpy(h->c.title, title);
   strcpy(h->c.categ,  categ);
   strcpy(h->c.dNunit, dNunit);
   h->c.pw = pw;
   h->c.logv = logv;

   strcpy(h->x.label, label);
   strcpy(h->x.unit, unit);
   return;
   /*
      h%c%title = title
      h%c%categ =  categ
      h%c%dNunit = dNunit
      h%c%pw = pw
      h%c%logv = logv
      h%x%label = label
      h%x%unit =  unit
      return
   */
 }
 void kwhistdir(struct histogram1 *h, char *dir) {
   /*
c     *******************
      entry kwhistdir(h, dir)
c     ******************must be called after kwhistai is called
   */
   strcpy(dirstr, dir);
   strcat(h->c.dir,h->c.categ);
   if( strcmp(h->c.categ, "") != 0 ) {
     strcat(h->c.dir,"/");
   }
   strcat(h->c.dir, dirstr);
   return;
 }

   /*
c
c     *********************
      entry kwhistpr( h, fno )
c     ****************print  hist
   */
void  kwhistpr(struct histogram1 *h, FILE *fno) {
  double isumw;

  int i;
  float dx;

  isumw = 0;
  for(i = h->x.imin; i<= h->x.imax; i++){
    isumw += *(h->dnw+i-1);
  }
  h->x.sumw = isumw ;

  if(h->c.norm == 0. && h->x.sumw > 0.)  normf = h->x.sumw; 
  else if(h->c.norm <= 0. ) normf = 1.0;
  else  normf = h->c.norm;




  float xx;

  if( h->x.tklg ) 
    xx = pow(10.0, (h->x.xm + h->x.bin/2.0))
      * pow(h->x.inc, (h->x.imin-1));
  else
    xx = h->x.xm + h->x.bin/2. + h->x.inc*(h->x.imin-1);

      //c        header
  fprintf(fno, "#hist1 %d %d %d %d %d %f\n", 
	  h->c.eventno, h->x.nhist, 
	  h->x.cent, h->x.ufl, h->x.ofl, h->x.bin);
  fprintf(fno, "#t   %s\n", h->c.title);
      //         write(fno, '(a,a)') '#c ', h%c%categ(1:klena(h%c%categ))
  fprintf(fno, "#c   %s\n", h->c.categ); 
   //         write(fno, '(a, a,1x, a)')  '#x ', h%x%label, h%x%unit
  fprintf(fno, "#x  %s   %s\n",  h->x.label, h->x.unit);
     //   write(fno, '(a,f10.2)') '#pw ', h%c%pw
  fprintf(fno, "#pw   %f\n", h->c.pw);
      //         write(fno, '(a,a)') '#dN ', h%c%dNunit(1:klena(h%c%dNunit))
  fprintf(fno, "#dN   %s\n", h->c.dNunit);
     //         write(fno, '(a,a)') '#k ', h%c%id(3:klena(h%c%id))
  fprintf(fno, "#k  %s\n", h->c.id);
  fprintf(fno, "#l  %d  %d\n", h->x.tklg, h->c.logv);
  fprintf(fno, "#n  %11.3e  %12.5e\n", isumw, normf);
  fprintf(fno, "#o %d %d %12.5e %12.5e\n",
	  h->x.imin, h->x.imax,
	  h->x.xm, h->x.inc);
      
     /*
      isumw = h%x%sumw

      if( h%x%tklg ) then
         xx = 10.0**(h%x%xm + h%x%bin/2.0) * h%x%inc**(h%x%imin-1)
      else
         xx = h%x%xm + h%x%bin/2. + h%x%inc*(h%x%imin-1)
      endif
c        header
      if(fno .lt. 0) then
         write(*, '(a, i3)') '#hist1 ', h%c%eventno
         write(*, '(a,a)') '#t ', h%c%title(1:klena(h%c%title))
         write(*, '(a,a)') '#c ', h%c%categ(1:klena(h%c%categ))

         write(*, '(a, a,1x, a)')  '#x ', h%x%label, h%x%unit
         write(*, '(a,f10.2)') '#pw ', h%c%pw

         write(*, '(a,a)') '#dN ', h%c%dNunit(1:klena(h%c%dNunit))
         write(*, '(a,a)') '#k ', h%c%id(1:klena(h%c%id))
         itempx  = 0
         if( h%x%tklg )  itempx = 1
         itempv = 0
         if( h%c%logv) itempv = 1
         write(*,'(a, 2i3)')'#l ',  itempx, itempv
         write(*, '(a, 1p2E15.6)')'#n ', isumw, normf
      else
         write(fno, '(a,i3)') '#hist1 ', h%c%eventno
         write(fno, '(a,a)') '#t ', h%c%title(1:klena(h%c%title))
         write(fno, '(a,a)') '#c ', h%c%categ(1:klena(h%c%categ))

         write(fno, '(a, a,1x, a)')  '#x ', h%x%label, h%x%unit
         write(fno, '(a,f10.2)') '#pw ', h%c%pw

         write(fno, '(a,a)') '#dN ', h%c%dNunit(1:klena(h%c%dNunit))
         write(fno, '(a,a)') '#k ', h%c%id(3:klena(h%c%id))
         itempx = 0
         if( h%x%tklg )  itempx = 1
         itempv = 0
         if( h%c%logv) itempv = 1
         write(fno,'(a,2i3)') '#l ',  itempx, itempv
         write(fno, '(a, 1p2E15.6)') '#n ', isumw, normf
      endif
     */

     //c          directory where data is saved
     //dirstr= h%c%dir(1:klena(h%c%dir))
  strcpy(dirstr, h->c.dir);
  fprintf(fno, "#d %s\n", dirstr);
     /*
      if(fno .lt. 0 ) then
         write(*,'(a, a)') '#d ', dirstr(1:nstr)
      else
         write(fno,'(a, a)') '#d ', dirstr(1:nstr)
      endif
     */


     //      do i = h%x%imin, h%x%imax
      

  for(i = h->x.imin; i<=h->x.imax; i++) {

    fprintf(fno, "%5d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n", 
	    i, xx, *(h->dndx+i-1), *(h->dnw+i-1), 
	    *(h->mean+i-1), *(h->xforI+i-1), *(h->integ+i-1) );

       //         isumw = isumw -  h%dnw(i)
       /*       
         if( h%x%tklg ) then
            xx =  xx * h%x%inc 
         else
            xx = xx + h%x%inc
         endif
       */
    if(h->x.tklg) xx *=h->x.inc;
    else xx += h->x.inc;
       //      enddo
  }
  fprintf(fno, "%5d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n", 
	h->x.imax+1, xx, 0., 0.,  xx,  *(h->xforI + h->x.imax),
	  *(h->integ + h->x.imax) );
     /*
c       trailer
      if(fno .lt. 0) then
         write(*,'(7i3)')  0,0,0,0,0,0, 0
      else
         write(fno,'(7i3)')  0,0,0,0,0,0, 0
      endif
     */
  fprintf(fno, "0  0  0  0  0  0  0\n");

     //     return
  return;
}
     /*
c     *********************
      entry kwhistw(h, bfnow)
c     ********************      
c       binary write of h to bfnow
     */
 void kwhistw(struct histogram1  *h, FILE *bfnow){
       /*
      write(bfnow) '#hist1'
      write(bfnow) h%x%nhist
      write(bfnow) h%x, h%c
      write(bfnow) h%xw, h%dnw, h%mean, h%dndx
       */
   size_t size;
   size=sizeof("#hist1");

   fwrite( "#hist1",  1, size,  bfnow ); 


   fwrite(&h->x.nhist, sizeof(short int), 1, bfnow);

   int nbin;
   nbin = h->x.nhist;

   fwrite( &h->x, sizeof(h->x), 1,  bfnow);
   fwrite( &h->c, sizeof(h->c), 1,  bfnow);

   fwrite( h->xw, sizeof(*h->xw), nbin,  bfnow);
   fwrite( h->dnw, sizeof(*h->dnw), nbin,  bfnow);
   fwrite( h->mean, sizeof(*h->mean), nbin,  bfnow);
   fwrite( h->dndx, sizeof(*h->dndx), nbin,  bfnow);
   fwrite( h->xforI, sizeof(*h->xforI), nbin,  bfnow);
   fwrite( h->integ, sizeof(*h->integ), nbin,  bfnow);

   return;
 }
 /*
c     *********************
      entry kwhistr(h, bfnor, icon)
c     ********************
c        #hist1 must be read outside
c
 */
 void kwhistr(struct histogram1 *h, FILE *bfnor, int icon){
   //      read(bfnor, end =222)  nbin 
   short int nbin;
   fread(&nbin, sizeof(short int), 1, bfnor);
   //      allocate( h%xw(nbin) )
   h->xw = (float *) malloc (nbin*sizeof (*h->xw));
   // allocate( h%dnw(nbin) )
   h->dnw = (float *) malloc (nbin*sizeof (*h->dnw));
   //      allocate( h%mean(nbin) )
   h->mean = (float *) malloc (nbin*sizeof (*h->mean));
   //      allocate( h%dndx(nbin) )
   h->dndx = (float *) malloc (nbin*sizeof (*h->dndx));
   h->xforI = (float *) malloc (nbin*sizeof (*h->xforI) );
   h->integ = (float *) malloc (nbin*sizeof (*h->integ) );

   //      read(bfnor, end=222)  h%x, h%c
   fread(&h->x, sizeof(h->x), 1, bfnor);

   fread(&h->c, sizeof(h->c), 1, bfnor);
   //      read(bfnor, end= 222)  h%xw, h%dnw, h%mean, h%dndx
   fread(h->xw, sizeof(*h->xw), nbin, bfnor);
   fread(h->dnw, sizeof(*h->dnw), nbin, bfnor);
   fread(h->mean, sizeof(*h->mean), nbin, bfnor);
   fread(h->dndx, sizeof(*h->dndx), nbin, bfnor);
   fread(h->xforI, sizeof(*h->xforI), nbin, bfnor);
   fread(h->integ, sizeof(*h->integ), nbin, bfnor);

   //icon = 0
   icon =0;
   //return
   return;
   /*
 222  continue
      write(0,*) ' kwhistr reached EOF unexpectedly'
      icon = 1
      return
   */
 }
 /*
c     ****************
      entry kwhistd(h)
c     ***************
c      deallocate histogram area
c
 */
 void kwhistd(struct histogram1 *h){
   //      h%c%init = ' '
   strcpy(h->c.init, " ");
   //      deallocate( h%xw, h%dnw,  h%mean,  h%dndx, stat=dealloc)
   free(h->xw);
   free(h->dnw);
   free(h->mean);
   free(h->dndx);
   free(h->xforI);
   free(h->integ);
   /*
      if(dealloc .ne. 0) then
         write(0,*) ' dealloc failed =',dealloc
         stop 12345
      endif
      return
   */
   return;
 }

 /*
c     ********************
      entry kwhista(h1, h2, h)
c     ******************
c      h = h1 + h2  of bin area. For others, h1 is inherited
c      h,  h1 and h2 must be the same size  histogram  of same 
c     type.  h can be h1
c
 */
 void kwhista(struct histogram1 *h1, struct histogram1 *h2, struct histogram1 *h){
   /*
      if( h1%x%nhist .ne. h2%x%nhist) then
         write(0, *)
     *    ' h1 and h2 diff. size histogram in kwhista'
         stop 9876
      endif
   */
   short int nbin;
   int i;

   if( h1->x.nhist != h2->x.nhist ){
     fprintf(stderr,"h1 != h2 in kwhista\n");
     exit;
   }

   /*
      if( h%c%init .ne. 'initend') then
c           not yet initialized.
         nbin = h1%x%nhist
         allocate( h%xw(nbin) )
         allocate( h%dnw(nbin) )
         allocate( h%mean(nbin) )
         allocate( h%dndx(nbin) )
         h%c%init = 'initend'
      endif
   */
   if( strcmp(h->c.init, "initend") !=0) {
     nbin=h1->x.nhist;
     h->xw = (float *) malloc(sizeof(*h->xw) *nbin);
     h->dnw = (float *) malloc(sizeof(*h->dnw) *nbin);
     h->mean = (float *) malloc(sizeof(*h->mean) *nbin);
     h->dndx = (float *) malloc(sizeof(*h->dndx) *nbin);
     h->xforI = (float *) malloc(sizeof(*h->xforI) *nbin);
     h->integ = (float *) malloc(sizeof(*h->integ) *nbin);

     strcpy(h->c.init, "initend");
   }
   /*
      h%x = h1%x
      h%c = h1%c
   */
   h->x = h1->x;
   h->c = h1->c;
   /*
      do i = 1, h%x%nhist
         h%xw(i) = h1%xw(i) + h2%xw(i)
         h%dnw(i) = h1%dnw(i) + h2%dnw(i)
      enddo
   */
   for(i=0; i< h->x.nhist; i++){
     *(h->xw+i) = *(h1->xw+i) + *(h2->xw+i);
     *(h->dnw+i) = *(h1->dnw+i) + *(h2->dnw+i);
   }
   //  end
   return;
 }
 /*
      subroutine kwhistso( binw )
c        specify output method
      implicit none
      include "Z90histCom.h"
      integer binw  ! input.  1--> ascii write
                    !         2--> binary write
 */
 void kwhistso(int binw){

   //      BinWrite = binw
   BinWrite = binw;
   /*
      if(binw .ne. 1 .and. binw .ne. 2) then
         write(0,*) 'binw=',binw,' for kwhistso is invalid'
         stop
      endif
   */
   if(binw != 1 && binw != 2 ) {
     fprintf(stderr, "binw=%d is invalid for kwhistso\n", binw);
     exit(1);
   }
       //      end
 }
 /*
      subroutine kwhistp( h, fno )
c
c         print or binary write histogram
c         kwhistso must be called to
c         fix binary write or print
c
      implicit none
      include 'Z90histc.h'
      include 'Z90histo.h'
      include "Z90hist1.h"
      include "Z90histCom.h"
      type(histogram1) h
      integer fno

 */
 void kwhistp(struct histogram1 *h, FILE *fno){
   /*

      if( BinWrite .eq. 2 ) then
         call kwhistw(h, fno)
      elseif(BinWrite .eq. 1 ) then
         call kwhistpr(h, fno)
      endif
      end
   */

   if( BinWrite == 2 )  kwhistw( h, fno);
   else  kwhistpr( h, fno);
   return;
 }
int kwhistIxy(struct histogram1 *h, double x[],
		double y[],int n) {
  int m;
  int i;

  m = h->x.imax-h->x.imin + 2;
  if( m > n)    m = -1;
  else {
    m = 0;
    for( i=h->x.imin; i<=h->x.imax+1; i++) {
      m++;
      x[m-1] = *(h->xforI +i-1);
      y[m-1] = *(h->integ +i-1);
    }
  }
  return m;
}
int kwhistxy(struct histogram1 *h, double x[],
		double y[], int n) {
  int m;
  int i;

  m = h->x.imax-h->x.imin + 2;
  if( m > n)    m = -1;
  else {
    m = 0;
    for( i=h->x.imin; i<=h->x.imax+1; i++) {
      m++;
      x[m-1] = *(h->xw +i-1);
      y[m-1] = *(h->dndx +i-1);
    }
  }
  return m;
}
