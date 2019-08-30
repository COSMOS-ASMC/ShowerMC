//  define next if cosmos library in $TAMCDB/F/lib/
//  is used for the random number generator.
//  If USECOSMOS is not defined, gnu random number 
//  generators are used.  gnu routines are assumed
//  to be in /usr/local/ 
//
//  see Readme

#undef   USECOSMOS


/* ******************************************* */ 

#include <stdio.h>
#include <stdlib.h>

#ifndef USECOSMOS
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

//#include <config.h>
#include <math.h>
#include "Z90histfuncdef.h"
//#include "Z90histc.h"
//#include "Z90histo.h"
//#include "Z90hist1.h"
int main(int argc, char *argv[]){
    //      integer nav, nsig, npw
    // parameter (nav=2, nsig=5, npw=3 )
#define NAV 2
#define NSIG 5
#define NPW 3
  //  type(histogram1)  h(nav, nsig)
  //  type(histogram1)  k(npw)
  //      save h
static  struct histogram1 h[NSIG][NAV];
struct  histogram1 k[NPW];

  //      real*8 av, sig, pw
  double av, sig, pw;
  //      real*8 x
  double x;

  //integer i, j, m, fno
  int i, j, m;
  FILE *fno;
  //      integer nc
  int nc;

  //      character*48 dirstr
  //    character*38 key
  char dirstr[48];
  char key[38];
  
  if( argc != 3 ) {
    fprintf(stderr,
	    "Usage: %s  asciiOrbin(1 or 2) and output file name\n", argv[0]);
    exit(1);
  }
  int binorascii;
  binorascii=atoi(argv[1]);

  fno=fopen(argv[2], "wb");
  if(fno == NULL ) {
    fprintf(stderr, "ouput file=%s cannot be opened\n",
	    argv[2]);
  }
  //***********************gsl
#ifndef USECOSMOS
  const gsl_rng_type *T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T=gsl_rng_default;
  r=gsl_rng_alloc(T);
#endif

  kwhistso(binorascii);
  /*
c       minimum calls
      pw = 0.8
  */
  pw = 0.8;

  //      do i = 1, npw
  for(i=1; i<=NPW; i++) {
    //         pw = pw + 0.2
    pw +=0.2;
    /*
c              init.
         call kwhisti(k(i), 1.5, 0.1, 30, b'01111' )
    */
    kwhisti( &k[i-1], (float) 1.5, (float ) 0.1, 30, 017);

    //c              clear
    //     call kwhistc(k(i))
    kwhistc(&k[i-1]);
    //         do j= 1, 1000000
    //
    for(j=1; j<=1000000; j++) {
      double x;
      //call rndc(x)
#ifdef USECOSMOS
      rndc_(&x);
#else
      //x = gsl_ran_flat(r, 0., 1.);
      x = gsl_rng_uniform(r);
#endif
      //            x = x**(-pw)
      x = pow(x, -pw);
      /*
c               take histo
            call kwhist( k(i), sngl(x), 1.0 )
      */
      kwhist(&k[i-1], (float) x, (float) 1.0);
      //    enddo
    }
    //enddo
  }
  /*
c         output 
      do i =1, npw
  */
  for(i=1; i<=NPW; i++) {
    /*
         call kwhists(k(i), 0.)
         call kwhistp(k(i), fno )
    */
    kwhists(&k[i-1], (float)0.);
    kwhistp(&k[i-1], fno);
    //   enddo
  }
  ///////////
  fprintf(stderr, "end of power spec. write\n");

  /*
c ++++++++++++++++++++++++++++++
c    some standard
c
  */
  //  av=0.
  av = 0.;
  //      do i = 1, nav
  for(i=1; i<=NAV; i++) {
    //av = av + 2. 
    av +=2.;
    // sig = 0.2 
    sig=0.2;
    //         do j = 1, nsig
    for(j=1; j<=NSIG; j++){
      //sig = sig + 0.2
      sig += 0.2;
      /*
c             init.
            call kwhisti(h(i, j), sngl(av-7*sig), sngl(av+7*sig),  
     *       500,  b'10000')
     */
      kwhisti(&h[j-1][i-1], (float)(av-7.0*sig), (float) (av+7.0*sig),
	      500, 020);
      /*
c             clear
            call kwhistc(h(i, j))
      */
      kwhistc(&h[j-1][i-1]);

      //            do m = 1, 1000000
      for(m=1; m<= 1000000; m++){

	//     call kgauss(av, sig, x)
#ifdef USECOSMOS
	kgauss_(&av, &sig, &x);
#else
	x =gsl_ran_gaussian(r, sig) + av;
#endif
	/*
c                 take histo
               call kwhist(h(i,j), sngl(x), 1.0 )
	*/
	kwhist( &h[j-1][i-1], (float) x, (float) 1.0);

	//            enddo
      }
      /*
c             give additional info.
            call kwhistai(h(i,j), 
     *           "Test Gaussian dist.",
     *           "gauss", "event", .false., 0., 
     *           "x", "m") 
     */
      kwhistai( &h[j-1][i-1], "Test Gaussian dist.",
	         "gauss", "event", 0,  (float) 0., 
                "x", "m") ;
      /*
c             make key for diff. parameers
            write(key,'(1p2E11.3)' ) av, sig
      */
      sprintf(key, "%11.3e %11.3e", av, sig);
      /*      
c              inform it
            call kwhistid(h(i,j), key)
      */
      kwhistid(&h[j-1][i-1], key);
      /*
c              make directory: maindir/gauss/{av1,av2}
            write(dirstr,'("av",i2,"/")')  i
      */

      sprintf(dirstr, "av%d/", i);
      /*
c                   next two is better to shrink the string length
c                   but only 3rd line can be  ok. (white blank will 
c                   be eliminated inside).
c            call kseblk(dirstr,"|", nc)
c            call kwhistdir(h(i,j), dirstr(1:nc))
            call kwhistdir(h(i,j), dirstr)
      */
      kwhistdir( &h[j-1][i-1], dirstr);
      //         enddo
    }
    // enddo
  }
  /*
c         output.
      do i = 1, nav
  */
  for(i=1; i<=NAV; i++){
    //         do j = 1, nsig
    for(j=1; j<=NSIG; j++){
      //            call kwhists(h(i,j), 0.)
      //    call kwhistp(h(i,j), fno)
      kwhists(&h[j-1][i-1], (float)0.);
      kwhistp(&h[j-1][i-1], fno);
      ///// just test ***
      if( i==1 && j==2 ){
	int nnn, kk;
	double xxx[500],  yyy[500];
	nnn=kwhistIxy(&h[j-1][i-1], xxx, yyy, 500); 
	fprintf(stderr,"nnn=%d\n", nnn);
	for(kk=0;kk<nnn; kk++) {
	  fprintf(stderr, "%f %f\n", xxx[kk], yyy[kk]);
	}
      }
      //////////***
      //         enddo
    }
    // enddo
  }
  //   end
}

