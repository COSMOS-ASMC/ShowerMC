#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Z90histfuncdef.h"
#include <math.h>
#define MAXLENG 121  
int  kwhistReadAscii(struct histogram1 *h, FILE *fno) {
  /*  h and fno are input
      return value: 
      -1: requested data not exist
      >=0: some data read. the number of bins
      with non zero data is given
  */
        //c        header
  char string[MAXLENG];
  char head[20];
  char temp[20];

  static  int MAXL=MAXLENG-1;
  int nbin;
  int i;
  float xx, xx2, xx3;
  float tt, tt1, tt2;
  int icon, jj, headl;
  
  icon = -1;
  
  while ( fgets(string, MAXL, fno) != NULL ) {
    jj=strlen(string);
    //  string[jj-1] has \n ; replace it by \0
    strcpy( &string[jj-1], "\0");
    
    sscanf(string, "%s", &head);
    headl = strlen(head);
    /////
    //    fprintf(stderr,"%s head=%s\n", string, head);
    ///////

    if( strcmp(head, "#hist1") == 0 ) {
      icon = 0;
      sscanf(string, "%s %d %d %d %d %d %f",
	     &temp, &h->c.eventno, &h->x.nhist,
	     &h->x.cent, &h->x.ufl, &h->x.ofl, &h->x.bin);
      /////////
      //      fprintf(stderr, " %d %d %d %d %d %f\n",
      //	     h->c.eventno, h->x.nhist,
      //     h->x.cent, h->x.ufl, h->x.ofl, h->x.bin);
      //////////
	      
      if( strcmp(h->c.init, "initend") !=  0 ) {
	nbin = h->x.nhist;
	h->xw = (float *) malloc (nbin*sizeof (*h->xw));
	h->dnw = (float *) malloc (nbin*sizeof (*h->dnw));
	h->mean = (float *) malloc (nbin*sizeof (*h->mean));
	h->dndx = (float *) malloc (nbin*sizeof (*h->dndx));
	h->xforI = (float *) malloc (nbin*sizeof (*h->xforI));
	h->integ = (float *) malloc (nbin*sizeof (*h->integ));
	// clear array
	kwhistc( h );

	strcpy(h->c.init, "initend");
      }
    }
    //  fprintf(fno, "#t   %s\n", h->c.title);
    else if( strcmp(head, "#t") == 0 ) {
      strcpy( h->c.title, "");
      strcpy( h->c.title, &string[headl]);
    //      sscanf(string,"%s %s", &temp,  &h->c.title);
      /////
      // fprintf(stderr, "%s\n", h->c.title);
      //////////
    }
    else if(strcmp(head, "#c") == 0 ) {
      strcpy( h->c.categ, "");
      strcpy( h->c.categ, &string[headl]);
      /////////
      // fprintf(stderr,"%s\n", h->c.categ);
      //////////////
    }
    //  fprintf(fno, "#x  %s   %s\n",  h->x.label, h->x.unit);
    else if(strcmp(head,"#x") == 0 ) {
      strcpy(h->x.unit,"");
      sscanf(string, "%s %s %s", &temp, 
	     &h->x.label, &h->x.unit);
      //////
      // fprintf(stderr, "%s %s \n",
      //     h->x.label, h->x.unit);
      ///////////
    }
    //  fprintf(fno, "#pw   %f\n", h->c.pw);
    else if(strcmp(head,"#pw") == 0) {
      sscanf(string, "%s %f",&temp,  &h->c.pw);
      //////////////
      // fprintf(stderr, "%f \n",
      //	      h->c.pw);
      ///////////
    }
    //  fprintf(fno, "#dN   %s\n", h->c.dNunit);
    else if(strcmp(head,"#dN") == 0) {
      strcpy( h->c.dNunit, "");
      strcpy( h->c.dNunit, &string[headl]);
      //////////
      // fprintf(stderr, "%s \n",
      //	      h->c.dNunit);
      ///////////
    }
    //fprintf(fno, "#k  %s\n", h->c.id);
    else if(strcmp(head,"#k") == 0 ){
      strcpy( h->c.id, "");
      strcpy( h->c.id, &string[headl]);
      //////////
      // fprintf(stderr, "%s \n",
      //      h->c.id);
      ///////////
    }
    //fprintf(fno, "#l  %3d  %3d\n", itempx, itempv);
    else if(strcmp(head,"#l") == 0 ) {
      sscanf(string,"%s %d %d", &temp, &h->x.tklg, &h->c.logv);
      //////////
      // fprintf(stderr, "%d %d \n",
      //      h->x.tklg, h->c.logv);
      ///////////
    }
    //fprintf(fno, "#n  %11.3e  %11.3e\n", isumw, normf);
    else if(strcmp(head, "#n") == 0) {
      sscanf(string,"%s %f %f", &temp, &h->x.sumw, &h->c.norm);
      //////////
      // fprintf(stderr, "%f %f \n",
      //	      h->x.sumw, h->c.norm);
      ///////////
    }
    else if(strcmp(head, "#o") == 0) {
      sscanf(string, "%s %d %d %f %f",
	&temp, &h->x.imin, &h->x.imax,
	     &h->x.xm, &h->x.inc);
      //////////
      // fprintf(stderr, "%d %d %f %f\n",
      //	 h->x.imin, h->x.imax,
      //     h->x.xm, h->x.inc);
      ///////////

    }
    //  strcpy(dirstr, h->c.dir);
    //  fprintf(fno, "#d %s\n", dirstr);
    else if(strcmp(head, "#d") == 0 ) {
      strcpy( h->c.dir, "");
      strcpy( h->c.dir, &string[headl]);
      /////
      // fprintf(stderr, "%sf\n",
      //      h->c.dir);
      ///////////
    }
    /*
      dx = h->x.bin;
      for(i = h->x.imin; i<=h->x.imax; i++) {
      if(h->x.tklg ) dx = pow(10.0, (h->x.xm + i*h->x.bin)) 
			- pow(10.0, (h->x.xm + (i-1)*h->x.bin));
			fprintf(fno, "%5d %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e\n",   i, xx, *(h->dndx+i-1), *(h->dnw+i-1), isumw/normf,
	     *(h->mean+i-1), dx);
	     */
    else if( strcmp(head, "0") == 0 ) {
      break;
    }
    else { 
      sscanf(string, "%d  %e  %e %e %e %e %e",
	     &i, &xx, &tt, &tt1, &tt2, &xx2, &xx3);
      if( i <= h->x.imax+1 ) { 
	*(h->dndx+i-1)=tt;
	*(h->dnw+i-1) = tt1;
	*(h->mean+i-1)= tt2;
	*(h->xforI+i-1) =xx2;
	*(h->integ+i-1) =xx3;
	icon ++;
      }
      ///// next one dose not work in sscanf:
     //	     &i, &xx, (h->dndx+i-1), (h->dnw+i-1), &xx2,
     //             (h->mean+i-1), &xx3);
      /////////////
      //fprintf(stderr, "%d  %e  %e %e %e %e %e\n",
      //      i, xx, *(h->dndx+i-1), *(h->dnw+i-1), xx2,
      //     *(h->mean+i-1), xx3);
      ///////////
    }
  }
  return icon;
}
