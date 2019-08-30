/*
  make some ascii histogram data file: 

  For that test purpose, you can use
  test1.c which generate 3 histograms of power spectra
  and 10 histograms by  gaussian distributions.
  you can create it as ascii or binary file.
  In case of binary file, you can use bin2ascii.c
  to convert it to ascii file.
  For usage see Readme.
*/

#include <stdio.h>
#include <stdlib.h>
#include "Z90histo.h"
#include "Z90histc.h"
#include "Z90hist1.h"
#include "Z90histfuncdef.h"

int main(int argc, char *argv[]){
  struct histogram1 h;
  FILE *fno, *fni;
  if(argc != 2 ) {
    fprintf(stderr, "Usage: %s <input_asciifile  output_bin_file_name\n",
	    argv[0]);
    exit(1);
  }
 
  /*
  fni=fopen("mytest2.ahist", "r");
  if(fni == NULL ){
    fprintf(stderr, "cannot open input\n");
    exit(1);
  }
  */
  fni=stdin;

  fno = fopen( argv[1], "wb");
  if(fno == NULL) {
    fprintf(stderr, "output file  %s cannot be opened\n",
	    argv[1]);
  }

  int icon;
  icon = 0;
  kwhistso(2);  //  binary output
  while (( icon = kwhistReadAscii( &h, fni) )>=0) {
    fprintf(stderr,"None zero hist rec=%d\n",icon);
    kwhistw( &h, fno);
    kwhistd( &h);
  }
  fprintf(stderr, "last icon=%d\n", icon);
}
