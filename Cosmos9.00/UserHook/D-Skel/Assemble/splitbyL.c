/* FDD source data is asuumed in some directory (*.dat files)
 * This reads  *.dat split *.dat files into xxx-Lyy.dat
 * where yy's are the FDD layer number.
 */
#include <stdlib.h>
#include <stdio.h>
char datfile[150];
int  nterm = 14;  // make this 13 for old FDD source data
#define  MaxLayers 70
struct gfile myfile[MaxLayers];
int ifopened[MaxLayers];

//  FILE *fdnewdat[MaxLayers][Nk];
FILE *fdnewdat[MaxLayers];
FILE *fddat;

void splitDatAndWrite();
void closeNewDatFiles();
void openDatFile();
void openLfile(int l);

struct gfile {
  char name[250];
};

  char inputdir[100];
  char base[50];
  char outputdir[100];


  int i;
  int ly, code, nr, nf;

  int  code, subcode, charge;
  int ly, ridx, faiidx;
  double r, fai, E, tm, w1, w2, w3;


  double E0;


int main(int argc, char *argv[]){
  fprintf(stderr, "checking command line arg\n");
  fprintf(stderr, "You gave %d command line args\n", argc-1); 
  if( argc !=  4 ) { 
    fprintf(stderr, "Usage: \n");
    fprintf(stderr, "  splitByL  inputfilePath file-body  output-dir \n");
    fprintf(stderr, "file-body is xxxx in xxxx.hyb, xxxx.nrfai\n");
    fprintf(stderr, "ouput-dir must have  last /\n");
    return 1;
  }
  fprintf(stderr, "getting input filepath\n");
  strcpy(datfile, argv[1]);
  fprintf(stderr, "dat file is %s\n", datfile);
  fprintf(stderr, "getting base file name\n");
  strcpy(base, argv[2]);
  fprintf(stderr, "base file name is %s\n", base);

  fprintf(stderr, "getting output dir\n");
  strcpy(outputdir, argv[3]);
  fprintf(stderr, "output dir is %s\n", outputdir);
 
 
    
  fprintf(stderr, " input comannd line arg. processing end\n");
  // form file names
  for(i=0; i<MaxLayers;  i++) ifopened[i]=0;

  openDatFile();

  fprintf(stderr, " openfiles end\n");

  splitDatAndWrite();
  fprintf(stderr, " split  end\n");
  //
   closeNewDatFiles();
  //----------
   fprintf(stderr, " closenewdat end\n");

}	    


/*    *****************************           */
 void openDatFile(){
  //  open files
  fddat=fopen(datfile, "r");
  if( fddat == NULL ) {
    fprintf(stderr, "open err for datfilen");
    exit(1);
  }
 }
void openLfile(int ly) {
  sprintf(myfile[ly-1].name,
	  "%s%s-L%2.2d.dat",outputdir,  base, ly);
  fdnewdat[ly-1]=fopen(myfile[ly-1].name,  "a");
 }
 /*        ********************* */

 void splitDatAndWrite(){
  // read .dat and convert
  double fais;
  int codef;
  double dummy;
#define  LINESIZE 200  
  char val[LINESIZE];
  while( fgets(val, sizeof(val), fddat) !=NULL) {

    //    nterm = kcountTerms(val);

    // in old FDD no 14th term (dummy below).
    // new one has weight there to be able to use
    // for thinnined FDD. (not used yet).  
    if(nterm == 14) {
       sscanf(val,
	      "%d   %d %d %d   %d %d  %lf %lf %lf %lf %lf %lf %lf %lf",
	      &ly, &code, &subcode, &charge, &ridx, &faiidx,
	      &r, &fai, &E, &tm, &w1, &w2, &w3, &dummy);
       
    }
    else {
       sscanf(val,
	      "%d   %d %d %d   %d %d  %lf %lf %lf %lf %lf %lf %lf",
	      &ly, &code, &subcode, &charge, &ridx, &faiidx,
	      &r, &fai, &E, &tm, &w1, &w2, &w3);
    }
   fais = fai;
    //////////
    //fprintf(stderr, "fais=%f  faiidx=%d\n",fais, faiidx);
    ///////////
    // update nrfai 

    if(code > 4) codef=4;
    else codef = code;

    //  1  1  3  0 21  2  9.837E-01  16.1  7.343E-04  6.311E+02 -0.3136 -0.7232  0.615272
    if( ifopened[ly-1] == 0 ) {
      ifopened[ly-1] = 1;
      openLfile(ly);
      fprintf(stderr, "ly =%d opened\n", ly);
    }

    fprintf(fdnewdat[ly-1],
	    "%3d%3d%3d%3d%3d%3d%11.3e%6.1f%11.3e%11.3e%9.4f%9.4f%9.6f\n",
	    ly, code, subcode, charge, ridx, faiidx, r, fais, E, tm,
	    w1, w2, w3);

  }
 }
 void closeNewDatFiles(){
       // close  .data-new
  for(ly=0; ly<MaxLayers; ly++) {
    if( ifopened[ly] == 1 )  fclose(fdnewdat[ly]);
  }
 }




