extern "C" {
  //#include <sstream>
#include "Zlogical.h"
#include "Zprivate.h"
#include "Ztrack.h"
#include <string.h>
#include <string.h>
#include <iostream.h>
#include <strstream.h>

  //extern int klena_(char *, int);
  // #include "cputNull.cc"

  extern void ksplit_(char *, int *,  int *, char *, int *, int, int);

  int main(){
    
    integer  i,  nr, ndev, withir, nev;
    struct memo {
      integer num;
      integer cumnum;
      integer irevent[2];
      struct track Zfirst;
    } memoinst;

    char str[160], msg[160];
    char skelin[80], stro[4][80], numbin[80];


    integer nlow;

    //    Mdev = 71;
    //    Rdev = 70;
    //    ndev = 69;


    strcat(str, " ");
    cin >> str;

    for(i = 0; i<=3; i++)  strcpy(&stro[i][0], " ");

    int temp1 =80;
    int temp2 =4;

    ksplit_(str, &temp1, &temp2, stro, &nr, 160, 80);
    //    call ksplit(str,  80, 4, stro,  nr)


    if( nr < 4 ) {
      cerr << " must give 3 files and flag " << endl;
      exit (1);
    };
    strncpy(seklin, &stro[0][0], 80 );    //      skelin = stro(1)
    strncpy(numbin, &stro[1][0], 80 );    //      numbin = stro(2)
    /*
       // read( stro(3), * )  withir  
       this is realized by the next:
    istringstream iss(stro[2]);
    iss >> withir ; //this is the same as  sscanf(stro[2], "%d", &withir);
    next is also the same as above;
    */

    istrstream iss(stro[2], sizeof stro[2]);
    iss >> withir;

    strncpy(Mskel, &stro[3][0], 80);     // Mskel = stro(4)

    //    cputNull(Mskel);
    //    cputNull(skelin);
    //    cputNull(numbin);
    
    ostrstream oss(msg, sizeof msg);

    fstream Rdev;
    fstream ndev;
    fstream Mdev;
    Rdev.open(skelin, ios::in  | ios::binary); 
    ndev.open(numbin, ios::in  | ios::binary);
    Mdev.open(Mskel,  ios::out | ios::binary);
    
    while( 1 ){
      if(withir == 0) {
	ndev.read( (char *) (&nev), sizeof(nev) );
	  //	  read(ndev, *, end= 1000) nev;
      }
      else {
	ndev.read( (char *) (irevent), sizeof(irevent));
	ndev.read( (char *) (&nev), sizeof(nev));
	//     read(ndev, *, end= 1000) irevent, nev;
      }
      if( ! ndev.fail()) {
	Copy =false;
	while ( !  Copy) {
	  Rdev.read( (char *) (& memoinst), sizeof(memoinst));
	  if(memoinst.cumnum >  nev) {
	    oss << "specified event #=" << nev << 
	      " not exist in skeleton-node file" << ends;
	    cerr << msg << endl;
	    oss << "skip to event=" << cumnum << ends;
	    cerr << msg << endl;
	    while( nev <  memoisnt.cumnum ) {
	      if( withir ==  0) ndev.read( (char *)(&nev), sizeof(nev));
	      else {
		ndev.read( (char *) irevent, sizeof irevent);
		ndev.read( (char *) (&nev), sizeof nev);
	      }
	      if(ndev.fail()) {
		cerr << " all events copied " <<  endl;
		exit (0);
	      }
	    }
	  }

	  if( memoinst.cumnum == nev ) Copy = true; 
	  else Copy = false;

	  if(Copy) Mdev.write( (char *) (&memoinst), sizeof memoinst);

	  fstream *pdev;
	  pdev =  &Rdev;
	  cgetHES(pdev);
	  nlow = 1;
	  while (nlow != -1) {
	    Rdev.read( (char *) (&nlow), sizeof nlow);
	    Rdev.read( (char *) (&p), sizeof p);
	    if(Copy) {
	      Mdev.write( (char *) (&nlow), sizeof nlow );
	      Mdev.write((char *) (&p), sizeof p);
	    }
	    for(i = 0; i< nlow; i++) {
	      Rdev.read( (char *) (&c), sizeof c);
	      if(Copy) Mdev.write((char *) (&c), sizeof c);
	    }
	  }
	}
	cerr << " all events copied" << endl; 
      }
    }
  }

  void cgetHES(fstream *from) {
    integer i;
    from.read( (char *)(&Np), sizeof Np);
    if(Copy)  Mdev.write((char *)(&Np), sizeof Np);
    for(int i = 0; i <  Np; i++) {
      from.read( (char *) (&o[i]), sizeof o[i]);
      if(Copy) Mdev.write((char *) (&o[i]), sizeof o[i]);
    }
  }

}
