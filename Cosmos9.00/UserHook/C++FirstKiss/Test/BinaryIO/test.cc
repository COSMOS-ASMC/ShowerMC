#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
using namespace std;
struct xyz {
    double a;
    float  b;
    int    x;
    char   z[5];
} ;
//
// subprog to read a binary file containing the struct xyz
//  

void subout(fstream *);
main(){

  fstream  fid;
  struct xyz  instxyz;
  instxyz.a=1.0;
  instxyz.b =13.0;
  instxyz.x = 5;
  strcpy(instxyz.z, "abcxy");


  //  fid.open("temp.data", ios::binary | ios::out | ios::app);
  //  out is not enough for "rewind"-equivalent operation 
  fid.open("temp.data", ios::binary | ios::out | ios::in );
  if( fid.fail() ) { 
    cerr<< " error "<< endl;  
    exit(1);
  };
  
  int nega=-1;
  fid.write( (char *)(&nega) , sizeof(nega) );
  fid.write( (char *)(&instxyz), sizeof( instxyz ));
  //  fid.close();
  //  fid.open("temp.data", ios::binary | ios::in);
  // instead of above two we can use rewind-equivalent operation
  fid.flush();
  fid.seekg(0);  // this is short for the next
  //  fid.seekg(0, ios::beg);
  //

  subout(&fid);
  /*
  while (1) {
     // reading like below, i.e, each element is read,
     // succeeds for the first time, but
     // 2nd time  on, data read is not correct.
     //	actually; first read is not strictly good;
     //	last character is strange.
     // However, if last 'char' is omitted, everything is ok.
    
    double a;
    float  b;
    int    x;
    char   z[5];


    fid.read( (char *) (&a), sizeof(a));
    if( fid.fail())  exit (0);
    fid.read( (char *) (&b), sizeof(b));  
    fid.read( (char *) (&x), sizeof(x));
    fid.read( z, sizeof(z) );

    cout << a << " " << b << " " << x <<" " << z << endl;

    struct xyz abc;
    fid.read((char *)(&abc), sizeof( abc ) );
    if( fid.fail())  exit (0);
    cout << abc.a << " " << abc.b << " " << abc.x << abc.z << endl;
  }
  */
}

void subout(fstream *pfid){
  struct xyz abc;
  int nega;
  (*pfid).read((char *)(&nega), sizeof( nega ) );
  (*pfid).read((char *)(&abc), sizeof( abc ) );
  if( (*pfid).fail())  exit (0);
  cout << nega << endl;
  cout << abc.a << " " << abc.b << " " << abc.x << abc.z << endl;
}  



