#include <iostream>
#include <strstream>
#include <string.h>
#include <fstream>
extern "C" {
extern void ksplit_(char *, int *, int *, char *, int *, int, int);
extern void cerrormsg_(char *, int *, int);
}
using namespace std;
int main(){
  char  msg[100];
  char  str[20][3];
  char  str20[20];
  char  stro[5][10]= { "123456789",
		       "234567890", "345678901", 
		       "456789012","567890123"};
  ostrstream mout(msg, sizeof msg);
  int x=100;
  mout << " This is my name " <<  x << ends;
  int temp = 1;
  cerrormsg_(msg,  &temp, strlen(msg) );
  x = -1;
  mout.seekp(0);
  mout << "Additional " << x  << ends;
  cout << " bef cerrormsg "<< endl;
  cerrormsg_(msg,  &temp, strlen(msg) );
  cout << " aft cerrormsg "<< endl;

  strcpy(str20, "abc de fgx");
  //  cout << str20 << endl;
  int temp1 = 10;
  int temp2 = 5;
  int nr;

  ksplit_(str20,  &temp1, &temp2, stro[0], &nr, 20, 10);
  cout << " nr=" << nr << endl;
  strncpy(str20, &stro[0][0],10);
  cout  << str20 << endl;

  strncpy(str20, &stro[1][0],10);
  cout << str20 << endl;
  cout  << stro[1] << endl;
  strncpy(str20, &stro[2][0],10);
  cout << str20 << endl;


  cin >>str20;


  cout << str20 << endl;

  ifstream fin("temp.data");
  if( fin.fail() ) {
    cerr<< " error "<< endl;
  }
  else {
    cerr << " ok" << endl;
    fin >>  x;
    cout << x << " x value " << endl;
  }
}


