#include <strstream>
#include <iostream>
#include <string.h>
using namespace std;
  main(){
    char  msg[100];
    ostrstream mout(msg, sizeof msg);
    mout << "mesage " << ends;
    cout << msg << endl;
    mout << "message 2" << ends;
    cout << msg << endl;
    //    mout.close();
    istrstream min(msg, sizeof msg);
    strcpy(msg , "123 3.5");
    int x;
    double y;
    min >> x;
    min >> y;
    cout << x << " " << y << endl;
  }

