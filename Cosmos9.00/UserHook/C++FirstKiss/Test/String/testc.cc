#include <iostream>
#include <string.h>
using namespace std;

extern "C" {
  extern int klena_(char *, int );
  void cputNull(char *str) {
    str[ klena_(str, strlen(str)) ] = '\0';
  }

  extern   void ddd_();
  int main(){
    ddd_();
  }

  void xxx_(char  * a, int len){
    char x[20];
    cout << " len is " << len << "\n";
    cout << " length of a=" << strlen(a) << "\n";
    cputNull(a);
    
    cout << " after nullput; length of a=" << strlen(a) << "\n";
    cout << " in xxx a is:  " << a << "\n";
    cout << " klena(a) = " << klena_( a,  strlen(a)) << "\n";
    strncpy(x, a, klena_(a, strlen(a)));
    cout <<" x=" << x << "\n";
    cout << " length of x is " << strlen(x) <<endl;
    if( a == "append") cout << "a is append" << endl;
    else  cout << a <<  "a is not append" << endl;
    cout << "result of strncmp "<< endl;
    cout <<  strncmp(a, "append", 6);
    strcpy(x,a);
    //    x[klena_(x, strlen(x))] = '\0';
    cputNull(x);
     /*
       cout <<" x=" << x << "\n";
       cout << " length of x is " << strlen(x) <<endl;
    */
     //    strncat( x, "987", strlen(x));
     strcat( x, "987");
     //    cout << " x=" << x << "\n";
     strcpy(a, x);
  }
}

  


