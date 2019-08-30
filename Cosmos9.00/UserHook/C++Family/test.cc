
#include <iostream.h>
#include "Zprivate.h"
void cccc(double x[][5], int, int);

int main(){
  double x[2][5] = {{1.0, 2.0, 3.0, 4.0, 5.0},
  {-2, -3, -1, -4, -5}};
  cccc( x, 2, 5);
  com.nnn = 0;
}
void cccc(double x[][5], int p, int q) {
  for(int i=0; i< q; i++) {
    for(int j=0; j< p; j++) {
      cout << x[j][i];
    }
    cout << endl;
  }
}

