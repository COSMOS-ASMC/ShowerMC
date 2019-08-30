#include <stdio.h>
const double  i=1000;
const double  j=1000e3;
/* next line cannot be compiled */
const int  k  = (int) (j/i);
int main(){
  char  zz[k];
  zz[0] = 'a';
  zz[1]= '\0';
  printf("%s\n", zz);
}

