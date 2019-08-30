struct xxxx {
  float z;
  int ar[5];
  int j;
  int tf; 
  char string[2][10];
} ;
struct yyyy {
  struct xxxx a1, a2;
  float a3;
};

extern struct {
  struct yyyy ppp;
} abc_;

int cleng(char *);

#include <stdio.h>  
void xyz_(float *x, float *y, int *i){
  /*  char  ps[10]; */
  char * xxx;
  int j; 
  printf(" in C  z j %f %d\n", abc_.ppp.a1.z, abc_.ppp.a2.j);
  printf(" in C   tf%d\n",  abc_.ppp.a2.tf); 
  for( j=0; j<=4; j++) printf(" ar %d %d\n", j, abc_.ppp.a2.ar[j]);
  printf(" a3 %f\n", abc_.ppp.a3);
  /*  for(j=0;j<=9; j++) ps[j]=abc_.ppp.a1.string[j];*/
  /*  ps[10] = '\n'; */
  /*  j = cleng_( &abc_.ppp.a1.j+4 );
  printf(" j = %d\n", j);
  */
  xxx = "ab\n";
  /*
  j= cleng_( xxx); 
  printf(" j = %d\n", j);
  */
  printf(" xxx %s\n", xxx);

  /*  abc_.ppp.a1.string[10] = '\n'; */
  /*  printf(" string %s\n", &ps);*/
  printf(" string %s\n", &abc_.ppp.a1.string[0][0]);
  printf(" string %s\n", &abc_.ppp.a1.string[1][0]);
  printf(" in C XYZ %f %f %d\n", *x, *y, *i);
  *x = 10.0;
  *y = 20.0;
  *i = 30.0;
  pqr_( x, i);
  printf(" in C  after pqr %f %d\n", *x, *i);
 }



      
