#include <stdio.h>
struct histc {
  float norm;
  float pw;
  short int eventno;
  unsigned char  logv;
  char  init[8];
  char  title[128];
  char  categ[8];
  char  id[96];
  char  dNunit[32];
  char  dir[128];
  char  dummy;
};

struct histogram1 {
  struct histc c;
  float *xw;
  float *dnw;
  float *dndx;
  float *mean;
};

void func(struct histogram1 *h);

int main(){
  struct histogram1 h, k;
  func( &h );
  int i;
  for(i=0; i<100; i++) {
    printf("h.xw=%f \n", *(h.xw+i));
  }
  printf("h.c.norm=%f\n", h.c.norm);
  printf("h.mean[0]=%f\n", *h.mean);
  printf("title is: %s\n", h.c.title);
  printf("eventno= %d\n",  h.c.eventno);
  printf("size of h.c=%d\n", sizeof(h.c));
  FILE *fp;
  fp = fopen("temp.dat", "wb");
  fwrite( &h.c.norm, sizeof(float), 1,  fp);
  fprintf(stderr, " size=%d %d\n", sizeof(h.xw), sizeof(*h.xw));
  fwrite(h.xw, sizeof(*h.xw), 100, fp);  
  fclose(fp);
  fp= fopen("temp.dat", "rb");
  float temp;
  fread(&temp, sizeof(float), 1, fp);
  fprintf(stderr, " bin dat read=%f\n", temp);
  size_t rsize;
  k.xw=(float *) malloc(100*sizeof(float));

  rsize=fread( k.xw, sizeof(*h.xw), 100, fp);
  fprintf(stderr, " rsize=%d\n", rsize);
  for(i=0; i<100; i++) {
    fprintf(stderr, "%d %f\n",i, *(k.xw+i));
  }
  k.xw = h.xw;
  fprintf(stderr, "copy: 4 %f\n", *(k.xw+3));
}
void  func(struct histogram1 *h){
  h->xw = (float *) malloc (100*sizeof (float));
  h->dnw = (float *) malloc (100*sizeof (float));
  h->dndx = (float *) malloc (100*sizeof (float));
  h->mean = (float *) malloc (100*sizeof (float));
  fprintf(stderr, "in func size of h->c=%d\n", sizeof(h->c));
  float *xw;
  int i;
  xw= h->xw;
  for(i=0; i<100; i++) {
    *xw=i;
    xw++;
  }
  float *mean;
  mean=h->mean;
  *mean=-5000.;
  char dummy[10];

  strcpy(h->c.title,"this is title");
  strcat(h->c.title, dummy);
  strcat(h->c.title," added script"); 

  h->c.norm=-1.0;

  int temp =100;
  h->c.eventno= (temp/2*2 -temp) !=0;
  printf(" event 100->%d\n", h->c.eventno);
  temp =101;
  h->c.eventno= (temp/2*2 -temp) !=0;
  return;
}



