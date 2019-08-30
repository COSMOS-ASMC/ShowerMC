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
  struct histogram1 h;
  func( &h );
  int i;
  for(i=0; i<100; i++) {
    printf("h.xw=%f \n", *(h.xw+i));
  }
  printf("h.c.norm=%f\n", h.c.norm);
  printf("h.mean[0]=%f\n", *h.mean);
  printf("title is: %s\n", h.c.title);
  printf("eventno= %d\n",  h.c.eventno);
}
void  func(struct histogram1 *h){
  h->xw = (float *) malloc (100*sizeof (float));
  h->dnw = (float *) malloc (100*sizeof (float));
  h->dndx = (float *) malloc (100*sizeof (float));
  h->mean = (float *) malloc (100*sizeof (float));
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



