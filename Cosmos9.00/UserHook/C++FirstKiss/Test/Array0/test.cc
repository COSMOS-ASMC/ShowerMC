extern "C" {
  extern void test1_(int * );
}
extern "C" {
  extern void test2_();
}


extern "C" 
 {
  extern struct {

    float  abc[11];
    int    xx[4];
    short int  yy[2];

  } fff_;
 }
void main(){
  int cc[3];
  int i;
  for(i=0; i<=2; i++) cc[i] = -i;
  test1_(cc);
  for(i=0; i<=2; i++) {
    fff_.abc[i] =i;
    fff_.xx[i] = i;
    fff_.yy[i] = i;
  };
  test2_();
}








