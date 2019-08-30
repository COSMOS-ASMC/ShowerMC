#include <stdio.h>

int main(){
  FILE *fp; 
  char buf[100];
  size_t size, rsize;
  fp=fopen("temp.dat", "wb");
  size=sizeof("xyz#xxxxxxxxx");

  fprintf(stderr,"size=%d\n", size);

  fwrite("xyz#", 1, size, fp);
  fclose(fp);

  fp=fopen("temp.dat", "rb");
  size--;
  rsize=fread(buf, 1, size, fp);
  // actually 5 chars are write/read (size=rsize=5)
  // it seems \0 is counted.
  // if we add \0 in the input string.  size becomes 6;
  //  why ? 
  fprintf(stderr, "size=%d rsize=%d\n", size, rsize);
  fprintf(stderr, "%s\n",buf);
}
