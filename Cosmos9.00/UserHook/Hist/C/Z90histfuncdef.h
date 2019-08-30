#ifndef Z90histodef_
#define Z90histodef_

#include "Z90histc.h"
#include "Z90histo.h"
#include "Z90hist1.h"
void kwhisti(struct histogram1 *h, float ixmin, float ibinORxmax,int inbin, int itklg);

void  kwhistc(struct histogram1 *h);

void kwhist(struct histogram1 *h, float x, float w);
void kwhists0(int fromwhich);
void kwhists(struct histogram1 *h, float inorm);
void kwhistev(struct histogram1 *h, int evno);
void kwhistid(struct histogram1 *h, char *id);
void kwhistai(struct  histogram1 *h, char *title, char *categ, char * dNunit,
	      int logv, float pw, char * label, char *unit);
void kwhistdir(struct histogram1 *h, char *dir) ;
void  kwhistpr(struct histogram1 *h, FILE *fno) ;
void kwhistw(struct histogram1  *h, FILE *bfnow);
void kwhistr(struct histogram1 *h, FILE *bfnor, int icon);
void kwhistd(struct histogram1 *h);
void kwhista(struct histogram1 *h1, struct histogram1 *h2, struct histogram1 *h);
void kwhistso(int binw);
void kwhistp(struct histogram1 *h, FILE *fno);
int  k90whistReadAscii(struct histogram1 *h, FILE *fno);
int  kwhistIxy(struct histogram1 *h, double x[], double y[], int n);
int  kwhistxy(struct histogram1 *h, double x[], double y[], int n);
#endif
