#ifndef ZsimHeader_
#define ZsimHeader_
struct simHeader {
   int event;
   int geom;
   int good;
   int code;
   int charge;
   int kerg;
   double wx, wy, wz;
   double x0, y0;
   double xbgot, ybgot;
   double xbgob, ybgob;
   double firstz;
   double sumEAnti, sumESi, sumESciFi, sumEBGOTop, sumEBGO, sumE;
};
#endif

