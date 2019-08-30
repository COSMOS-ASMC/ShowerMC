include $(COSMOSTOP)/site.config
#
#
#vpath %.h  $(COSMOSINC)
# vpath k%.f $(COSMOSTOP)/KKlib
#vpath c%.f  $(COSMOSTOP)/Manager:$(COSMOSTOP)/Primary
#vpath testPrim.f $(COSMOSTOP)/Primary/Util

#  OBJS	      = cerrorMsg.o \
		copenf.o    \
		cskipComment.o \
	        testPrim.o  \
		cmkPrimSTbl.o \
		cprintPrim.o \
		cprocPrimDt.o \
		crdPrimData.o \
		ciniSPrim.o \
		csampPrimary.o
 OBJS = testPrim.o

a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME)

