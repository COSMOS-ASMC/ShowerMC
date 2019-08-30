include $(COSMOSTOP)/site.config
#
#
#vpath %.h  $(COSMOSINC)
#vpath k%.f $(COSMOSTOP)/KKlib
#vpath c%.f  $(COSMOSTOP)/Manager:$(COSMOSTOP)/Primary
# vpath testSampPrim.f $(COSMOSTOP)/Primary/Util

# OBJS	      = cerrorMsg.o \
		testSampPrim.o \
		cmkPrimSTbl.o \
		cprintPrim.o \
		cprocPrimDt.o \
		crdPrimData.o \
		ciniSPrim.o \
		csampPrimary.o
OBJS = testSampPrim.o

a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME)

testSampPrim.o: $(COSMOSINC)/Zptcl.h \
        $(COSMOSINC)/Zprimary.h 




