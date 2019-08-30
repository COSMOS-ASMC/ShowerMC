include $(COSMOSTOP)/site.config

objs = reanal.o

reanal$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)


reanal.o:	Zprivate.h \
	$(COSMOSINC)/Ztrack.h



