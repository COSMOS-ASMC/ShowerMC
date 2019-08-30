include $(COSMOSTOP)/site.config

objs = randomsel.o

randomsel$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)


randomsel.o: $(COSMOSINC)/Ztrack.h \
	Zprivate.h


