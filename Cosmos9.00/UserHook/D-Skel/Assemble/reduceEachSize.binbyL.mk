include $(COSMOSTOP)/site.config

objs = reduceEachSize.binbyL.o

reduceEachSize.binbyL.$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
