include $(COSMOSTOP)/site.config

objs = reduceEachSize.bin.o

reduceEachSize.bin.$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
