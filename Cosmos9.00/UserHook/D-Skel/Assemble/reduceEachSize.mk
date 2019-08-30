include $(COSMOSTOP)/site.config

objs = reduceEachSize.o

reduceEachSize$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
