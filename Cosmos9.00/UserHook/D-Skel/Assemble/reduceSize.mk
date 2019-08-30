include $(COSMOSTOP)/site.config

objs = reduceSize.o

reduceSize$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
