include $(COSMOSTOP)/site.config

objs = reducebinSize.o

reducebinSize$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
