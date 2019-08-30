include $(COSMOSTOP)/site.config

objs = bin2bin.o

bin2bin$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
