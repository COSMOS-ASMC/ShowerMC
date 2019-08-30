include $(COSMOSTOP)/site.config

objs = bin2ascii.o

bin2ascii$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
