include $(COSMOSTOP)/site.config

objs = mknrfai.o

mknrfai$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
