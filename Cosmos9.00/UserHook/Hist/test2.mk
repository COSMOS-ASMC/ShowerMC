include $(COSMOSTOP)/site.config

objs = test2.o

test2$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

