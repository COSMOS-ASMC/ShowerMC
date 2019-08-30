include $(COSMOSTOP)/site.config

objs = test3.o 

test3$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
