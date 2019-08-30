include $(COSMOSTOP)/site.config

objs = test1.o 

test1$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

