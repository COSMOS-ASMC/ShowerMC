include $(COSMOSTOP)/site.config

objs = add2hist.o 

add2hist$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

