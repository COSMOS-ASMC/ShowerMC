include $(COSMOSTOP)/site.config

objs = add2hyb.o 

add2hyb$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

