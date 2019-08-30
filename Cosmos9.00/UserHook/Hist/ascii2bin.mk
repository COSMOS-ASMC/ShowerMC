include $(COSMOSTOP)/site.config

objs = ascii2bin.o 


ascii2bin$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
