include $(COSMOSTOP)/site.config

objs = add2nrfai.o 

add2nrfai$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

