include  $(EPICSTOP)/site.config

objs =  readcpy.o 

readcpy$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos


