include  $(EPICSTOP)/site.config

objs =  testGetXXsec.o 

testGetXXsec$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos





