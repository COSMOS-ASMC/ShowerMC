include  $(EPICSTOP)/site.config

objs =  convEposRes2.o

convEpos2$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos


