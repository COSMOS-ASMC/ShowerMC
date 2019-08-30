include  $(LIBLOFT)/site.config

objs =  testBrem.o

testBrem(ARCH): $(objs)
	$(LD) $(LDFLAGS)  -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(LIBLOFT)/lib/$(ARCH) -lloft







