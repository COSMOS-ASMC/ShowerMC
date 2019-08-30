include  $(LIBLOFT)/site.config

objs =  user3.0-6.o	

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(LIBLOFT)/lib/$(ARCH) -lcosmos




