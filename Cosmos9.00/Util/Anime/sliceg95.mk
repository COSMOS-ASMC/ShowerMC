include $(COSMOSTOP)/site.config

objs = timeslice.o 

slice$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME) -L$(LIBLOFT)/lib/$(ARCH) -lloft






