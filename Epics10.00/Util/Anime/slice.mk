include $(EPICSTOP)/site.config

objs = timeslice.o 

slice$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)   -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos





