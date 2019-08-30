include $(EPICSTOP)/site.config

objs = timeslice.o 

slice$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v   -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos


###	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME)




