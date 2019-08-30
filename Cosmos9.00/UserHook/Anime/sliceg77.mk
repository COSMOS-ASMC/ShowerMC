include $(COSMOSTOP)/site.config

objs = timesliceg77.o 

slice$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

