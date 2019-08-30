include $(COSMOSTOP)/site.config

objs = timeslice2.o 

slice2$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

