include $(COSMOSTOP)/site.config

objs = csoftenPiK.o


test$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME) $(LDFLAGS) 
