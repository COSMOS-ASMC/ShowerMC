include $(COSMOSTOP)/site.config


objs = procLat0.o procLat.o 

procLat$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
