include $(COSMOSTOP)/site.config

objs = testAAXs.o 
AAxs$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME) $(LDFLAGS)
