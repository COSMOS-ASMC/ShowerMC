include $(COSMOSTOP)/site.config

objs = testgeomagMap.o 

testgeomagMap$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
