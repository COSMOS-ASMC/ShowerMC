include $(COSMOSTOP)/site.config

objs = testReadAscii.o 

testReadAscii$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
