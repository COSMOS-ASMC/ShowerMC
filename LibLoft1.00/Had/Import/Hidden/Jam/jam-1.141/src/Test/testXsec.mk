include $(COSMOSTOP)/site.config

objs = testXsec.o

testXsec$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME) $(LDFLAGS)

