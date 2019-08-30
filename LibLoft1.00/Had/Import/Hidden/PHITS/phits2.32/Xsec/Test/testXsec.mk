include $(COSMOSTOP)/site.config

objs = testXsec.o

xs$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)
