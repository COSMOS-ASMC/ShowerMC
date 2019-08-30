include $(COSMOSTOP)/site.config

objs = ceposXsecB.o

ceposXsecB$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)
