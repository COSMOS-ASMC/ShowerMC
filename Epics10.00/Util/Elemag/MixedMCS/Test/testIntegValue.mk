include $(COSMOSTOP)/site.config

objs =   testIntegValue.o  cixsec2.o 

Inteval$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)

