include $(COSMOSTOP)/site.config

objs = NRLGenData.o
	
NRLGenData$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)

