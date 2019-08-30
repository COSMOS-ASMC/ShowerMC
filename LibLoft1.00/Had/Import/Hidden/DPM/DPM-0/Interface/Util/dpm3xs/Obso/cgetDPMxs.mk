include $(COSMOSTOP)/site.config

objs = cgetDPMxs.o

a.out: $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)

