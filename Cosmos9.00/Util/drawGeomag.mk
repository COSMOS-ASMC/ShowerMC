include $(COSMOSTOP)/site.config

objs = drawGeomag.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L/$(LIBLOC) -l$(LIBNAME)


a.out: $(COSMOSINC)/Zcoord.h \
	$(COSMOSINC)/Zmagfield.h
