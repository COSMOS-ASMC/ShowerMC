include $(COSMOSTOP)/site.config
# vpath %.h  $(COSMOSINC)
# vpath Geomag.f:$(COSMOSTOP)/Util
objs = Geomag.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L/$(LIBLOC) -l$(LIBNAME)


Geomag.o: $(COSMOSINC)/Zcoord.h \
	$(COSMOSINC)/Zmagfield.h
