include $(COSMOSTOP)/site.config
#

OBJS = getST2.o cprimFlux0.o intePrim2.o

a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME) -L$(LIBLOFT)/lib/$(ARCH) -lloft




