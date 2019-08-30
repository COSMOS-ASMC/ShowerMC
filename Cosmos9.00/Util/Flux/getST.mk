include $(COSMOSTOP)/site.config
#

OBJS = getST.o cprimFlux0.o intePrim2.o

a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME)





