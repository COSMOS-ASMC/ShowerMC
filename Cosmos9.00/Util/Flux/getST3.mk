include $(COSMOSTOP)/site.config
#

OBJS = getST3.o cprimFlux0.o intePrim2.o

a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME)





