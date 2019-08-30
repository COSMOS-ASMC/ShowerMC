include $(COSMOSTOP)/site.config
#

OBJS = getST.o cprimFlux.o cprimAcceptance.o

a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME)





