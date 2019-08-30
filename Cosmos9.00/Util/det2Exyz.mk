include $(COSMOSTOP)/site.config

OBJS = det2Exyz.o

a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME)




