include $(COSMOSTOP)/site.config

OBJS = det2prim.o
a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME)




