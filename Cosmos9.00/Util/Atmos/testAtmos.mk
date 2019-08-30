include $(COSMOSTOP)/site.config

OBJS = testAtmos.o

a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o  $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME)  -L$(LIBLOFT)/lib/$(ARCH) -lloft

