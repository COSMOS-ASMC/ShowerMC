include $(COSMOSTOP)/site.config
#
#


 OBJS = atmosd.o

a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME)

