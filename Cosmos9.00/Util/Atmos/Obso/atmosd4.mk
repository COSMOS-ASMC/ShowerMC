include $(COSMOSTOP)/site.config
#
#


 OBJS = atmosd4.o

a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME)

