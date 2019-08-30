include $(COSMOSTOP)/site.config
#
#


 OBJS = atmosd2.o

a.out: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L$(LIBLOC) -l$(LIBNAME)

