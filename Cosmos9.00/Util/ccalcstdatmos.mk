include  $(COSMOSTOP)/site.config


objs	      = ccalcstdatmos.o


SRCS	      = ccalcstdatmos.f



a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

