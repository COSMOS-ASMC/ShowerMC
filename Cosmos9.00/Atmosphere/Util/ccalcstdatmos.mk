include $(COSMOSTOP)/site.config
#
objs =  ccalcstdatmos.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)









