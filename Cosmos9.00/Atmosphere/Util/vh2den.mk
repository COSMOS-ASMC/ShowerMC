include $(COSMOSTOP)/site.config
#
objs =  vh2den.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)









