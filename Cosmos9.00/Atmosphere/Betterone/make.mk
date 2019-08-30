include $(COSMOSTOP)/site.config
#
objs =  creadAtmosD.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)









