include $(COSMOSTOP)/site.config

objs = testK+pXsec.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)



