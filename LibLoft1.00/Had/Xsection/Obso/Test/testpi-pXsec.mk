include $(COSMOSTOP)/site.config

objs = testpi-pXsec.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)



