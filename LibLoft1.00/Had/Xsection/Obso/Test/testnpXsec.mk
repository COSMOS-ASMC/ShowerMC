include $(COSMOSTOP)/site.config

objs = testnpXsec.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)



