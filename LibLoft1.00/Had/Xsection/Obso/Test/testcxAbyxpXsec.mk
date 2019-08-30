include $(COSMOSTOP)/site.config

objs = testcxAbyxpXsec.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)



