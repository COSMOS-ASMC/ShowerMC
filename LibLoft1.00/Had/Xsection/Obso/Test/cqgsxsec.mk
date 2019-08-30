include $(COSMOSTOP)/site.config

objs = testQGSxsec.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)



