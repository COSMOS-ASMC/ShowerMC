include $(COSMOSTOP)/site.config

objs = testpbarpXsec.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)



