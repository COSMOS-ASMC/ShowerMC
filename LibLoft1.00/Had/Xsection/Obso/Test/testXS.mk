include $(COSMOSTOP)/site.config

objs = testXS.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)



