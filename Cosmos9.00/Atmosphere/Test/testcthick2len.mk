include $(COSMOSTOP)/site.config

objs = testcthick2len.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
