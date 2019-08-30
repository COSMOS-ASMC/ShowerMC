include $(COSMOSTOP)/site.config
#
objs =  testcstdatmos0.o

a.out: $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)  -L$(LIBLOFT)/lib/$(ARCH) -lloft


