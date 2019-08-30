include $(EPICSTOP)/site.config

objs = testSampAF2.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v -O1 -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos

