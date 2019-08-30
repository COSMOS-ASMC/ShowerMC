include  $(EPICSTOP)/site.config

objs = tree2index.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos


clean:;		@rm -f $(OBJS) core *~ a.out

