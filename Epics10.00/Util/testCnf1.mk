include  $(EPICSTOP)/site.config

objs = testCnf1.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos  

clean:;		@rm -f $(objs) core *~ a.out temp*.f
