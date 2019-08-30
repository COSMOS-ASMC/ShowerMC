include  $(EPICSTOP)/site.config



SRCS = testMuBrSamp.f epMuBrEcheck.f
objs = testMuBrSamp.o epMuBrEcheck.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




