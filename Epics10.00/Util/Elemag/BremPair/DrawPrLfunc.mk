include  $(EPICSTOP)/site.config

SRCS = DrawPrLfunc.f epPairLowE.f
objs = DrawPrLfunc.o epPairLowE.o




a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




