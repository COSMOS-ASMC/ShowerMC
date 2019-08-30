include  $(EPICSTOP)/site.config

objs = drawConfig.o epdrawComp.o mediaCount.o  epDrawNewVol.o

drawconfig$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos


clean:;		@rm -f $(OBJS) core *~ a.out

