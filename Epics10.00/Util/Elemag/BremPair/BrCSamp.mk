include  $(EPICSTOP)/site.config


objs = testBrCSamp.o epGetEffZA.o epReadMTbl.o epX0.o epCoulombC.o epBrCSamp.o epCompScrBr.o epSetSTblCns.o


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




