include  $(EPICSTOP)/site.config


objs = testPrCSamp.o epGetEffZA.o epReadMTbl.o epX0.o epCoulombC.o epPrCSamp.o epCompScrPr.o epSetSTblCns.o epCompScrBr.o 


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




