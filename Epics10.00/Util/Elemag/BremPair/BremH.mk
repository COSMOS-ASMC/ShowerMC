include  $(EPICSTOP)/site.config


objs = testBremH.o  epReadMTbl.o  psiim.o epBPfuncH.o epGetEffZA.o epX0.o epCoulombC.o epX0Old.o epCompScrBr.o epCompScrPr.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




