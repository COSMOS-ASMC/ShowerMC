include  $(EPICSTOP)/site.config


objs = epBPfuncH.o epReadMTbl.o  epGetEffZA.o epExpot.o testCrePrH.o epCrePrSTblH.o  kbchop.o epSetSTblCns.o epX0Old.o epX0.o epCoulombC.o epwtSmpTbl.o epCompScrPr.o epCompScrBr.o psiim.o


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




