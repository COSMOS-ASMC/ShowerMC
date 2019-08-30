include  $(EPICSTOP)/site.config


objs = epBPfuncH.o epReadMTbl.o  epGetEffZA.o epExpot.o testCreBrH.o epCreBrSTbH.o  kbchop.o epSetSTblCns.o  epX0.o epCoulombC.o epwtSmpTbl.o psiim.o epX0Old.o epCompScrBr.o epCompScrPr.o


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




