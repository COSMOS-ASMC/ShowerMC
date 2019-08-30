include  $(EPICSTOP)/site.config

objs = epBrSeltzer.o testCreBrS.o epwtSmpTbl.o eprdSmpTbl.o epReadMTbl.o epSetSTblCns.o  kbchop.o epGetEffZA.o epX0.o epCoulombC.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




