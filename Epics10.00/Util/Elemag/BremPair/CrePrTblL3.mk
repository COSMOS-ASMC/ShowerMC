include  $(EPICSTOP)/site.config


objs = epBPfunc3.o epReadMTbl.o  epGetEffZA.o epExpot.o testCrePair.o epCrePrSTbl1.o  epBPAux.o kbchop.o epSetSTblCns.o exX0Old.o epX0.o epCoulombC.o epwtSmpTbl.o


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




