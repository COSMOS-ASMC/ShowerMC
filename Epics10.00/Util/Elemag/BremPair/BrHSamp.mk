include  $(EPICSTOP)/site.config


objs = testBrHSamp.o  eprdSmpTbl.o epBrHSamp.o epGetEffZA.o epReadMTbl.o epX0.o epCoulombC.o


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




