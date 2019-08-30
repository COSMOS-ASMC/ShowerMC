include  $(EPICSTOP)/site.config



objs = epBPfunc3.o  epBPAux.o  epwtSmpTbl.o DrawBrems2.o epBrSeltzer.o epBPfuncH.o  epBPAuxH.o epPairLowE.o epBrgeneric.o  epNormLPM.o epBPgeneini.o 



drawbrems2_3.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




