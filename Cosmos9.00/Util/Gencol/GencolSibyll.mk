include  $(EPICSTOP)/site.config

objs =  Gencol.o  sibyllinit.o  sibyll2.1.o getDiffCode.o getImpactParam.o

Gencol$(ARCH).out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos

clean:;		@rm -f $(OBJS) core *~ a.out

