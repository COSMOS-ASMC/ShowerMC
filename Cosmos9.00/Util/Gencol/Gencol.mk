include  $(COSMOSTOP)/site.config

objs =  Gencol.o getDiffCode.o getImpactParam.o

Gencol$(ARCH).out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos

clean:;		@rm -f $(OBJS) core *~ a.out

