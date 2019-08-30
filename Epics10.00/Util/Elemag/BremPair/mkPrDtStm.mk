include  $(EPICSTOP)/site.config


objs = epmkPairDtStm.o epwtStabInData.o  kmkDataStm2a.o  kmkDataStm1.o



a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




