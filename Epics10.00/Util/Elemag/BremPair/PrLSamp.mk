include  $(EPICSTOP)/site.config


objs = testPrLSamp.o  eprdSmpTbl.o epPrLSamp.o 


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




