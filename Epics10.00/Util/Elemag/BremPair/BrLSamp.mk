include  $(EPICSTOP)/site.config


objs = testBrLSamp.o  eprdSmpTbl.o epBrLSamp.o


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




