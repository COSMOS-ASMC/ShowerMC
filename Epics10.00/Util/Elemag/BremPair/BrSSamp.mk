include  $(EPICSTOP)/site.config


objs = testBrSSamp.o  eprdSmpTbl.o epBrSSamp.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




