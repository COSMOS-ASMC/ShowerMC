include  $(EPICSTOP)/site.config


objs = testUrban.o  


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos  -J$(COSMOSTOP)/lib/$(ARCH)



clean:;		@rm -f $(OBJS) core *~ a.out




