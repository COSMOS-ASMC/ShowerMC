include  $(EPICSTOP)/site.config



objs =  drawCompFunc.o 


drawCompFunc.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 


clean:;		@rm -f $(OBJS) core *~ a.out


