include  $(EPICSTOP)/site.config

objs = testNewVol.o    episoDirec.o



a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 


clean:;		@rm -f $(objs) core *~ a.out temp*.f



