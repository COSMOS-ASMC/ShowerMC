include  $(EPICSTOP)/site.config
vpath %.f $(EPICSTOP)/prog

objs = mkNewVolume.o  


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(objs) core *~ a.out temp*.f




