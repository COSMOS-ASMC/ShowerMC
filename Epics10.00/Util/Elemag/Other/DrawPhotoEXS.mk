include  $(EPICSTOP)/site.config



objs =  DrawPhotoEXS.o 


drawphotoe.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 


clean:;		@rm -f $(OBJS) core *~ a.out



