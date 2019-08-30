include $(COSMOSTOP)/site.config

objs = timefnc.o timeFit.o

timeFit$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -lminuit 
	ln -sf timeFit$(ARCH) timeFit
