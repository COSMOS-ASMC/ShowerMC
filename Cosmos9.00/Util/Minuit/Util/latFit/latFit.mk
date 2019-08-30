include $(COSMOSTOP)/site.config

objs = latfnc.o latFit.o

latFit$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -lminuit 
	ln -sf latFit$(ARCH) latFit
