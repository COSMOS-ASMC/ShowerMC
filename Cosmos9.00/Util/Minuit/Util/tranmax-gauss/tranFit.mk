include $(COSMOSTOP)/site.config

objs = tranfnc.o tranFit.o

tranFit$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -lminuit 
	ln -sf tranFit$(ARCH) tranFit
