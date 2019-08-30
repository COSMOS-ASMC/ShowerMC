include $(COSMOSTOP)/site.config

objs = getBasicHistoInfo.o

getBasicHistoInfo${ARCH}: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)
