include $(COSMOSTOP)/site.config

objs = ranseeascii.o  asciiprint.o

ranseeascii$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)


ranseeascii.o:	$(COSMOSINC)/Ztrack.h \
	Zprivate.h
asciiprint.o:	$(COSMOSINC)/Ztrack.h \
	Zprivate.h
