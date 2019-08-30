include $(COSMOSTOP)/site.config

objs = seeascii.o  asciiprint.o

seeascii$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)


seeascii.o:	$(COSMOSINC)/Ztrack.h \
	Zprivate.h
asciiprint.o:	$(COSMOSINC)/Ztrack.h \
	Zprivate.h
