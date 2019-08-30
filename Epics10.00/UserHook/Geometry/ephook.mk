include  $(EPICSTOP)/site.config

objs =   epGeom.o ephook.o eppos2B.o

sepics$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos


ephook.o: $(EPICSINC)/ZepTrack.h \
	$(EPICSINC)/Zmove.h \
	$(EPICSINC)/BlockData/epblkepi.h \
	$(EPICSINC)/BlockData/epblksepi.h \
	$(EPICSINC)/BlockData/epblkCnfig.h 



