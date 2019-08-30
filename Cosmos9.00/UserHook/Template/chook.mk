include $(COSMOSTOP)/site.config

objs = chook.o

cosmos$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)  -L$(LIBLOFT)/lib/$(ARCH) -lloft


chook.o: $(LIBLOFTINC)/Zmanagerp.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Ztrackp.h
