include $(COSMOSTOP)/site.config

objs = chook.o

cosmos$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)


chook.o: $(COSMOSINC)/Zmanagerp.h \
	$(COSMOSINC)/Ztrack.h \
	$(COSMOSINC)/Zptcl.h \
	$(COSMOSINC)/Zcoord.h \
	$(COSMOSINC)/Zpos.h \
	$(COSMOSINC)/Zdirec.h \
	$(COSMOSINC)/Zmagfield.h \
	$(COSMOSINC)/Zobs.h \
	$(COSMOSINC)/Zobsp.h \
	$(COSMOSINC)/Zobsv.h \
	$(COSMOSINC)/Ztrackv.h \
	$(COSMOSINC)/Ztrackp.h
