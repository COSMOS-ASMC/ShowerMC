include $(COSMOSTOP)/site.config

objs = chookSkel.o

skel$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME) $(LDFLAGS)


chookSkel.o: $(COSMOSINC)/Zmanagerp.h \
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


