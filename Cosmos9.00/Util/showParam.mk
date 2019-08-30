include $(COSMOSTOP)/site.config
#
#
#vpath %.h  $(COSMOSINC)
#vpath showParam.f:$(COSMOSTOP)/Util
OBJS	      = showParam.o	


showParam: $(OBJS)
	$(LD) $(LDFLAGS) -o $@  $(OBJS) -L/$(LIBLOC) -l$(LIBNAME)


showParam.o: $(COSMOSINC)/Zmanagerp.h \
	$(COSMOSINC)/BlockData/cblkEvhnp.h \
	$(COSMOSINC)/Zevhnp.h \
	$(COSMOSINC)/BlockData/cblkHeavy.h \
	$(COSMOSINC)/Zcode.h \
	$(COSMOSINC)/Zheavyp.h \
	$(COSMOSINC)/BlockData/cblkIncident.h \
	$(COSMOSINC)/Zincidentp.h \
	$(COSMOSINC)/BlockData/cblkElemag.h \
	$(COSMOSINC)/Zelemagp.h \
	$(COSMOSINC)/BlockData/cblkObs.h \
	$(COSMOSINC)/Zcoord.h \
	$(COSMOSINC)/Zobs.h \
	$(COSMOSINC)/Zobsp.h \
	$(COSMOSINC)/BlockData/cblkTracking.h \
	$(COSMOSINC)/Ztrackp.h \
	$(COSMOSINC)/BlockData/cblkXsec.h \
	$(COSMOSINC)/Zxsectionp.h
