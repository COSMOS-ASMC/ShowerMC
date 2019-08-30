include $(COSMOSTOP)/site.config
#
#
#vpath %.h  $(COSMOSINC)

#vpath c%.f  $(COSMOSTOP)/Tracking/Util:$(COSMOSTOP)/Tracking:$(COSMOSTOP)/Tracking/Geomag:$(COSMOSTOP)/Tracking/AS



src = testInitObs.f cinitObs.f 


objs := $(src:.f=.o)


initObs :  $(objs)  
	$(LD) $(LDFLAGS) -o $@  $(objs) -L$(LIBLOC) -l$(LIBNAME)

###


cinitObs.o: $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/BlockData/cblkObs.h

