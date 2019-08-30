include $(COSMOSTOP)/site.config
#
#
#vpath %.h  $(COSMOSINC)

#vpath c%.f  $(COSMOSTOP)/Tracking



src = testMkInc.f

objs := $(src:.f=.o)


mkInc:  $(objs)  
	$(LD) $(LDFLAGS) -o $@  $(objs) -L$(LIBLOC) -l$(LIBNAME)
	touch mkInc
###


cmkIncident.o: $(LIBLOFTINC)/Ztrack.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/Zincidentp.h \
	$(LIBLOFTINC)/Zincidentv.h 



