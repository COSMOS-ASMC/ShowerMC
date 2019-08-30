include $(COSMOSTOP)/site.config
#
#


src = testSPrimAng.f csPrimAng.f

objs := $(src:.f=.o)


sPrimAng :  $(objs)  
	$(LD) $(LDFLAGS) -o $@  $(objs) -L$(LIBLOC) -l$(LIBNAME)
	touch sPrimAng
###

testSPrimAng.o: $(LIBLOFTINC)/Ztrack.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zincidentp.h \
	$(LIBLOFTINC)/Zincidentv.h \
	$(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/BlockData/cblkElemag.h \
	$(LIBLOFTINC)/Zelemagp.h

$(LIBRARY)(csPrimAng.o): $(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zincidentp.h \
	$(LIBLOFTINC)/Ztrack.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zincidentv.h \
	$(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/BlockData/cblkIncident.h


