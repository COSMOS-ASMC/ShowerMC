include $(LIBLOFT)/site.config

objs = test.o

test$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)


test.o: $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Ztrack.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/Ztrackp.h
