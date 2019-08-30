include $(COSMOSTOP)/site.config

objs = select.o

select$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)


select.o: $(COSMOSINC)/Ztrack.h \
	Zprivate.h


