include $(COSMOSTOP)/site.config

objs = read_show_skel.o 

read_show_skel$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)


read_show_skel.o:	$(COSMOSINC)/Ztrack.h \
	Zprivate.h

