include  $(COSMOSTOP)/site.config


objs	      = drawRigCut2.o


SRCS	      = drawRigCut2.f


EXTHDRS	      = $(COSMOSINC)/Zmanagerp.h \
		$(COSMOSINC)/ZrigCut.h


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

