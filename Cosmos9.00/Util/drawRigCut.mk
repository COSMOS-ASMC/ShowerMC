include  $(COSMOSTOP)/site.config


objs	      = drawRigCut.o


SRCS	      = drawRigCut.f


EXTHDRS	      = $(COSMOSINC)/Zmanagerp.h \
		$(COSMOSINC)/ZrigCut.h


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

###
drawRigCut.o: $(COSMOSINC)/Zmanagerp.h \
	$(COSMOSINC)/ZrigCut.h
