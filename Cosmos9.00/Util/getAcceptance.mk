include  $(COSMOSTOP)/site.config


objs	      = getAcceptance.o cprimAcceptance.o 


SRCS	      = getAcceptance.f cprimAcceptance.f


EXTHDRS	      = $(COSMOSINC)/Zmanagerp.h \
		$(COSMOSINC)/ZrigCut.h


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)


