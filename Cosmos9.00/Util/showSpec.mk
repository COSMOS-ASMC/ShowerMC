include  $(COSMOSTOP)/site.config


objs	      = cprimFlux.o cprimAcceptance.o  showSpec.o cprimFlux0.o 


SRCS	      = cprimFlux.f  showSpec.f cprimAcceptance.f cprimFlux0.f  


EXTHDRS	      = $(COSMOSINC)/Zmanagerp.h \
		$(COSMOSINC)/ZrigCut.h


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

