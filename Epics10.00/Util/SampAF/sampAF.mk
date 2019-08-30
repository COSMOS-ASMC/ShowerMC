include $(EPICSTOP)/site.config

objs = testSampAF.o

sampAF$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v -O1 -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos

#  -I$(DEST) -I$(COSMOSTOP)/lib/$(ARCH)




# not needed
# testSampAF.o: $(COSMOSINC)/csampAF.f90
