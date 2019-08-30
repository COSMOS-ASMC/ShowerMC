include  $(EPICSTOP)/site.config

objs =  epXCOM.o epZ2ChemSymb.o

epXCOM$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos





