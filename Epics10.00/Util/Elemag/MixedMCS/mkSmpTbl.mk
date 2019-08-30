include  $(EPICSTOP)/site.config

objs =  cmkSmpTbl.o  cgetMuc.o csetSmpTblconst.o  mkSmpTblMain.o  readMedia.o getConvFac.o


mksmptab$(ARCH).out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v   -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 

