include  $(EPICSTOP)/site.config
vpath %.f  $(EPICSTOP)/Util/Elemag

#objs =   epmuBPNgeneI.o epmuAuxP.o epmuSetCnst.o  epCreMuPTab.o epwtSmpTbl.o   epmudsdvdr.o  epmudsdv.o testCrePr.o epresetmuVmin.o epmuvmax.o
objs =  epCreMuPTab.o epwtSmpTbl.o  testCrePr.o epresetmuVmin.o


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




