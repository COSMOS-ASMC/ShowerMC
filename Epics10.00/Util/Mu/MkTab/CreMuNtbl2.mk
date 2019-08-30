include  $(EPICSTOP)/site.config
vpath %.f  $(EPICSTOP)/Util/Elemag

#objs =   epmuBPNgeneI.o epmuAuxNBB.o epmuSetCnst.o epmuNdsdvKK.o testCreN.o epCreMuNTab.o epwtSmpTbl.o epresetmuVmin.o
objs =  epmuAuxNBB.o  epmuNdsdvKK.o testCreN.o epCreMuNTab.o epwtSmpTbl.o epresetmuVmin.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




