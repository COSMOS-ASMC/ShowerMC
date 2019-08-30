include  $(EPICSTOP)/site.config
vpath %.f  $(EPICSTOP)/Util/Elemag


#objs =   epmuBPNgeneI.o epmuAuxB.o epmuSetCnst.o epmuBrem.o testCreBr.o epCreMuBTab.o epwtSmpTbl.o  epresetmuVmin.o epmuvmax.o

objs = epmuAuxB.o  epmuBrem.o testCreBr.o epCreMuBTab.o  epresetmuVmin.o epwtSmpTbl.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




