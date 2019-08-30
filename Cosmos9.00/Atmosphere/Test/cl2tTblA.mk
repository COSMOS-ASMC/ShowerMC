include $(COSMOSTOP)/site.config
#
objs =  testcl2tTblA.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)









