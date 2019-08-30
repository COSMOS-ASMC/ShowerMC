include $(COSMOSTOP)/site.config
#
objs =  testcl2tTbl.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)









