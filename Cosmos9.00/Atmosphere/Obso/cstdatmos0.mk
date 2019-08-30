
#
#
vpath %.h .
vpath copenf.f ../../../Sysdep
vpath k%.f ../../../KKlib
vpath c%.f  .
#       KKlib
src0 =  cskipComment.f  copenf.f
#
#src1 = cstdatmos0-multi-seg.f testcstdatmos0.f
src1 = cstdatmos0.f testcstdatmos0.f
obj0 := $(src0:.f=.o)
obj1 := $(src1:.f=.o)
objs = $(obj0) $(obj1)

a.out: $(objs)
	$(FC) $(LDFLAGS) -o $@ $(objs) 











