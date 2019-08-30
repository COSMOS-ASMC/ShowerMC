objs = testReadAscii.o k90whistReadAscii.o k90whist1.o

testReadAscii: $(objs)
	cc -o $@  $(objs) -lm 
