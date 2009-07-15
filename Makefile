CXX=g++ -g -Wall -I/usr/local/include/colib $(OPT)
CC=$(CXX)
OPT=-O4 # -DUNSAFE
LDLIBS=-lm 

all: rast-test rast cedges

rast: rast.o librast.a
rast-test: rast-test.o librast.a
cedges: cedges.cc
	$(CXX) -o cedges -DMAIN cedges.cc -DUNSAFE -O4

LIBRAST=cedges.o calignmentp2d.o cinstancep2d.o \
	clinesp2d.o cliness2d.o crastp2d.o crastss2d.o crasts2d.o \
	crastrs2d.o
librast.a: $(LIBRAST)
	ar cr $@ $^

clean:
	rm -f *.o *.a

# $(LIBRAST): h/rast.h h/misc.h h/struct.h h/geo.h h/trie.h h/vecmat.h




