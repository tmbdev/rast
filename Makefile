PYTHON=python2.6
PYINC=/usr/include/$(PYTHON)
CXX=g++ -g -Wall -I/usr/local/include/colib $(OPT)
CC=$(CXX)
OPT=-O4 # -DUNSAFE
LDLIBS=-lm 

all: rast-test rast cedges _rast.so

rast: rast.o librast.a
rast-test: rast-test.o librast.a
cedges: cedges.cc
	$(CXX) -o cedges -DMAIN cedges.cc -DUNSAFE -O4

LIBRAST=cedges.o calignmentp2d.o cinstancep2d.o \
	clinesp2d.o cliness2d.o crastp2d.o crastss2d.o crasts2d.o \
	crastrs2d.o
librast.a: $(LIBRAST)
	ar cr $@ $^
_rast.so: rast.i librast.a
	swig -python -c++ rast.i 
	g++ -g -fPIC -I$(PYINC) -shared rast_wrap.cxx -o _rast.so librast.a

clean:
	rm -f *.so *wrap.cxx *.o rast.py *.so


