PYTHON=python3
PYINC=/usr/include/$(PYTHON)
CXX=g++ -g -Wall $(OPT)
CC=$(CXX)
OPT=-O3 -fPIC # -DUNSAFE
LDLIBS=-lm

VPATH = src:tests:bindings/python

all: rast-test rast cedges _rast.so

rast: rast.o librast.a
rast-test: rast-test.o librast.a
cedges: cedges.cc
	$(CXX) -o cedges -DMAIN $< -DUNSAFE -O3

LIBRAST=cedges.o calignmentp2d.o cinstancep2d.o \
	clinesp2d.o cliness2d.o crastp2d.o crastss2d.o crasts2d.o \
	crastrs2d.o
librast.a: $(LIBRAST)
	ar cr $@ $^

_rast.so: rast.i librast.a
	swig -python -c++ -Isrc -outdir . -o rast_wrap.cxx $<
	$(CXX) -fPIC -I$(PYINC) -shared rast_wrap.cxx -o _rast.so librast.a

install:
	cp _rast.so rast.py /usr/local/lib/$(PYTHON)/dist-packages/.
	chmod ugo+rX /usr/local/lib/$(PYTHON)/dist-packages/_rast.so
	chmod ugo+rX /usr/local/lib/$(PYTHON)/dist-packages/rast.py
	cp librast.a /usr/local/lib
	chmod ugo+rX /usr/local/lib/librast.a
	cp rast rast-test cedges /usr/local/bin
	chmod ugo+rX /usr/local/bin/rast
	chmod ugo+rX /usr/local/bin/cedges
	chmod ugo+rX /usr/local/bin/rast-test
	cp src/rast.h /usr/local/include
	chmod ugo+rX /usr/local/include/rast.h

clean:
	rm -f *.so *wrap.cxx *.o rast.py *.a rast rast-test cedges
