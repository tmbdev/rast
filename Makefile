PYTHON=python3
PYINC=/usr/include/$(PYTHON)
CXX=g++ -g -Wall $(OPT)
CC=$(CXX)
OPT=-O3 -fPIC
LDLIBS=-lm

VPATH = src:tests:bindings/python
SRCDIR = src

all: rast-test rast cedges rast.so tests

rast: rast.o librast.a
rast-test: rast-test.o librast.a
cedges: cedges.cc
	$(CXX) -Wno-error -o cedges -DMAIN $< -DUNSAFE -O3

LIBRAST=cedges.o calignmentp2d.o cinstancep2d.o \
	clinesp2d.o cliness2d.o crastp2d.o crastss2d.o crasts2d.o \
	crastrs2d.o
librast.a: $(LIBRAST)
	ar cr $@ $^

# cedges.cc is legacy code that doesn't conform to the strict warning set
# (heavy float<->double mixing, internal shadowing). Compile it without
# -Werror so the rest of the strict checks still apply project-wide.
cedges.o: cedges.cc
	$(CXX) -Wno-error -c $< -o $@

# Python bindings via pybind11. The pybind11 headers are in the conda
# environment's include path; PYINC is the matching Python.h dir.
rast.so: rast_pybind.cc librast.a
	$(CXX) -fPIC -I$(SRCDIR) -I$(PYINC) -shared $< librast.a -o $@

DOCTESTS = calignmentp2d_test.o cinstancep2d_test.o clinesp2d_test.o cliness2d_test.o \
	crastp2d_test.o crastrs2d_test.o crasts2d_test.o crastss2d_test.o
tests: test_main.o $(DOCTESTS) librast.a
	$(CXX) -o tests test_main.o $(DOCTESTS) librast.a -lm

install:
	cp rast.so /usr/local/lib/$(PYTHON)/dist-packages/.
	chmod ugo+rX /usr/local/lib/$(PYTHON)/dist-packages/rast.so
	cp librast.a /usr/local/lib
	chmod ugo+rX /usr/local/lib/librast.a
	cp rast rast-test cedges /usr/local/bin
	chmod ugo+rX /usr/local/bin/rast
	chmod ugo+rX /usr/local/bin/cedges
	chmod ugo+rX /usr/local/bin/rast-test
	cp $(SRCDIR)/rast.h /usr/local/include
	chmod ugo+rX /usr/local/include/rast.h

clean:
	rm -f *.so *.o *.a rast rast-test cedges tests
