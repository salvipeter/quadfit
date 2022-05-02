all: quadfit-test

EIGEN=/usr/include/eigen3
LIBGEOM=../libgeom
JETWRAP=../jet-wrapper

CXXFLAGS=-Wall -pedantic -std=c++17 -O0 -g -I$(LIBGEOM) -I$(JETWRAP) -I$(EIGEN)
LDFLAGS=-L$(LIBGEOM)/debug -lgeom -lomp

quadfit-test: quadfit-test.o quadfit.o bspline-fit.o $(JETWRAP)/jet-wrapper.o
	$(CXX) -o $@ $^ $(LDFLAGS)
