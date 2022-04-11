all: quadfit-test

EIGEN=/usr/include/eigen3
LIBGEOM=../libgeom

CXXFLAGS=-Wall -pedantic -std=c++17 -I$(LIBGEOM) -I$(EIGEN)
LDFLAGS=-L$(LIBGEOM)/debug -lgeom

quadfit-test: quadfit-test.o quadfit.o
	$(CXX) -o $@ $^ $(LDFLAGS)
