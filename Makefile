all: quadfit-test

LIBGEOM=../libgeom

CXXFLAGS=-Wall -pedantic -std=c++17 -I$(LIBGEOM)
LDFLAGS=-L$(LIBGEOM)/debug -lgeom

quadfit-test: quadfit-test.o quadfit.o
	$(CXX) -o $@ $^ $(LDFLAGS)
