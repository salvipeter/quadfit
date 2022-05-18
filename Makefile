all: quadfit-test ribbon-test connect-test

EIGEN=/usr/include/eigen3
LIBGEOM=../libgeom
JETWRAP=../jet-wrapper

CXXFLAGS=-Wall -pedantic -std=c++17 -O0 -g -I$(LIBGEOM) -I$(JETWRAP) -I$(EIGEN)
LDFLAGS=-L$(LIBGEOM)/debug -lgeom -lomp

quadfit-test: quadfit-test.o quadfit.o bspline-fit.o io.o $(JETWRAP)/jet-wrapper.o
	$(CXX) -o $@ $^ $(LDFLAGS)

ribbon-test: ribbon-test.o fit-ribbon.o io.o
	$(CXX) -o $@ $^ $(LDFLAGS)

connect-test: connect-test.o connect-g1.o io.o
	$(CXX) -o $@ $^ $(LDFLAGS)
