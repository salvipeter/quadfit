all: quadfit-test ribbon-test connect-test

EIGEN=/usr/include/eigen3
LIBGEOM=../libgeom
JETWRAP=../jet-wrapper
EXTRACT=../bezier-extractions

CXXFLAGS=-fsanitize=address -Wall -pedantic -std=c++17 -O0 -g \
		 -I$(LIBGEOM) -I$(JETWRAP) -I$(EXTRACT) -I$(EIGEN)
LDFLAGS=-lasan -L. -L$(LIBGEOM)/debug -lgeom -lomp

quadfit-test: quadfit-test.o libquadfit.a
	$(CXX) -o $@ $< -lquadfit $(LDFLAGS)

libquadfit.a: quadfit.o bspline-fit.o io.o discrete-mask.o \
              fit-ribbon.o connect-g1.o multiply.o \
              $(EXTRACT)/bezier-extractions.o $(JETWRAP)/jet-wrapper.o
	$(AR) r -o $@ $^

ribbon-test: ribbon-test.o fit-ribbon.o multiply.o io.o $(EXTRACT)/bezier-extractions.o
	$(CXX) -o $@ $^ $(LDFLAGS)

connect-test: connect-test.o connect-g1.o multiply.o io.o $(EXTRACT)/bezier-extractions.o
	$(CXX) -o $@ $^ $(LDFLAGS)
