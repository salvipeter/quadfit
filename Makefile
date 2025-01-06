all: quadfit-test ribbon-test connect-test knots-test curves2obj

EIGEN=/usr/include/eigen3
LIBGEOM=../libgeom
JETWRAP=../jet-wrapper
EXTRACT=../bezier-extractions
CONNECT=../bspline-connect

# Release
CXXFLAGS=-Wall -pedantic -std=c++17 -O3 -g -DNDEBUG \
		 -I$(LIBGEOM) -I$(JETWRAP) -I$(EXTRACT) -I$(EIGEN) -I$(CONNECT)
LDFLAGS=-L. -L$(LIBGEOM)/release -lgeom -lomp -lgmp

# Debug
# CXXFLAGS=-Wall -pedantic -std=c++17 -O0 -g -fsanitize=address -DDEBUG \
# 		 -I$(LIBGEOM) -I$(JETWRAP) -I$(EXTRACT) -I$(EIGEN) -I$(CONNECT)
# LDFLAGS=-L. -L$(LIBGEOM)/debug -lasan -lgeom -lomp -lgmp

quadfit-test: quadfit-test.o libquadfit.a
	$(CXX) -o $@ $< -lquadfit $(LDFLAGS)

libquadfit.a: quadfit.o bspline-fit.o io.o discrete-mask.o closest-point.o knots.o \
              fit-ribbon.o connect-g1.o multiply.o switches.o \
              $(EXTRACT)/bezier-extractions.o $(JETWRAP)/jet-wrapper.o \
              $(CONNECT)/bsconnect.o
	$(AR) r -o $@ $^

ribbon-test: ribbon-test.o fit-ribbon.o multiply.o io.o $(EXTRACT)/bezier-extractions.o
	$(CXX) -o $@ $^ $(LDFLAGS)

connect-test: connect-test.o connect-g1.o multiply.o io.o $(EXTRACT)/bezier-extractions.o
	$(CXX) -o $@ $^ $(LDFLAGS)

knots-test: knots-test.o knots.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

curves2obj: curves2obj.cc
	$(CXX) -o $@ $< $(CXXFLAGS) $(LDFLAGS)

.PHONY: clean
clean:
	$(RM) *.o libquadfit.a quadfit-test ribbon-test connect-test knots-test curves2obj
