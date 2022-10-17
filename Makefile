all: quadfit-test ribbon-test connect-test

EIGEN=/usr/include/eigen3
LIBGEOM=../libgeom
JETWRAP=../jet-wrapper
EXTRACT=../bezier-extractions

# Release
CXXFLAGS=-Wall -pedantic -std=c++17 -O3 -DNDEBUG \
		 -I$(LIBGEOM) -I$(JETWRAP) -I$(EXTRACT) -I$(EIGEN)
LDFLAGS=-L. -L$(LIBGEOM)/release -lgeom -lomp

# Debug
# CXXFLAGS=-Wall -pedantic -std=c++17 -O0 -g -DDEBUG \
# 		 -I$(LIBGEOM) -I$(JETWRAP) -I$(EXTRACT) -I$(EIGEN)
# LDFLAGS=-L. -L$(LIBGEOM)/debug -lgeom -lomp

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

.PHONY: clean
clean:
	$(RM) *.o libquadfit.a quadfit-test ribbon-test connect-test
