ALL=vocab docpd-hacoo doccoo
LD=g++
CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -Wall -Wextra -O2 -I../ -g
LDLIBS=-lm -fopenmp ../../cunit-local/lib/libcunit.a -lopenblas

all: $(ALL)
vocab: vocab.o tokenize.o
	$(CXX) $(CXXFLAGS) -o $@ $^
docpd-hacoo: docpd-hacoo.o tokenize.o ../alto.o ../hacoo.o ../matrix.o ../cpd.o ../mttkrp.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)
doccoo: doccoo.o tokenize.o ../alto.o ../hacoo.o ../matrix.o ../cpd.o ../mttkrp.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(ALL) *.o