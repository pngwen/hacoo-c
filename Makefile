all: main

CC = gcc
override CFLAGS += -g -fopenmp -pthread -Wno-everything -I../cunit-local/include -I. -mbmi2
CXXFLAGS := $(CFLAGS)
LDLIBS=-lm -fopenmp ../cunit-local/lib/libcunit.a -lopenblas

SRCS = $(shell find . -name '.ccls-cache' -type d -prune -o -type f -name '*.c' -print)
HEADERS = $(shell find . -name '.ccls-cache' -type d -prune -o -type f -name '*.h' -print)

main: $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) $(SRCS) -o "$@" $(LDLIBS)

main-debug: $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) -O0 $(SRCS) -o "$@" $(LDLIBS)

mttkrp_test: hacoo.o mttkrp_test.o matrix.o mttkrp.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

candecomp: candecomp.o hacoo.o matrix.o cpd.o mttkrp.o alto.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

hacoo_test: hacoo_test.o hacoo.o matrix.o cpd.o mttkrp.o alto.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

matrix_op_test: matrix_op_test.o matrix.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

alto_encode_test: alto_encode_test.o alto.o hacoo.o matrix.o cpd.o mttkrp.o
	g++ $(CFLAGS) -o $@ $^ $(LDLIBS)

hacoo_mttkrp: hacoo_mttkrp.o alto.o hacoo.o matrix.o cpd.o mttkrp.o
	g++ $(CFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f main main-debug mttkrp_test *.o
