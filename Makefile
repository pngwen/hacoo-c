TARGETS = hacoo_test candecomp

CC ?= clang
CFLAGS += -g -Wno-everything -pthread -lm -fopenmp
CFLAGS += `pkg-config cunit --cflags`
LDLIBS=-lm 
LDLIBS+= `pkg-config cunit --libs`

HACOO = hacoo.o cpd.o matrix.o mttkrp.o

all: hacoo_test candecomp

hacoo_test: hacoo_test.o $(HACOO)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

candecomp: candecomp.o $(HACOO) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(TARGETS) *.o
