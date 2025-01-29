TARGETS = hacoo_test candecomp

CC ?= clang
override CFLAGS += -g -Wno-everything -pthread -lm `pkg-config --cflags cunit` -fopenmp
LDLIBS+=-lm -fopenmp `pkg-config --libs cunit`

HACOO = hacoo.o cpd.o matrix.o mttkrp.o

all: hacoo_test candecomp

hacoo_test: hacoo_test.o $(HACOO)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

candecomp: candecomp.o $(HACOO) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(TARGETS) *.o
