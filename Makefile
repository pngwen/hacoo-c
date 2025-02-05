TARGETS = hacoo_test candecomp

CC = clang
override CFLAGS += -g -Wno-everything -pthread -lm -fopenmp
CFLAGS += -I/home/linuxbrew/.linuxbrew/Cellar/cunit/2.1-3/include/
LDLIBS=-lm -lcunit
LDLIBS += -L/home/linuxbrew/.linuxbrew/Cellar/cunit/2.1-3/lib

HACOO = hacoo.o cpd.o matrix.o mttkrp.o

all: hacoo_test candecomp

hacoo_test: hacoo_test.o $(HACOO)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

candecomp: candecomp.o $(HACOO) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(TARGETS) *.o