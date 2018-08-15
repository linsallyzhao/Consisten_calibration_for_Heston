CC = g++
CFLAGS = -Wall -Wextra -Wno-unused-parameter  -O2 -flto -march=native -mtune=native
LIBS = -lm -llevmar -lfaddeeva
DEBUGFLAGS = -fno-omit-frame-pointer -g3
SANITISERS = -fsanitize=undefined -fsanitize=address -fsanitize=pointer-subtract -fsanitize=pointer-compare

all: HC CFG DIS

# -fprofile-* currently off due to wall time of run\

HC: HC.cpp
	$(CC) $(CFLAGS) $(LIBS) HC.cpp -o HC

CFG: CFG.cpp
	$(CC) $(CFLAGS) $(LIBS) CFG.cpp -o CFG

DIS: DIS.cpp
	$(CC) $(CFLAGS) $(LIBS) DIS.cpp -o DIS

clean:
	rm -f *~ *.o HC CFG DIS perf.diff

perf:
	$(CC) $(CFLAGS) $(DEBUGFLAGS) $(LIBS) DIS.cpp -o DIS
	rm perf.data
	perf record -g ./DIS >/dev/null
	perf script | c++filt | gprof2dot -f perf --colour-nodes-by-selftime -n 0.1 | dot -Tpdf -o /tmp/cfg.pdf
	mupdf /tmp/cfg.pdf
	rm /tmp/cfg.pdf
	#perf report -g

diff: all
	./compare_outputs

sanitise:
	$(CC) $(CFLAGS) $(DEBUGFLAGS) $(SANITISERS) $(LIBS) DIS.cpp -o DIS
