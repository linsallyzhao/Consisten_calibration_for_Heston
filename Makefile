CC = g++
CFLAGS = -Wall -Werror -Wextra -O2 -fno-omit-frame-pointer -g3
LIBS = -lm -llevmar

all: HC CFG

HC: HC.cpp
	$(CC) $(CFLAGS) $(LIBS) HC.cpp -o HC

CFG: CFG.cpp
	$(CC) $(CFLAGS) $(LIBS) CFG.cpp -o CFG

clean:
	rm -f *~ *.o HC CFG perf.diff

perf: CFG
	rm perf.data
	perf record -g ./CFG >/dev/null
	perf script | c++filt | gprof2dot -f perf | dot -Tpdf -o /tmp/cfg.pdf
	mupdf /tmp/cfg.pdf
	rm /tmp/cfg.pdf

diff: all
	./compare_outputs
