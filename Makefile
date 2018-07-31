CC = g++
CFLAGS = -Wall -Werror -Wextra
LIBS = -lm -llevmar

all: HC CFG

HC: HC.cpp
	$(CC) $(CFLAGS) $(LIBS) HC.cpp -o HC

CFG: CFG.cpp
	$(CC) $(CFLAGS) $(LIBS) CFG.cpp -o CFG

clean:
	rm -f *~ *.o HC CFG

diff: all
	./compare_outputs
