LINKLIB=  -L/usr/local/lib/ -lm

#LINKLIB= -lm -lc -L/usr/local/lib/ -lfftw3

INCLUDE=-I/usr/local/include/


CFLAGS=-g
CC=gcc


random: random.o  
	$(CC) $(CFLAGS) -o random random.o $(LINKLIB)
	rm -rf *.o
	rm -rf *~

clean:
	rm -rf *.o
	rm -rf *~
