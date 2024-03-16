LINKLIB=  -L/usr/local/lib/  -fopenmp -lfftw3f_omp -lfftw3f -lm

#LINKLIB= -lm -lc -L/usr/local/lib/ -lfftw3

INCLUDE=-I/usr/local/include/


CFLAGS=-g
CC=gcc

lightcone: lightcone.o allotarrays.o lightcone_funcs.o funcs.o
	$(CC) $(CFLAGS) -o lightcone  lightcone.o allotarrays.o lightcone_funcs.o funcs.o $(LINKLIB)
	rm -rf *.o
	rm -rf *~

lightcone.o:	lightcone.c
	$(CC) -c $(CFLAGS) $(INCLUDE) lightcone.c $(LINKLIB)

lightcone_funcs.o:	lightcone_funcs.c
	$(CC) -c $(CFLAGS) $(INCLUDE) lightcone_funcs.c $(LINKLIB) 

allotarrays.o:	allotarrays.c
	$(CC) -c $(CFLAGS) $(INCLUDE) allotarrays.c $(LINKLIB)

funcs.o:	funcs.c
	$(CC) -c  $(CFLAGS)   $(INCLUDE) funcs.c $(LINKLIB)

clean:
	rm -rf *.o
	rm -rf *~
