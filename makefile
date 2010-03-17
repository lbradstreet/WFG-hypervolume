EXE = wfg
OBJ = read.o avl.o

OPT = -O3 -march=nocona
#OPT = -g
#CC = gcc -std=c99 -Wall -pedantic -Werror -lm -O5
CC = gcc -std=c99 -Wall  $(OPT)

$(EXE): $(OBJ) wfg.c
#	$(CC) -Dopt=0 -Dmode=0 -o wfg0 wfg.c  $(OBJ)
#	$(CC) -Dopt=1 -Dmode=0 -o wfg1 wfg.c $(OBJ)
#	$(CC) -Dopt=2 -Dmode=0 -o wfg2 wfg.c $(OBJ)
#	$(CC) -Dopt=3 -Dmode=0 -o wfg3 wfg.c $(OBJ)
#	$(CC) -Dopt=4 -Dmode=0 -o wfg4 wfg.c $(OBJ)
	$(CC) -Dopt=0 -Dmode=1 -o exc1 wfg.c $(OBJ)
	$(CC) -Dopt=0 -Dmode=2 -o exc2 wfg.c $(OBJ)
	$(CC) -Dopt=0 -Dmode=3 -o exc3 wfg.c $(OBJ)
#	$(CC) -Dopt=0 -Dmode=4 -o exc4 wfg.c $(OBJ)

%.o: %.c
	$(CC) -c $<

clean: 
	rm wfg-* *.o

