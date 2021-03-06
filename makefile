EXE = wfg
OBJ = read.o

OPT = -O3 



march_error = $(error please define an architecture, e.g., 'make march=pentium')
ifdef ARCH
OPT=$(ARCH) $(COPT)
else
ifndef march 
  $(march_error)
  endif
  OPT += -march=$(march)
endif


#OPT = -g
#CC = gcc -std=c99 -Wall -pedantic -Werror -lm -O5
CC = gcc -std=c99 -Wall  $(OPT)

$(EXE): $(OBJ) wfg.c
	$(CC) -Dopt=0 -o wfg0 wfg.c $(OBJ)
	$(CC) -Dopt=1 -o wfg1 wfg.c $(OBJ)
	$(CC) -Dopt=2 -o wfg2 wfg.c $(OBJ)

%.o: %.c
	$(CC) -c $<

clean: 
	rm -f wfg[0-4] exc[0-4] *.o 
	rm -rf *.dSYM

VERSION=$(shell git tag)

release: clean
	rm -rf '/tmp/WFG_'$(VERSION)
	mkdir -p '/tmp/WFG_'$(VERSION)
	cp -rf * '/tmp/WFG_'$(VERSION)
	cd /tmp/; tar czvf $(PWD)/../'WFG_'$(VERSION).tar.gz 'WFG_'$(VERSION)
