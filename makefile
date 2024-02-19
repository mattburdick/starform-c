# Makefile for "starform", a star system and planet generator

CFLAGS = -g
OBJS = starform.o accrete.o enviro.o stars.o display.o utils.o
LIBS = -lm
SHARFILES = README makefile.msc makefile.tc makefile starform.c accrete.c \
	enviro.c stars.c display.c utils.c const.h structs.h config.h protos.h


.c: const.h config.h structs.h protos.h
	$(CC) $(CFLAGS) -c $(<).c

starform: $(OBJS)
	$(CC) $(LDFLAGS) -o starform $(OBJS) $(LIBS)
	@echo "starform made"

clean:
	rm -f *.o *.ln

clobber:
	rm -f *.o *.ln starform

lint:
	lint -abchp starform.c accrete.c enviro.c stars.c display.c utils.c

shar: $(SHARFILES)
	shar -abcCs $(SHARFILES) >starform.shar

