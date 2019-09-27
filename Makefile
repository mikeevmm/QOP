OBJS	= main.o gate.o util.o option.o circuit.o vector.o
SOURCE	= main.c gate.c util.c option.c circuit.c vector.c
HEADER	= main.h gate.h util.h option.h circuit.h vector.h
OUT	= main
CC	 = gcc
FLAGS	 = -g -c -Wall
LFLAGS	 = -lm -lblas -llapack -llapacke

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

main.o: main.c
	$(CC) $(FLAGS) main.c

gate.o: gate.c
	$(CC) $(FLAGS) gate.c

util.o: util.c
	$(CC) $(FLAGS) util.c


clean:
	rm -f $(OBJS) $(OUT)

run: $(OUT)
	./$(OUT)