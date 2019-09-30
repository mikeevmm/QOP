OBJS	= main.o gate.o util.o option.o circuit.o vector.o iter.o
SOURCE	= main.c gate.c util.c option.c circuit.c vector.c iter.c
HEADER	= main.h gate.h util.h option.h circuit.h vector.h iter.h
OUT	= main.out
CC	 = gcc
FLAGS	 = -g -O1 -c -Wall -Wextra -Wconversion -pedantic
LFLAGS	 = -lm -lblas -llapack -llapacke

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

clean:
	rm -f $(OBJS) $(OUT)

run: $(OUT)
	./$(OUT)

debug: $(OUT)
	gdb ./$(OUT)