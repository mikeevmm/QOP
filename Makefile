OBJS	= main.o gate.o util.o option.o circuit.o vector.o
SOURCE	= main.c gate.c util.c option.c circuit.c vector.c
HEADER	= main.h gate.h util.h option.h circuit.h vector.h
OUT	= main.out
CC	 = gcc
FLAGS	 = -g -O0 -c -Wall -Wextra -Wconversion -pedantic
LFLAGS	 = -lm -lblas -llapack -llapacke

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

clean:
	rm -f $(OBJS) $(OUT)

run: $(OUT)
	./$(OUT)

debug: $(OUT)
	gdb ./$(OUT)