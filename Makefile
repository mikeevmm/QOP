EXE	= main.out
DBG = gdb

SRC_DIR = src
OBJ_DIR = obj

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

CPPFLAGS += -Iinclude
CFLAGS += -g -O0 -Wall -Wextra -Wconversion -pedantic
LDFLAGS += -Llib
LDLIBS += -lm -lblas -llapack -llapacke

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJ)

fresh:
	make clean
	make all

run:
	./$(EXE)

debug:
	$(DBG) ./$(EXE)

check:
	make all
	valgrind --leak-check=full ./$(EXE)

.PHONY: all clean