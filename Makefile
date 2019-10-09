EXE	= main.out
DBG = gdb

SRC_DIR = src
OBJ_DIR = obj
PYTHON_DIR = python
TEST_DIR = tests

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

CPPFLAGS += -Iinclude/.. 
CDBGFLAGS = -g -O0
CRLSFLAGS = -O3
CFLAGS += -Wall -Wextra -Wconversion -Wno-unused -pedantic
LDFLAGS += -Llib
LDLIBS += -lm

ifndef no_user
PYFLAGS = --user
endif

all: release

release:
	OPT="$(CRLSFLAGS)" make $(EXE)

$(EXE): $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CPPFLAGS) $(OPT) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJ)

fresh:
	make clean
	make all

run:
	./$(EXE)

debug:
	make clean
	OPT="$(CDBGFLAGS)" make $(EXE)
	$(DBG) -q ./$(EXE)

check:
	make all
	valgrind --leak-check=full ./$(EXE)

py_build: $(SRC_DIR)/*.c $(PYTHON_DIR)/*.c
	python3 buildext.py build
	python3 buildext.py install $(PYFLAGS)
	python3 -q -X faulthandler $(TEST_DIR)/test_ext.py

py_debug:
	gdb -q -ex start --args python3 $(TEST_DIR)/test_ext.py

.PHONY: all clean