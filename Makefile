EXECUTABLE=DoAnalysis

SRC_DIR=src
OBJ_DIR=obj

SRC=$(wildcard $(SRC_DIR)/*.cc)
OBJ=$(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

INC=include analyses
INC_PARAMS=$(foreach d, $(INC), -I$d)

CC=`root-config --cxx`
CFLAGS=-c -g -Wall `root-config --cflags` $(INC_PARAMS)
LDFLAGS=-lTreePlayer `root-config --ldflags`
LDLIBS=-lTreePlayer `root-config --glibs`

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)
