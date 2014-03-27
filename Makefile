CC = g++
CFLAGS = -O2 -Wall
LIB = -lm

SRC_DIR = src
OBJ_DIR = obj

EXEC = SolarSystem

OBJ_FILES = $(addprefix $(OBJ_DIR)/, solarsystem.o orbital.o)

$(EXEC): $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $(OBJ_FILES)

$(OBJ_DIR)/orbital.o: $(SRC_DIR)/orbital.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJ_DIR)/solarsystem.o: $(SRC_DIR)/solarsystem.cpp $(SRC_DIR)/parameter_readin.h $(SRC_DIR)/solarsystem.h $(SRC_DIR)/orbital.o
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f *.o