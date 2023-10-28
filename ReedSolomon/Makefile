
CC := gcc
CFLAGS := -Wall -Wextra -g
SRC_DIR := src
OBJ_DIR := obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.c)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRC_FILES))
EXECUTABLES := genRS_LUT benchmark test_math
COMMON_OBJS := obj/rs_gf8.o obj/gf8.o obj/rs_gf16.o obj/gf16.o
COMMON_HEADERS := src/rs_gf8.h src/gf8.h src/rs_gf16.h src/gf16.h
COMMON_DEPS := $(COMMON_OBJS) $(COMMON_HEADERS)

.PHONY: all
all: $(EXECUTABLES)

genRS_LUT: $(OBJ_DIR)/genRS_LUT.o $(COMMON_DEPS)
	$(CC) $(CFLAGS) -o $@ $^

test_math: $(OBJ_DIR)/test_math.o $(COMMON_DEPS)
	$(CC) $(CFLAGS) -o $@ $^

benchmark: $(OBJ_DIR)/benchmark.o $(COMMON_DEPS)
	$(CC) $(CFLAGS) -O3 -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -f *.bin $(OBJ_DIR)/* *.exe