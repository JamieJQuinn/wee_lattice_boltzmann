CC=clang++
CFLAGS=-std=c++17 -O3 -fopenmp -march=native -DNDEBUG -I$(INCLUDE_DIR)
LDFLAGS=-fopenmp
SRC_DIR=src
BUILD_DIR=build
INCLUDE_DIR=include

SOURCES=$(wildcard $(SRC_DIR)/*.cpp)
HEADERS=$(wildcard $(INCLUDE_DIR)/*.hpp)
OBJECTS=$(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
EXECUTABLE=exe

.PHONY: all
all: $(BUILD_DIR)/$(EXECUTABLE)

# .PHONY: flags
# flags:
# 	echo $(CFLAGS) > compile_flags.txt

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/$(EXECUTABLE): $(BUILD_DIR) $(OBJECTS) $(HEADERS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
