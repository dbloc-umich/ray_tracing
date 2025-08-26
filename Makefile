# Top-level directories
BIN_DIR = bin
DEP_DIR = dep
INC_DIR = include
OBJ_DIR = obj
SRC_DIR = src
# TEST_DIR = test

# Subdirecotires
INC_SUBDIRS = $(shell find $(INC_DIR) -type d)
SRC_SUBDIRS = $(shell find $(SRC_DIR) -type d)

SRC = $(shell find $(SRC_DIR) -name '*.cpp')
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
DEP = $(SRC:$(SRC_DIR)/%.cpp=$(DEP_DIR)/%.d)
# TEST = $(wildcard $(TEST_DIR)/*.cpp)
TARGET = $(BIN_DIR)/main.exe

INC = $(addprefix -I,$(INC_SUBDIRS))
CXX = g++
CXXFLAGS = -Wall -O2 -MMD -MP $(INC)

all: $(TARGET)

$(TARGET): $(OBJ)
	@mkdir -p $(BIN_DIR)
	@$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@) $(DEP_DIR)/$(dir $*)
	@$(CXX) $(CXXFLAGS) -c $< -o $@ -MF $(DEP_DIR)/$*.d

-include $(DEP)

clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR) $(DEP_DIR)

run: $(TARGET)
	@./$(TARGET)