CXX = g++
CXXFLAGS = -Wall -Wextra --std=c++14 -MMD -MP -Iinclude

FILES = main Point Direction Shape Sphere Box BoundingBox Octree RayTracing Vector
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
DEP_DIR = dep

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC))
DEP = $(patsubst $(SRC_DIR)/%.cpp, $(DEP_DIR)/%.d, $(SRC))
TARGET = $(BIN_DIR)/main.exe

all: $(TARGET)

$(TARGET): $(OBJ)
	@mkdir -p $(BIN_DIR)
	@$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR) $(DEP_DIR)
	@$(CXX) $(CXXFLAGS) -c $< -o $@
	@mv $(basename $@).d $(DEP_DIR)/

-include $(DEP)

clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR) $(DEP_DIR)

run: $(TARGET)
	@./$(TARGET)