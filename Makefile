# ============================================================================
# Makefile — 多平台构建
# ----------------------------------------------------------------------------
#   make            构建（-O3 -march=native）
#   make debug      调试构建（-O0 -g，带 AddressSanitizer）
#   make run        构建并运行（默认参数）
#   make clean      清理产物
# ============================================================================

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -march=native -pthread -Wall -Wextra
LDFLAGS  ?= -pthread

TARGET   := rare_sample_efficient
SRCS     := main.cpp sampler.cpp recorder.cpp thresholds.cpp shared_state.cpp
HDRS     := config.hpp thresholds.hpp rng.hpp scompute.hpp \
            shared_state.hpp sampler.hpp recorder.hpp
OBJS     := $(SRCS:.cpp=.o)

.PHONY: all debug run clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

debug: CXXFLAGS := -std=c++17 -O0 -g -pthread -Wall -Wextra -fsanitize=address,undefined
debug: LDFLAGS  := -pthread -fsanitize=address,undefined
debug: clean $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
ifeq ($(OS),Windows_NT)
	-del /Q $(OBJS) $(TARGET) $(TARGET).exe 2>NUL
else
	rm -f $(OBJS) $(TARGET) $(TARGET).exe
endif
