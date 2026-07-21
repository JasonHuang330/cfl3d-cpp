CXX = g++
CGNS_INC = /opt/homebrew/include
CGNS_LIB = /opt/homebrew/lib

# Sources use flat  #include "foo.h" , so every source folder must be on the
# include path. Root (.) is needed for  #include "runtime/..." .
INCDIRS  = . runtime $(shell find src -type d) $(CGNS_INC)
CXXFLAGS = -O2 -std=c++17 $(addprefix -I,$(INCDIRS)) -MMD -MP
LDFLAGS  = -L$(CGNS_LIB) -lcgns

TARGET = main
SRCS   = $(shell find src runtime -name '*.cpp')
OBJS   = $(SRCS:.cpp=.o)
DEPS   = $(OBJS:.o=.d)

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)
