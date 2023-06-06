#https://gist.github.com/zachwhaley/9458612
CXXFLAGS = -g -O3 -Wall -Wextra -pedantic -funroll-loops -finline-functions #-march=native -std=c++14 
LDLIBS =

PRGM  = tensorReducer
SRCS := fullgeneration.cpp vertexData.cpp #$(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)

.PHONY: all clean

all: $(PRGM)

$(PRGM): $(OBJS)
	$(CXX) $(OBJS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf $(OBJS) $(DEPS) $(PRGM)

-include $(DEPS)
