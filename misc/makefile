#https://gist.github.com/zachwhaley/9458612
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic -funroll-loops -finline-functions #-march=native -std=c++14 
LDLIBS =

PRGM  = tensorReducer
SRCS := fullgeneration.cpp vertexData.cpp #$(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)

.PHONY: all clean

# By default, make begins by processing the first target that does not begin with a . aka the default goal
all: $(PRGM)

debug: CXXFLAGS += -pg -ggdb3
debug: $(PRGM)
debug:
	$(info Now optionally run:)
	$(info time valgrind --tool=callgrind --dump-instr=yes --collect-jumps=yes ./tensorReducer 7-25-8-26-9-27-10-28-11-29-12-30 gh /tmp/binbry.test)
	$(info callgrind_annotate --show-percs=yes --tree=both --auto=yes callgrind.out.20269 > tree.20269)
	$(info WHEREVER/gprof2dot.py -f callgrind callgrind.out.20269 | dot -Tsvg -o output.svg)
	$(info kcachegrind callgrind.out.WHATEVERNUMBER)

$(PRGM): $(OBJS)
	$(CXX) $(OBJS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf $(OBJS) $(DEPS) $(PRGM)

-include $(DEPS)
