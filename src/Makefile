# Usage:
# make        # compile all binary
# make clean  # remove ALL binaries and objects

.PHONY: all clean

CC = g++

LINKERFLAG = -lm

SRCS := HCV_Model.cpp main.cpp
BINS := $(SRCS:%.cpp=%)

all: ${BINS}

%: %.o
	@${CC} ${LINKERFLAG} $< -o $@

%.o: %.cpp 
	@${CC} -c $< 

clean:
	@rm *.o ${BINS}