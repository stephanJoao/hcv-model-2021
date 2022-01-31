# Usage:
# make        # compile all binary
# make clean  # remove ALL binaries and objects

# .PHONY: all clean

#* Variables

# Name of the project
PROJ_NAME = program

# .cpp files
CPP_SOURCE = ${wildcard *.cpp}

# .hpp files
HPP_SOURCE = ${wildcard *.hpp}

# Object files
OBJ = ${CPP_SOURCE:.c=.o}

# Compiler
CC = g++

# Compiler flags
CC_FLAGS = -c -W -Wall -ansi -pedantic


# SRCS := main.cpp HCV_Model.cpp
# BINS := ${SRCS:%.cpp=%}

#* Compilation and linking

all: ${PROJ_NAME} clean

${PROJ_NAME}: ${OBJ}
	@${CC} -o $@ $^

%.o: %.cpp %.hpp
	@${CC} -o $@ $< ${CC_FLAGS}

main.o: main.cpp ${H_SOURCE}
	@${CC} -o $@ $< ${CC_FLAGS}

# main.o: main.cpp HCV_Model.hpp
# 	@${CC} -o main.o main.cpp -c -W -Wall -ansi -pedantic

# HCV_Model.o: HCV_Model.cpp HCV_Model.hpp GaussLegendreQuadrature.h
# 	@${CC} -o HCV_Model.o HCV_Model.cpp -c -W -Wall -ansi -pedantic

clean:
	@rm -rf *.o