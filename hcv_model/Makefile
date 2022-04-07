# Usage:
# make        # compile all binary
# make clean  # remove all binaries and objects

#* Variables

# Name of the project
PROJ_NAME = program

# .cpp files
CPP_SOURCE = ${wildcard ./source/*.cpp}

# .hpp files
HPP_SOURCE = ${wildcard ./source/*.hpp}

# Object files
OBJ = ${subst .cpp,.o,${subst source,build,${CPP_SOURCE}}}

# Compiler
CC = g++

# Compiler flags
CC_FLAGS = -c -W -Wall -ansi -pedantic -O3


#* Compilation and linking

all: buildFolder ${PROJ_NAME}

${PROJ_NAME}: ${OBJ}
	@${CC} $^ -o $@

./build/%.o: ./source/%.cpp ./source/%.hpp
	@${CC} $< ${CC_FLAGS} -o $@

./build/main.o: ./source/main.cpp ${H_SOURCE}
	@${CC} $< ${CC_FLAGS} -o $@

buildFolder:
	@mkdir -p build

clean:
	@rm -rf ./build/*.o

run:
	@./program

.PHONY: all clean