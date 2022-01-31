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
OBJ = ${subst .cpp,.o,${subst source,objects,${CPP_SOURCE}}}

# Compiler
CC = g++

# Compiler flags
CC_FLAGS = -c -W -Wall -ansi -pedantic


#* Compilation and linking

all: objectsFolder ${PROJ_NAME}

${PROJ_NAME}: ${OBJ}
	@${CC} $^ -o $@

./objects/%.o: ./source/%.cpp ./source/%.hpp
	@${CC} $< ${CC_FLAGS} -o $@

./objects/main.o: ./source/main.cpp ${H_SOURCE}
	@${CC} $< ${CC_FLAGS} -o $@

objectsFolder:
	@mkdir -p objects

clean:
	@rm -rf ./objects/*.o

.PHONY: all clean