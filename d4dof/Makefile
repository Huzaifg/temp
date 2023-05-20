# Specify compiler and flags
CXX = g++
CXXFLAGS = -Wall -O3 -DADEPT_RECORDING_PAUSABLE
#Run time seach path for the shared library of adept
LDFLAGS = -Wl,-rpath -Wl,/usr/local/lib
# Normal object files to be created
OBJECTS = fourDof.o main.o
# Program name
PROGRAM = fourDof
# Include-file location
INCLUDES = -I/usr/local/include
# Library location and name, plus the math library
LIBS = -L/usr/local/lib -lm -ladept
# Rule to build the program (typing "make" will use this rule)
$(PROGRAM): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $(PROGRAM)
# Rule to build a normal object file (used to compile all objects in OBJECTS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) -c $<