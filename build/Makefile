# Directories
ROOTDIR := ..
BUILDDIR := $(ROOTDIR)/build
SRCDIR := $(ROOTDIR)/src
BINDIR := $(ROOTDIR)/build
INCDIR := $(ROOTDIR)/inc
OBJDIR := $(BUILDDIR)/obj

# Compiler
CXX := g++
MPICXX := mpic++

# Compiler flags
CXXFLAGS := -std=c++14 -Wall -Wextra -Wpedantic -Wconversion -Wunused-parameter -Wunused-but-set-parameter -O2 -I$(INCDIR)

# Files
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

# Program name
PROGRAM_NAME := test
#ifdef NAME
#	PROGRAM_NAME := $(NAME)
#endif

# Target binary
TARGET := $(BINDIR)/$(PROGRAM_NAME)
MPI_TARGET := $(BINDIR)/$(PROGRAM_NAME)_mpi

# Default target
all: $(TARGET)

# MPI target
mpi: CXX := $(MPICXX)
mpi: $(MPI_TARGET)

# Linking
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(MPI_TARGET): $(OBJECTS)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

# Compilation
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create build directory
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Clean
clean:
	rm -rf $(OBJDIR) $(TARGET) $(MPI_TARGET)

# Phony targets
.PHONY: all clean mpi
