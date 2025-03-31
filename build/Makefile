# Directories
ROOTDIR := ..
BUILDDIR := $(ROOTDIR)/build
SRCDIR := $(ROOTDIR)/src
BINDIR := $(ROOTDIR)/build
INCDIR := $(ROOTDIR)/inc
OBJDIR := $(BUILDDIR)/obj
OUTDIR := $(ROOTDIR)/output

# Compiler
CXX := g++
MPICXX := mpic++

# Compiler flags
CXXFLAGS := -std=c++14 -Wall -Wextra -Wpedantic -Wconversion -Wunused-parameter -Wunused-but-set-parameter -O2 -I$(INCDIR)

# other flag to be considered -fsanitize=address,undefined

# Files
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

# Program name
PROGRAM_NAME := simulation

# Target binary
TARGET := $(BINDIR)/$(PROGRAM_NAME).exe
MPI_TARGET := $(BINDIR)/$(PROGRAM_NAME)_mpi.exe

# Default target
all: $(TARGET)

# MPI target
mpi: CXX := $(MPICXX)
mpi: $(MPI_TARGET)

# Linking
$(TARGET): $(OBJECTS) | $(OUTDIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(MPI_TARGET): $(OBJECTS) | $(OUTDIR)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

# Compilation
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

	
# Create build directory
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OUTDIR):
	mkdir -p $(OUTDIR)

# Clean
clean:

	rm -rf $(OBJDIR) $(TARGET) $(MPI_TARGET) $(OUTDIR)

# Phony targets
.PHONY: all clean mpi
