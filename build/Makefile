# Compiler
CXX := g++
CXXFLAGS := -std=c++14 -fsanitize=address -Wall -Wextra -Wpedantic -Wconversion -Wunused-parameter -Wunused-but-set-parameter -O2

#-fsanitize=address,undefined to check memory leaks

# Directories
ROOTDIR := ..
BUILDDIR := $(ROOTDIR)/build
SRCDIR := $(ROOTDIR)/src
BINDIR := $(ROOTDIR)/build
INCDIR := $(ROOTDIR)/inc
OBJDIR := $(BUILDDIR)/obj
OUTDIR := $(ROOTDIR)/output

# Files
SOURCES := $(wildcard $(SRCDIR)/*.cpp) #$(wildcard $(INCDIR)\chi2/*.cpp)
OBJECTS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

$(info $(SRCDIR))

# Program name
PROGRAM_NAME := test.exe

# Target binary
TARGET := $(BINDIR)/$(PROGRAM_NAME)

# Default target
all: $(TARGET)

# Linking
$(TARGET): $(OBJECTS) | $(OUTDIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

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
	rm -rf $(OBJDIR) $(TARGET) $(OUTDIR)

# Phony targets
.PHONY: all clean
