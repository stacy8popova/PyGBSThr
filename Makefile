CC = g++

CFLAGS = -lm -pthread -O3 -march=native -funroll-loops -Wall -Wextra -Wpedantic
BUILDDIR := build
OBJDIR := $(BUILDDIR)

OBJ := $(OBJDIR)/Minors.o $(OBJDIR)/Exact.o
HEADERS := headers.h matrix.h
MODULES := Minors Exact

all: dir $(OBJ) $(MODULES)
dir :
	mkdir -p $(BUILDDIR)
Minors : $(OBJDIR)/Minors.o 
	$(CC) $^ -o $@ $(CFLAGS)
Exact : $(OBJDIR)/Exact.o 
	$(CC) $^ -o $@ $(CFLAGS)
$(OBJDIR)/%.o : %.cpp $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)
.PHONY: clean
clean:
	rm -rf $(BUILDDIR)
	rm $(MODULES)