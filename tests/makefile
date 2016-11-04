# Makefile for the TESTS folder.

include ../makefile_inc

all: t01

EXE = \
	mtk_1d_div_test

# Compiling statements:

%.o: %.cc $(INCLUDE)/mtk.h
	$(CC) $(CCFLAGS) $<

# Linking statements:
t01: mtk_1d_div_test.o mtk_1d_div_test.cc $(INCLUDE)/mtk.h
	$(LINKER) mtk_1d_div_test.o $(LIBS) -o mtk_1d_div_test
	@echo ' '

clean:
	rm -f *~ *.o $(EXE)
