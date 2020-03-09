LLC:=/usr/local/opt/llvm@6/bin/llc
SEQLIB:=/Users/inumanag/Projekti/seq/devel/build

all: haptreex

haptreex: haptreex.o
	clang -L$(SEQLIB) -lomp -lseqrt -o haptreex haptreex.o

haptreex.o: haptreex.bc
	$(LLC) haptreex.bc -filetype=obj -o haptreex.o

haptreex.bc: main.seq
	DYLD_LIBRARY_PATH=$(SEQLIB); SEQ_PATH=$(SEQLIB)/../stdlib; $(SEQLIB)/seqc -o haptreex.bc main.seq
