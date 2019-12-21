LLC:=/usr/local/opt/llvm@6/bin/llc
SEQLIB:=/Users/inumanag/Desktop/Projekti/seq/devel/build 

all:
	seqc -o haptreex.bc main.seq
	$(LLC) haptreex.bc -filetype=obj -o haptreex.o
	clang -L$(SEQLIB) -lomp -lseqrt -o haptreex haptreex.o
