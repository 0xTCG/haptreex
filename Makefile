all:
	seqc -o haptreex.bc main.seq
	/Users/inumanag/Projekti/deps/Tapir-LLVM/build/bin/llc haptreex.bc -filetype=obj -o haptreex.o
	clang -L /Users/inumanag/Projekti/seq/devel/build -lomp -lseqrt -o haptreex haptreex.o
