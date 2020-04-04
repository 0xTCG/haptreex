LLC := llc
CLANG := clang

all: haptreex

haptreex: haptreex.o
	$(CLANG) -lomp -lseqrt -o build/haptreex build/haptreex.o

haptreex.o: haptreex.bc
	$(LLC) build/haptreex.bc -filetype=obj -o build/haptreex.o

haptreex.bc: src/main.seq
	mkdir -p build
	seqc -o build/haptreex.bc src/main.seq

clean:
	rm -rf build
