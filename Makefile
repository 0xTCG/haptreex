LLC := llc
CLANG := clang

all: haptreex

haptreex:
	seqc build -release -o build/haptreex src/main.seq

clean:
	rm -rf build
