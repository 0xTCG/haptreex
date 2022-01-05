LLC := llc
CLANG := clang

all: haptreex

haptreex:
	SEQ_LIBRARY_PATH=~/.seq/lib/seq ~/.seq/bin/seqc build -release -o haptreex src/main.seq
	patchelf --set-rpath '$$ORIGIN' haptreex

clean:
	rm -rf build
