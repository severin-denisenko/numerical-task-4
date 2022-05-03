all: build $(wildcard *.dat)
	./build -o

build: $(wildcard *.c *.h)
	clang -O3 *.c -o build
