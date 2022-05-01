all: build $(wildcard *.dat)
	./build -a

build: $(wildcard *.c *.h)
	clang -O3 *.c -o build