all: build $(wildcard *.dat)
	./build -a

build: $(wildcard *.c *.h)
	g++-11 -O3 *.c -o build

clean:
	rm -rf build result.dat