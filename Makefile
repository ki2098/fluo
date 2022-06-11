.alloc:
	g++ -g -c src/alloc.cpp -o obj/alloc.o

.setup:
	g++ -g -c src/setup.cpp -o obj/setup.o

.main:
	g++ -g -c src/main.cpp -o obj/main.o

.yyjson:
	g++ -g -c lib/yyjson.c -o obj/yyjson.o

all: .alloc .setup .main .yyjson
	g++ -g obj/main.o obj/alloc.o obj/setup.o obj/yyjson.o -o bin/fluo++.exe

clean:
	rm bin/*
	rm obj/*