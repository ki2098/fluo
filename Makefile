.yyjson:
	g++ -c ./lib/yyjson.c -o ./lib/yyjson.o

.jreader:
	g++ -c ./src/jreader.cpp -o ./obj/jreader.o

.alloc:
	g++ -c ./src/alloc.cpp -o ./obj/alloc.o

.boundary:
	g++ -c ./src/boundary.cpp -o ./obj/boundary.o

.domain:
	g++ -c ./src/domain.cpp -o ./obj/domain.o

.util:
	g++ -c ./src/util.cpp -o ./obj/util.o

.main:
	g++ -c ./src/main.cpp -o ./obj/main.o

.driver:
	g++ -c ./src/driver.cpp -o ./obj/driver.o

all: .yyjson .jreader .alloc .boundary .domain .util .main .driver
	g++ obj/util.o obj/jreader.o obj/main.o obj/b.o ./obj/domain.o ./obj/alloc.o ./obj/driver.o lib/yyjson.o -o bin/fluo++