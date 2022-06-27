.yyjson:
	g++ -c ./lib/yyjson.c -o ./lib/yyjson.o

.boundary:
	nvc++ -c ./src/boundary.cpp -acc -Minfo -o ./obj/boundary.o

.domain:
	nvc++ -c ./src/domain.cpp -acc -Minfo -o ./obj/domain.o

.fluo:
	nvc++ -c ./src/fluo.cpp -acc -Minfo -o ./obj/fluo.o

.jreader:
	g++ -c ./src/jreader.cpp -o ./obj/jreader.o

.main:
	nvc++ -c ./src/main.cpp -acc -Minfo -o ./obj/main.o

.mmac:
	nvc++ -c ./src/mmac.cpp -acc -Minfo -o ./obj/mmac.o

.poisson:
	nvc++ -c ./src/poisson.cpp -acc -Minfo -o ./obj/poisson.o

all: .yyjson .boundary .domain .fluo .jreader .main .mmac .poisson
	nvc++ ./lib/yyjson.o ./obj/boundary.o ./obj/domain.o ./obj/fluo.o ./obj/jreader.o ./obj/main.o ./obj/mmac.o ./obj/poisson.o -acc -Minfo -o ./bin/fluo++

clean:
	rm ./obj/*
	rm ./bin/*