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

.mmac_poisson:
	nvc++ -c ./src/mmac_poisson.cpp -acc -Minfo -o ./obj/mmac_poisson.o

.turbulence:
	nvc++ -c ./src/turbulence.cpp -acc -Minfo -o ./obj/turbulence.o

all: .yyjson .boundary .domain .fluo .jreader .main .mmac .mmac_poisson .turbulence
	nvc++ ./lib/yyjson.o ./obj/boundary.o ./obj/domain.o ./obj/fluo.o ./obj/jreader.o ./obj/main.o ./obj/mmac.o ./obj/mmac_poisson.o ./obj/turbulence.o -acc -Minfo -o ./bin/fluo++

clean:
	rm ./obj/*
	rm ./bin/*