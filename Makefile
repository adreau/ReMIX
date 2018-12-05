all: bin/filterFastq_by_bc bin/constructHbopFiles_harMean bin/XO_fromHBOP_withh0h1parity

bin/filterFastq_by_bc:
	g++ -std=c++11 -Wall -o bin/filterFastq_by_bc src/filterFastq_by_bc.cpp

bin/constructHbopFiles_harMean:
	g++ -std=c++11 -Wall -o bin/constructHbopFiles_harMean src/constructHbopFiles_fromRepMol_harMean.cpp src/Molecule.cpp

bin/XO_fromHBOP_withh0h1parity:
	g++ -std=c++11 -Wall -o bin/XO_fromHBOP_withh0h1parity src/XO_fromHBOP_withh0h1parity.cpp
	
clean:
	rm -f bin/*
