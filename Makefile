main.out: main.cpp Protein.o AminoAcidChain.o Structure.o AminoAcid.o Sequence.o LinearProbing.o
	g++ -g -Wall -o main.out main.cpp Protein.o AminoAcidChain.o Structure.o AminoAcid.o Sequence.o LinearProbing.o

Protein.o: Protein.h Protein.cpp AminoAcidChain.h Structure.h AminoAcid.h Sequence.h LinearProbing.h
	g++ -g -Wall -c Protein.cpp

AminoAcidChain.o: AminoAcidChain.h AminoAcidChain.cpp Structure.h AminoAcid.h Sequence.h LinearProbing.h
	g++ -g -Wall -c AminoAcidChain.cpp

Structure.o: Structure.h Structure.cpp
	g++ -g -Wall -c Structure.cpp

AminoAcid.o: AminoAcid.h AminoAcid.cpp LinearProbing.h
	g++ -g -Wall -c AminoAcid.cpp

Sequence.o: Sequence.h Sequence.cpp LinearProbing.h
	g++ -g -Wall -c Sequence.cpp

LinearProbing.o: LinearProbing.h LinearProbing.cpp
	g++ -g -Wall -c LinearProbing.cpp