CC = g++
SAMLIB = samtools-0.1.18
SAMOBJLIBS = $(SAMLIB)/sam.o $(SAMLIB)/bam.o $(SAMLIB)/bgzf.o $(SAMLIB)/kstring.o $(SAMLIB)/bam_import.o $(SAMLIB)/faidx.o $(SAMLIB)/bam_pileup.o $(SAMLIB)/bam_aux.o $(SAMLIB)/sam_header.o $(SAMLIB)/razf.o
CFLAGS = -O3 -std=c++11
ZFLAGS = -lz -lpthread

all : cor

cor : fileHandler.cpp fileHandler.h main.cpp mylib.cpp mylib.h readUpdate.cpp readUpdate.h definitions.h $(SAMOBJLIBS)
	g++ -o cor fileHandler.cpp main.cpp mylib.cpp readUpdate.cpp $(SAMOBJLIBS) $(ZFLAGS) $(CFLAGS)

clean :
	rm -rf cor $(SAMOBJLIBS)
