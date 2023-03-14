# Software: COR (Consensus Overrides Reads)

COR (Consensus Overrides Reads) aims to reduce the sequencing errors and remove the heterozygous information on the reads in order to reduce the complexity in the subsequent assembly process. COR updates the read sequences, which should be from the same species, according to the consensus of the read alignments against a reference sequence.

## Installation

The software was written in C++, and it has been tested under linux and MacOS platform. You need
to have C++ compiler installed in the machine in order to compile the source codes. The compilation
steps are shown as follows:

```
tar -zxvf cor-x.x.tar.gz
cd cor-x.x
make
```

Then an executable file named *cor* will appear

## Usage

Syntax:

```
./cor [in: sorted SAM file] [out: sorted SAM file]
```

Please note that the input SAM file has to be sorted.
