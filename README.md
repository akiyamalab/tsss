# TSSS
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

TSSS is a homology search tool based on two-step seed search strategy.

## Requirements
gcc version 4.3 or more  
boost version 1.32.0 or more

## Installation
Clone this repository and run make command.

```sh
$ git clone https://github.com/akiyamalab/tsss.git
$ cd tsss
$ make -j
```

## Usage

TSSS requires specifically formatted database files for homology search. These files can be generated from FASTA formatted protein sequence files. 
- Users have to prepare a database file in FASTA format and convert it into TSSS format database files by using `tsss db` command at first.
  - `tsss db` command requires 2 args ([-i database file] and [-o output file]).
  - `tsss db` command divides a database FASTA file into several database chunks and generates several files.
  - All generated files are needed for the search. Users can specify the size of each chunk. Smaller chunk size requires smaller memory, but efficiency of the search will decrease. 
- For executing homology search, `tsss aln` command is used and that command requires at least 3 args(`[-i query file]`, `[-d database file]` and `[-o output file]`).

## Example
```
$ ./tsss db  -i database.fasta -o db
$ ./tsss aln -i query.fasta -d db -o output
```
## Command and Options
`tsss db`: convert a protein fasta file to a TSSS format database file
```
usage: tsss db [options]

(Required):
  -i [ --database ] arg database file
  -o [ --output ] arg   output file

(Optional):
  -c [ --chunk_size ] arg (=1073741824) chunk size (bytes) [1073741824 (=1GB)]
  -t [ --threads ] arg (=1)             num of threads [1]
  -l [ --seed1_length ] arg (=3)        length of first seed [3]
  -a [ --seed1_amino ] arg (=14)        type of amino in first seed [14]
  -L [ --seed2_length ] arg (=5)        length of second seed [5]
  -A [ --seed2_amino ] arg (=8)         type of amino in second seed [8]
```

`tsss aln`: search homologues of queries from database
```
usage: tsss aln [options]

(Required):
  -i [ --query ] arg    query file
  -d [ --database ] arg database file
  -o [ --output ] arg   output file

(Optional):
  -q [ --query_type ] arg (=p)         query sequence type, p(protein) or d(dna) [p]
  -f [ --query_filter ] arg (=1)       filter query sequence, 1(enable) or 0(disable) [1]
  -n [ --alignments ] arg (=10)        number of outputs for each query [10]
  -t [ --threads ] arg (=1)            number of threads [1]]
  -h [ --seed1_hamming ] arg (=0)      allowed hamming distance of first seed in seed search [0]
  -H [ --seed2_hamming ] arg (=1)      allowed hamming distance of second seed in seed search [1]
  -F [ --outfmt ] arg (=1)             output format, 0(pairwise) or 1(tabular)[1]
  -e [ --evalue ] arg (=10)            evalue threshold [10]
  -c [ --chunk_size ] arg (=134217728) query chunk size [128MB]
```

## Search results
There are two output format.

### pairwise format
Pairwise format (`-F 0 or --outfmt 0`) shows alignment results like BLAST pairwise format.

### tabular format
Tabular format (`-F 1 or --outfmt 1`) shows alignment results like BLAST tabular format.
```
query database  100 20  0 0 1 20  1 20  1.00e-5 50.00
```
Each column shows:
1.  Name of a query sequence
2.  Name of a homologue sequence (subject)
3.  Sequence Identity
4.  Alignment length
5.  The number of mismatches in the alignment
6.  The number of gap openingsin the alignemt
7.  Start position of the query in the alignment
8.  End position of the query in the alignemnt
9.  Start position of the subject in the alignment
10. End position of the subject in the alignment
11. E-value
12. Normalized score

## Reference
Takabatake K, Izawa K, Akikawa M, Yanagisawa K, Ohue M, Akiyama Y. [**Improved large-scale homology search by two-step seed search using multiple reduced amino acid alphabets**](https://doi.org/10.3390/genes12091455). _Genes_, 12(9): 1455, 2021. doi:10.3390/genes12091455
