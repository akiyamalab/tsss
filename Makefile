MODE:=Release
#Release or Debug

CXX=g++
CXX_FLAGS=-Wall -std=c++11 -fopenmp

ifeq ($(MODE),Release)
	CXX_FLAGS += -O3
else ifeq ($(MODE),Debug)
	CXX_FLAGS += -g -DDEBUG -O3 #-O0 本来debugなら最適化レベルを下げる
endif

C_SRC=./ext/karlin/src/karlin.c \

CPP_SRC=./src/aligner.cpp \
	./src/alignment.cpp \
	./src/alphabet_coder.cpp \
	./src/alphabet_type.cpp \
	./src/chain_filter.cpp \
	./src/database.cpp \
	./src/database_build.cpp \
	./src/database_chunk.cpp \
	./src/dna_sequence.cpp \
	./src/dna_type.cpp \
	./src/edit_blocks.cpp \
	./src/fasta_sequence_reader.cpp \
	./src/gapped_extender.cpp \
	./src/main.cpp \
	./src/protein_query.cpp \
	./src/protein_sequence.cpp \
	./src/protein_type.cpp \
	./src/queries.cpp \
	./src/query.cpp \
	./src/reduced_alphabet_coder.cpp \
	./src/reduced_amino.cpp \
	./src/score_matrix.cpp \
	./src/score_matrix_reader.cpp \
	./src/seed_searcher.cpp \
	./src/sequence.cpp \
	./src/sequence_no_filter.cpp \
	./src/sequence_codon_filter.cpp \
	./src/statistics.cpp \
	./src/translated_dna_query.cpp \
	./src/translator.cpp \
	./src/ungapped_extender.cpp \
	./src/utils.cpp \

OBJS =
OBJS += $(C_SRC:%.c=%.o)
OBJS += $(CPP_SRC:./src/%.cpp=objs/%.o)

.PHONY: all
all:dir tsss

tsss: $(OBJS)
	$(CXX) $(CXX_FLAGS) -o $@ $(OBJS) -lboost_program_options
objs/%.o: src/%.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@

dir:
	mkdir -p objs

.PHONY: clean
clean:
	rm -rf objs tsss
