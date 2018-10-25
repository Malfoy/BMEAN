# CC=/usr/bin/g++
CXX= g++
CFLAGS = -O3 -std=c++11 -lpthread  -IBOA
LFLAGS =-IBOA
EXEC=testLR
OBJS := *.o BOA/align_lpo2.o  BOA/align_lpo_po2.o  BOA/align_score.o  BOA/black_flag.o  BOA/buildup_lpo.o  BOA/create_seq.o  BOA/fasta_format.o  BOA/heaviest_bundle.o  BOA/lpo_format.o  BOA/lpo.o   BOA/msa_format.o  BOA/numeric_data.o  BOA/remove_bundle.o  BOA/seq_util.o  BOA/stringptr.o SSW/src/*.o

all: $(EXEC)

ifeq ($(prof),1)
 CFLAGS+= -pg
endif
ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g
endif

ifeq ($(sani),1)
 CFLAGS= -std=c++11 -lpthread -fsanitize=address -fno-omit-frame-pointer -O1 -fno-optimize-sibling-calls -g
endif



test:
	./testLR

all: $(EXEC)



testLR:  testLR.o bmean.o utils.o
	$(CXX) $(LFLAGS)  $(OBJS) -o $@

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm testLR

rebuild: clean $(EXEC)
