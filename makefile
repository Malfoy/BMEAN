# CC=/usr/bin/g++
CXX= g++
CFLAGS = -O3 -std=c++11 -lpthread -Wall
EXEC=testLR
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



testLR:  testLR.cpp bmean.o utils.o
	$(CXX) -o $@  $^ $(CFLAGS)

%.o: %.cpp %.h
	$(CXX) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm testLR

rebuild: clean $(EXEC)
