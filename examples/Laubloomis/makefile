CC       	  = g++
CXXFLAGS 		= -Wall -Wextra -O3

# Please edit --->
GINACINC = ../../../../Downloads/cplusplus_packages/ginac-1.7.4/ginac
EIGENINC = ../../../../Downloads/cplusplus_packages/eigen-eigen-5a0156e40feb
LIBSINCLUDE = C:/cygwin/lib
#-----------------<

ROOT		= ../..

INC		= -I$(ROOT)/src -I$(ROOT)/external -I$(EIGENINC) -I$(GINACINC)

LIBINC = -L$(GINACINC)/.libs -L$(GINACINC) -L$(LIBSINCLUDE)

LIBS = -lboost_filesystem -lboost_iostreams -lboost_system -lginac -lcln

TARGET = laubloomis 

.PHONY: $(TARGET) 

all: $(TARGET) run

clear: 
	clear
	clear

call: clear all

laubloomis.o: laubloomis.cpp $(ROOT)/src/ReachableSet2cpp.h
	$(CC) $(INC) $(CXXFLAGS) laubloomis.cpp -c
	
laubloomis: laubloomis.o laubloomis.cpp
	$(CC) $(LIBINC) -o $(TARGET) ./laubloomis.o $(LIBS) 

run: 
	./$(TARGET)

clean:
	rm $(TARGET) datastorage.txt func_from_file.cpp func_from_file.so laubloomis.o
