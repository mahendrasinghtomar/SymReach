CC       	  = g++
CXXFLAGS 		= -Wall -Wextra -O3

# Please edit --->
GINACINC = ../../../../Downloads/cplusplus_packages/ginac-1.7.4/ginac
EIGENINC = ../../../../Downloads/cplusplus_packages/eigen-eigen-5a0156e40feb
LIBSINCLUDE = C:/cygwin/lib
#-----------------

ROOT		= ../..

INC		= -I$(ROOT)/src -I$(ROOT)/external -I$(EIGENINC) -I$(GINACINC)

LIBINC = -L$(GINACINC)/.libs -L$(GINACINC) -L$(LIBSINCLUDE)

LIBS = -lboost_filesystem -lboost_iostreams -lboost_system -lginac -lcln

.PHONY: vanderpol 

TARGET = vanderpol 

all: $(TARGET) run

vanderpol.o: vanderpol.cpp
	$(CC) $(INC) $(CXXFLAGS) vanderpol.cpp -c
	
vanderpol: vanderpol.o vanderpol.cpp
	$(CC) $(LIBINC) -o vanderpol ./vanderpol.o $(LIBS) 

run: 
	./vanderpol

clean:
	rm vanderpol 