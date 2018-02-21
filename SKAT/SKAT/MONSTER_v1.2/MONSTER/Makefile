#Makefile

#Change this variable, if GSL is not installed under the default library path for the compiler. Example: -L /usr/local/lib
LIBS_PATH = -L /usr/lib

LIBS_LNX_D_GSL = -lgsl -lgslcblas -lm
PATH_TO_EIGEN = ./src/Eigen

SRC_DIR  = ./src

CPP = g++

CPPFLAGS = -w -O3

OUTPUT = ./MONSTER

all: 
	$(CPP) ./src/MONSTER.cpp $(LIBS_PATH) $(LIBS_LNX_D_GSL) -I $(PATH_TO_EIGEN) $(CPPFLAGS) -o $(OUTPUT) 

debug: 
	$(CPP) $(LIBS_LNX_D_GSL) -I $(PATH_TO_EIGEN) -O0 -g ./src/MONSTER.cpp -o $(OUTPUT)

clean:
	rm -f $(OUTPUT)

static:
	$(CPP) ./src/MONSTER.cpp -static $(LIBS_PAGH) $(LIBS_LNX_D_GSL) -I $(PATH_TO_EIGEN) $(CPPFLAGS) -o ./staticMONSTER
