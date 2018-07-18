# Only works on 64-bit.
# Edit these variables.
MEXSUFFIX   = mexa64
MEX         = mex
CXX         = g++
#PARDISOHOME = /home/simonsch/homework/pardiso_example/lib ##$(HOME)/pardiso/lib
#PARDISOHOME = ../lib
PARDISOHOME = /project2/osmlab/pardiso_example/pardiso-matlab 

# Do not edit anything past here (unless you know what you're doing).
MEXFLAGS   = -cxx CXX=$(CXX) CC=$(CXX) LD=$(CXX) -largeArrayDims \
             -L$(PARDISOHOME) -lpardiso500-GNU481-X86-64 -lmwlapack -lmwblas \
              -lgfortran -lpthread -lm 

TARGETS_WO_MEXSUFFIX = pardisoinit pardisoreorder \
          pardisofactor pardisosolve \
          pardisofree

TARGETS = pardisoinit.$(MEXSUFFIX) pardisoreorder.$(MEXSUFFIX) \
          pardisofactor.$(MEXSUFFIX) pardisosolve.$(MEXSUFFIX) \
          pardisofree.$(MEXSUFFIX)

all: $(TARGETS_WO_MEXSUFFIX)

%: %.cpp common.cpp matlabmatrix.cpp sparsematrix.cpp \
                pardisoinfo.cpp
	$(MEX) $(MEXFLAGS) -output $@ $^

clean:
	rm -f $(TARGETS) 
