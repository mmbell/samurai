.SUFFIXES: .cu .cu.cpp

CUDALIBPATH =  /usr/local/cuda/lib
LIBCUBLAS = $(CUDALIBPATH)/libcudart.dylib -L/Developer/CUDA/lib -lcutil
LIBCUBLASEMU = $(CUDALIBPATH)/libcudart.dylib -L/Developer/CUDA/lib -lcutil

# uncomment to run on emulator instead of the device
#DEVICEEMU       = -DDEVICEEMU
#DEVICEEMU_NVCC  = -deviceemu $(DEVICEEMU) -D_DEBUG
#LIBCUBLAS = $(LIBCUBLASEMU)

# uncomment to use faster (less precise) math functions
FASTMATH = --use_fast_math

NVCC    = /usr/local/cuda/bin/nvcc -DCUDA -ccbin /usr/bin --ptxas-options=-v
PHASE   = -ptx
PHASE   = -cuda
NVOPT   = $(DEVICEEMU_NVCC) $(PROMOTE) $(DEBUGDEBUG) $(DEBUGOUTPUT) $(FASTMATH) \
          --host-compilation 'C++' -I/Developer/CUDA/common/inc

SRCS = GPUdriver.cu.cpp

all : $(SRCS)

.cu.cu.cpp : 
	$(NVCC) $(PHASE) $(NVOPT) $<

clean::
	rm -f $(SRCS)


