.SUFFIXES: .cu .cu.cpp

CUDALIBPATH =  /usr/local/cuda/lib
LIBCUBLAS = $(CUDALIBPATH)/libcudart.dylib -L/Developer/CUDA/lib -lcutil
LIBCUBLASEMU = $(CUDALIBPATH)/libcudart.dylib -L/Developer/CUDA/lib -lcutil

# uncomment to run on emulator instead of the device
#DEVICEEMU       = -DDEVICEEMU
#DEVICEEMU_NVCC  = -deviceemu $(DEVICEEMU) -D_DEBUG
#LIBCUBLAS = $(LIBCUBLASEMU)

NVCC    = /usr/local/cuda/bin/nvcc -DCUDA -ccbin /usr/bin --ptxas-options=-v
PHASE   = -ptx
PHASE   = -cuda
NVOPT   = $(DEVICEEMU_NVCC) $(PROMOTE) $(DEBUGDEBUG) $(DEBUGOUTPUT) \
          --host-compilation 'C++' --use_fast_math -I/Developer/CUDA/common/inc

SRCS = GPUdriver.cu.cpp

all : $(SRCS)

.cu.cu.cpp : 
	$(NVCC) $(PHASE) $(NVOPT) $<

clean::
	rm -f $(SRCS)


