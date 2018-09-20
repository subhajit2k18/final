# makefile for elasticity_2d problem
#-lineinfo ia used to facilitate profiling in nvidia visual profilor
CXX	 :=g++
CXXFLAGS :=-pg
NVCC 	 :=nvcc
NVCCFLAGS :=

EXE	 :=topOp
#DIRECTORY1 :=./CSR_SRC
#DIRECTORY1 :=./No_format
#DIRECTORY2 :=./common_SRC

CUDA_SRC :=main.cu assembly_parallel.cu boundary.cu brick_element.cu Gausspoints_wghts.cu topOp.cu
CPP_SRC1 := mesh.cpp boundary.cpp pre-processing.cpp data_read.cpp brick_element.cpp main.cpp 
#CPP_SRC2 :=brick_element.cpp Gausspoints_wghts.cpp mesh.cpp cgradient_v2.cpp
CUDA_H   :=header.h
CPP_H	 :=header_cpp.h

CUDA_OBJ :=$(patsubst %.cu,%.cu.o,$(CUDA_SRC))
CPP_OBJS1 :=$(patsubst %.cpp,%.cpp.o,$(CPP_SRC1))
#CPP_OBJS2 :=$(patsubst %.cpp,$(DIRECTORY2)/%.cpp.o,$(CPP_SRC2))


$(EXE):$(CPP_OBJS1) $(CUDA_OBJ)
	$(NVCC) -o $@ $^ -arch=sm_30 -pg

%.cpp.o:%.cpp $(CPP_H)
	$(CXX) -c -o $@ $< -O2 -g
%.cu.o:%.cu $(CUDA_H)
	$(NVCC) -dc -o $@ $< -lineinfo -pg -arch=sm_30 -I /home/pareto/cuda_lib/ #-maxrregcount 32

.PHONY:clean

clean:
	rm $(EXE) $(CPP_OBJS1) $(CUDA_OBJ)
