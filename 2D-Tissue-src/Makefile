# Makefile
# Haibo Ni <haibo.ni02@gmail.com> // 15:53:54, Thu, 23-January-2020, By Haibo
# Grandi model make
CXX=mpiicpc #cxx  #mpicxx
# CXXFLAGS = -O3 -std=c++11  -qopenmp -std=c++11 -g -lz   
#if mpicxx
# change to 
CXXFLAGS = -O3 -std=c++11  -qopenmp -std=c++11 -g -lz # (if g++)


INCLUDE=-I/usr/local/include
LIB    =-L/usr/local/lib

LDFLAGS= -lz  -lsundials_cvode ${INCLUDE} ${LIB} #  -vec-report -profile-functions -profile-loops=all -profile-loops-report=2 

INC_PARAMS=-Ilib -Ilib_tissue -I/usr/local/include
# CXX = g++
srcs = $(wildcard lib_tissue/*.cpp lib/*.cpp)

BUILD_DIR=build

objs := $(srcs:%=$(BUILD_DIR)/%.o)
deps := $(objs:.o=.d)


# Atria_3D_OMP_2D_2d: Atria_3D_OMP_2D_2d.cpp $(objs)  
Atria_2D_MPI_Ghost: Atria_2D_MPI_Ghost.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@ $(LDFLAGS)

# Benchmark_1D_Diffusion:	Benchmark_1D_Diffusion.cpp  $(objs)  
# 	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@ $(LDFLAGS)

# GB_Main_Vm_as_para:	GB_Main_Vm_as_para.cpp $(objs) 
# 	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@ $(LDFLAGS)
# GB_Main_lsoda_Vm_as_para:	GB_Main_lsoda_Vm_as_para.cpp $(objs) 
# 	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@ $(LDFLAGS)


# GB_Main_Vm_as_para_test:	GB_Main_Vm_as_para_test.cpp $(objs) 
# 	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(INC_PARAMS)   $(CXXFLAGS)  -MMD -MP -c $< -o $@ 

.PHONY: clean

# $(RM) is rm -f by default
clean:
	$(RM) $(objs) $(deps) GB_Main

-include $(deps)
MKDIR_P ?= mkdir -p