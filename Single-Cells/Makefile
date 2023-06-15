# Makefile
# Haibo Ni <haibo.ni02@gmail.com> // 15:53:54, Thu, 23-January-2020, By Haibo
# Grandi model make
CXX=mpiicpc

CXXFLAGS = -O3 -std=c++11  -qopenmp -std=c++11 -g -lz

#CFLAGS = -g -Ilib -ipo -O3 -openmp
INCLUDE=-I/usr/local/include
LIB    =-L/usr/local/lib

LDFLAGS= -lz  -lsundials_cvode ${INCLUDE} ${LIB} #  -vec-report -profile-functions -profile-loops=all -profile-loops-report=2 
# LDLIBS= -lz  -lsundials_cvode#-lopenblas -lpthread
# MY_APP=gfb

# INC=$(lib_cell) $(lib)
# INC_PARAMS=-Ilib -I/usr/local/include
# CXX = g++
# INC=$(lib_cell) $(lib)
INC_PARAMS=-Ilib -Ilib_tissue -I/usr/local/include
# CXX = g++
srcs = $(wildcard lib_tissue/*.cpp lib/*.cpp)
# 	   $(wildcard lib/*.cpp) 
BUILD_DIR=build

objs := $(srcs:%=$(BUILD_DIR)/%.o)
deps := $(objs:.o=.d)


# # Atria_3D_OMP_2D_2d: Atria_3D_OMP_2D_2d.cpp $(objs)  
# Atria_2D_MPI_Ghost: Atria_2D_MPI_Ghost.cpp $(objs)  
# 	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@ $(LDFLAGS)

# # Atria_3D_OMP_2D_2d: Atria_3D_OMP_2D_2d.cpp $(objs)  
# Atria_2D_MPI_Ghost_for_0DCell: Atria_2D_MPI_Ghost_for_0DCell.cpp $(objs)  
# 	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@ $(LDFLAGS)

# # Atria_3D_OMP_2D_2d: Atria_3D_OMP_2D_2d.cpp $(objs)  
# Two_D_diffusion_test: Two_D_diffusion_test.cpp $(objs)  
# 	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@ $(LDFLAGS)

# # Atria_3D_OMP_2D_2d: Atria_3D_OMP_2D_2d.cpp $(objs)  
# Three_D_diffusion_test: Three_D_diffusion_test.cpp $(objs)  
# 	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@ $(LDFLAGS)

GB_Main:	GB_Main.cpp  $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@ $(LDFLAGS)
# GB_Main_VC:	GB_Main_VC.cpp  $(objs)  
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