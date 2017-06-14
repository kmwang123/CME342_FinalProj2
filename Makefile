#CXX := g++
MPICC := mpic++
#CXXFLAGS := -g -Wall -Wextra -Wconversion -std=c++11 #-DNDEBUG 
CXXFLAGS := -O3 -Wall -Wextra -Wconversion -std=c++11 #-DNDEBUG 
NPROCS := 1

#HYPRE_DIR := /home/kmwang14/packages/hypre-2.10.0b/src/hypre
HYPRE_DIR := /home/kmwang14/packages/hypre-2.11.2/src/hypre
BOOST_DIR := /home/kmwang14/packages/boost_1_61_0

TARGET := main
OBJS := main.o mesh.o Array.o ionTransportEqns.o NSEqns.o mpiWrapper.o hypreSolver.o hypreSolverSStruct.o
INCS := mesh.hpp Array.hpp ionTransportEqns.hpp NSEqns.hpp mpiWrapper.hpp hypreSolver.hpp hypreSolverSStruct.hpp

$(TARGET): $(OBJS)
	$(MPICC) -o $(TARGET) $(OBJS) -I$(HYPRE_DIR)/include/ -I$(BOOST_DIR) -L$(HYPRE_DIR)/lib -lHYPRE -lm -lstdc++ 

%.o: %.cpp $(INCS)
	$(MPICC) -c -o $@ $< $(CXXFLAGS) -I$(HYPRE_DIR)/include/ -I$(BOOST_DIR) -L$(HYPRE_DIR)/lib -lHYPRE -lm -lstdc++

# use .PHONY for targets that do not produce a file
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET) *~

.PHONY: run
run: $(TARGET)
	mpirun -np $(NPROCS) ./$(TARGET)
