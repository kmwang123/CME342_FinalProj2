#CXX := g++
MPICC := mpic++
CXXFLAGS := -O3 -Wall -Wextra -Wconversion -std=c++11 #-DNDEBUG 
NPROCS := 1

TARGET := main
OBJS := main.o mesh.o Array.o ionTransportEqns.o NSEqns.o
INCS := mesh.hpp Array.hpp ionTransportEqns.hpp NSEqns.hpp

$(TARGET): $(OBJS)
	$(MPICC) -o $(TARGET) $(OBJS)

%.o: %.cpp $(INCS)
	$(MPICC) -c -o $@ $< $(CXXFLAGS)

# use .PHONY for targets that do not produce a file
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET) *~

.PHONY: run
run: $(TARGET)
	mpirun -np $(NPROCS) ./$(TARGET)
