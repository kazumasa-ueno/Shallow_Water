#environment
FC := gfortran

#options
FFLAGS := -O3
DEBUG_FLAGS := -g

#source(s)
SRC := constant.f90 boundary.f90 transfer.f90 calc_variables.f90 initialize.f90 solver.f90 main.f90

#object(s)
OBJ := $(patsubst %.f90, %.o, $(notdir $(SRC)))

#executable(s)
EXEC := main

.PHONY: all
all: $(EXEC)

.PHONY: debug
debug: FFLAGS += $(DEBUG_FLAGS)
debug: $(EXEC)

#rule(s)
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

# clean up
.PHONY: clean
clean:
	rm -rf $(OBJ)
	rm -rf *.mod
