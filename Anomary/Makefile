#environment
FC			:= gfortran

#options
FFLAGS	:= 

#source(s)
SRC			:= boundary.f90 calc_variables.f90 transfer.f90 main.f90
# SRC			:= boundary.f90 translation.f90 test.f90

#object(s)
OBJ			:= $(patsubst %.f90, %.o, $(notdir $(SRC)))

#executable(s)
EXEC		:= main
# EXEC		:= test
.PHONY: all
all:		$(EXEC)

#rule(s)
$(EXEC):	$(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
%.o:		%.f90
	$(FC) $(FFLAGS) -o $@ -c $<

# clean up
.PHONY: clean
clean:
	rm -rf $(OBJ)
	rm -rf *.mod
