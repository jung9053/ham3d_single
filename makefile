                     #Makefile for ham2d

F90 = ifort
#F90 = gfortran
CC  = icc
#CC  = gcc
FFLAGS =  -w -O2 -r8 
#FFLAGS  =  -w -O2 -fdefault-real-8 # uncomment this for gfortran
CFLAGS =  -w -O2 #-g #-traceback  -check uninit  #-warn-all#-g #-O2 #-g

#FFLAGS =  -w -O2 -lgfortran -fdefault-real-8 # uncomment this for gfortran
#CFLAGS =  -w -O2 -fbacktrace -fbounds-check # -g #-traceback  -check uninit  #-w


LFLAGS =$(CFLAGS) -lm
OBJECTS = ham3d.o readGrid.o preprocess.o find_faces.o\
	  initflow.o stepSolution.o updateSoln.o computeRHS.o computeRHSv.o\
          computeRHSk.o computeRHSkv.o computeLinearRHS.o flux_roe3d.o wallFlux.o flux_roe2d.o periodic_bc.o apply_periodic.o apply_periodic_LHS.o\
	  muscld.o muscld_deriv.o roeflx.o computeForce.o outputSolution.o wrest.o\
	  smoothGrid.o jac_roe.o flux_visc.o triSolvers.o mathops.o findDiagonals.o\
	  ADI.o DADI.o jac_visc.o gaussSeidel.o lineGaussSeidel.o \
     weno.o weno_deriv.o\

# Link Instruction
ham2d: $(OBJECTS)
	$(CC) $(OBJECTS) $(LFLAGS) -o ham3d

clean:
	@rm -rf *.o *.mod *.*~ ham3d

%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o

%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.c
	$(CC) $(CFLAGS) -c $< -o $*.o

