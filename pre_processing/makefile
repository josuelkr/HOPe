
OBJ1 = constants.f90 solvers.f90 derivative.f90 at_wall.f90 coef_generation.f90 

OBJ2 = constants.f90 beta_gen.f90

coef : $(OBJ1)
	ifort $(OBJ1) -o coef

beta : $(OBJ2)
	ifort $(OBJ2) -o beta

clean :
	rm *.bin *.mod *.o *~ coef fs *.dat form integ beta tecplot.phy

fs:	
	./fs.sh

.f.o:
	ifort -c -O2 $<
