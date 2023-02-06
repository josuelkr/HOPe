# makefile for baseflow
OBJ1 = main_par_2D.o derivs_par_2D.o loop_par_2D.o penta_par_2D.o poisson_par_2D.o residual_2D.o
OBJ2 = integ_2D.o
OBJ3 = format_base.o

base : $(OBJ1)
	mpiifort $(OBJ1) -o basef

int : $(OBJ2)
	ifort $(OBJ2) -o integ

form : $(OBJ3)
	ifort $(OBJ3) -o form

.f.o:
	mpiifort -c -O2 $<
#	mpiifort -O0 -g -Wall -fbounds-check -c $<

clean:
	rm *~ *.o *.bin *.dat basef form integ saida_base tecplot.phy
