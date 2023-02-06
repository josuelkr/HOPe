# makefile for gvts
OBJ1 = main_par.o derivs_par.o loop_par.o nlterms_par.o escreve_par.o penta_par.o poisson_par.o filter_par.o fft.o
OBJ2 = ftanalysis_par.o fft.o
OBJ3 = isoq_par.o fft.o
OBJ4 = formatted_par.o fft.o
OBJ5 = integ.o fft.o

dns : $(OBJ1)
	mpiifort $(OBJ1) -o dns

ft : $(OBJ2)
	ifort $(OBJ2) -o ft

iso : $(OBJ3)
	mpiifort $(OBJ3) -o isoq

form : $(OBJ4)
	ifort $(OBJ4) -o form

int : $(OBJ5)
	ifort $(OBJ5) -o int

.f.o:
	mpiifort -c -O3 -mcmodel=medium $<
#	mpif77 -O0 -g -Wall -fbounds-check -c $<

clean:
	rm *~ *.o *.dat *.plt fs prog ft isoq form int saida pert_* tecplot.phy
