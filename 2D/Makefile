# Global Variables controlling the compilation 
# h5fc is the fortran compiler with implied hdf5 linking

FC = h5fc
FFLAGS = -march=native -Wall -Ofast -fPIC -g #-fopenmp
F2PYFL = --f90flags=-fopenmp -lgomp 



consts.mod: consts.f95
	$(FC) $(FFLAGS) -c consts.f95

eno.mod: eno.f95
	$(FC) $(FFLAGS) -c eno.f95

fluxes.mod: fluxes.f95 consts.mod
	$(FC) $(FFLAGS) -c fluxes.f95

solvers.mod: solvers.f95 consts.mod fluxes.mod gtk.mod eno.mod eiep.mod supmoll.mod fluxmatsparam.mod
	$(FC) $(FFLAGS) -c solvers.f95  

supmoll.mod: supmoll.f95 consts.f95 consts.mod
	$(FC) $(FFLAGS) -c supmoll.f95

eiep.mod: consts.f95 fluxes.f95 eiep.f95
	$(FC) $(FFLAGS) -c eiep.f95

tc.mod: tc.f95 consts.mod
	$(FC) $(FFLAGS) -c tc.f95 

ffsdomain.mod: consts.mod ffsdomain.f95 tc.mod
	$(FC) $(FFLAGS) -c ffsdomain.f95

fluxmatsparam.mod: fluxmat6.f95
	$(FC) $(FFLAGS) -c fluxmat6.f95

gtk.mod: gtk.f95 fluxmatsparam.mod consts.mod fluxes.mod
	$(FC) $(FFLAGS) -c gtk.f95

ntffs: ffs
	ulimit -Ss unlimited
	./ffs


clean:
	rm *.o *.mod *.so *.h5

ffs: ffs.f95 consts.mod fluxes.mod solvers.mod eno.mod supmoll.mod eiep.mod ffsdomain.mod gtk.mod fluxmatsparam.mod tc.mod
	$(FC) $(FFLAGS) -c ffs.f95
	$(FC) $(FFLAGS) consts.o fluxes.o solvers.o eno.o supmoll.o eiep.o ffsdomain.o gtk.o fluxmat6.o ffs.o tc.o -o ffs


solvers: ffs
	f2py3  $(F2PYFL) -c  consts.f95 eno.f95 tc.f95 supmoll.f95 fluxes.f95 solvers.f95 gtk.f95  ffsdomain.f95 fluxmat6.f95 eiep.f95 -m solvers
