
########################################
# to run the program
# mpirun -np <nproc> <exename> <nodefile> <elemfile> <DirichletBCfile> <ForceBCfile>
#
# For example,
# to run 3D Poisson problem with 2 processors
# mpirun -np 2 ./tetrapoissonparallelimpl1 tet50-nodes.dat tet50-elems.dat tet50-DirichBC.dat

# to run 2D Elasticity problem with 4 processors
# mpirun -np 4 ./triaelasticityparallelimpl1 cookmembranetria-nodes.dat cookmembranetria-elems.dat cookmembranetria-DirichBC.dat cookmembranetria-ForceBC.dat


# --- compiler settings -----------------------------------

VPATH  =  ../src


VTKLIB62 = /usr/lib/x86_64-linux-gnu  -lvtkCommonCore-6.3  -lvtkIOCore-6.3  -lvtkIOGeometry-6.3  -lvtkIOXML-6.3 -lvtkIOImage-6.3 -lvtkIOExport-6.3 \
         -lvtkRenderingCore-6.3   -lvtkFiltersCore-6.3  -lvtkFiltersGeneral-6.3    -lvtkFiltersModeling-6.3  -lvtkFiltersSources-6.3  \
         -lvtkRenderingAnnotation-6.3  -lvtkRenderingVolume-6.3  -lvtkRenderingVolumeOpenGL-6.3 \
         -lvtkInteractionStyle-6.3   -lvtkIOParallel-6.3  -lvtkIOParallelXML-6.3  -lvtkCommonExecutionModel-6.3  -lvtkCommonDataModel-6.3 \
	 -lvtkRenderingOpenGL-6.3  -lvtkRenderingFreeType-6.3 -lvtkIOLegacy-6.3

CC = gcc 

COPT  = -std=c++11 -O3 -Wno-uninitialized -Wno-sign-compare -Wno-write-strings  -Wno-deprecated  -Wno-unused  -Wno-format -Wno-narrowing  -Wno-reorder  -Wreturn-type -fpermissive -frounding-math


CINCL = -I /usr/include/vtk-6.3 \
	-I /usr/include/petsc \
	-I /usr/include/openmpi  

CLIB  = -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0 -lstdc++ -lgfortran -fopenmp -lgomp -lpthread -lm  \
        -L$(VTKLIB62) \
		-L/usr/lib/gcc/x86_64-linux-gnu -lmetis \
		-L/usr/lib -lpetsc \
 		-L/usr/lib -lmpi

CF    = gfortran
FOPT  = -w -O3 -cpp -dM -fimplicit-none -ffixed-line-length-none

FINCL = -I /usr/include/openmpi  \
	-I /usr/include/petsc

FLIB  = -L /usr/include/petsc \
		-L/usr/lib -lmpi  -lmpi_mpifh \
		-L/usr/lib -lpetsc \
		-L/usr/lib/gcc/x86_64-linux-gnu -lmetis 

# --- obj file list

OBJECTS = writervtk.o solverpetsc.o elemutilitiesquadrature.o elementutilitiesbasisfuncs.o elementutilitiespoisson.o elementutilitieselasticity2D.o  elementutilitieselasticity3D.o

#OBJECTS=$(patsubst %.f,%.o,$(wildcard *.f))

# --- commands -------------------------------------------------------

#include makedepslist

%.o ../obj/%.o: %.cpp
	$(CC) -c $(COPT) $(CINCL) $< -o $@

%.o ../obj/%.o: %.f90
	$(CF) -c $(FOPT) $(FINCL) $< -o $@

%.o ../obj/%.o: %.F90
	$(CF) -c $(FOPT) $(FINCL) $< -o $@

%.o ../obj/%.o: %.f
	$(CF) -c $(FOPT) $(FINCL) $< -o $@

%.o ../obj/%.o: %.F
	$(CF) -c $(FOPT) $(FINCL) $< -o $@

%.o %.mod: %.F
#%.o: %.F
	$(CF) -c $(FOPT) $(FINCL) $< -o $@
	mv $@ ../obj
	mv *.mod ../obj



vtktovtu: vtktovtu.o
	${CC} -o vtktovtu vtktovtu.o ${CLIB}
	rm  vtktovtu.o

writervtk: writervtk.o
	${CF} -o writervtk writervtk.o ${FLIB}
	rm  writervtk.o

solverpetsc: solverpetsc.o
	${CC} -o solverpetsc solverpetsc.o ${FLIB}
	rm  solverpetsc.o

elementutilitiesbasisfuncs: elementutilitiesbasisfuncs.o
	${CF} -o elementutilitiesbasisfuncs elementutilitiesbasisfuncs.o ${FLIB}
	rm  elementutilitiesbasisfuncs.o

elementutilitieselasticity2D: elementutilitieselasticity2D.o
	${CF} -o elementutilitieselasticity2D elementutilitieselasticity2D.o ${FLIB}
	rm  elementutilitieselasticity2D.o

elementutilitieselasticity3D: elementutilitieselasticity3D.o
	${CF} -o elementutilitieselasticity3D elementutilitieselasticity3D.o ${FLIB}
	rm  elementutilitieselasticity3D.o

elementutilitiespoisson: elementutilitiespoisson.o
	${CF} -o elementutilitiespoisson elementutilitiespoisson.o ${FLIB}
	rm  elementutilitiespoisson.o

meshpartitioncpp: writervtk.o meshpartitioncpp.o
	${CC} -o meshpartitioncpp meshpartitioncpp.o writervtk.o ${CLIB}
	rm  meshpartitioncpp.o

meshpartitionfort: writervtk.o meshpartitionfort.o
	${CF} -o meshpartitionfort meshpartitionfort.o writervtk.o ${FLIB}
	rm  meshpartitionfort.o

triapoissonserialimpl1: writervtk.o triapoissonserialimpl1.o ${OBJECTS}
	${CF} -o triapoissonserialimpl1 triapoissonserialimpl1.o writervtk.o ${FLIB}
	rm  triapoissonserialimpl1.o

triapoissonserialimpl2: writervtk.o triapoissonserialimpl2.o
	${CF} -o triapoissonserialimpl2 triapoissonserialimpl2.o writervtk.o ${FLIB}
	rm  triapoissonserialimpl2.o

triapoissonserialimpl3: writervtk.o solverpetsc.o triapoissonserialimpl3.o
	${CF} -o triapoissonserialimpl3 triapoissonserialimpl3.o ${OBJECTS} ${FLIB}
	rm  triapoissonserialimpl3.o

triapoissonparallelimpl1: ${OBJECTS} triapoissonparallelimpl1.o 
	${CF} -o triapoissonparallelimpl1 triapoissonparallelimpl1.o ${OBJECTS} ${FLIB}
	rm  triapoissonparallelimpl1.o

triapoissonparallelimpl2: ${OBJECTS} triapoissonparallelimpl2.o 
	${CF} -o triapoissonparallelimpl2 triapoissonparallelimpl2.o ${OBJECTS} ${FLIB}
	rm  triapoissonparallelimpl2.o

tetrapoissonparallelimpl1: ${OBJECTS} tetrapoissonparallelimpl1.o 
	${CF} -o tetrapoissonparallelimpl1 tetrapoissonparallelimpl1.o ${OBJECTS} ${FLIB}
	rm  tetrapoissonparallelimpl1.o

triaelasticityparallelimpl1: ${OBJECTS} triaelasticityparallelimpl1.o
	${CF} -o triaelasticityparallelimpl1 triaelasticityparallelimpl1.o ${OBJECTS} ${FLIB}
	rm  triaelasticityparallelimpl1.o

tetraelasticityparallelimpl1: ${OBJECTS} tetraelasticityparallelimpl1.o
	${CF} -o tetraelasticityparallelimpl1 tetraelasticityparallelimpl1.o ${OBJECTS} ${FLIB}
	rm  tetraelasticityparallelimpl1.o

#meshpartitioncpp: meshpartitioncpp.o
#	${CC} -o meshpartitioncpp meshpartitioncpp.o ${CLIB}
#	rm  meshpartitioncpp.o

clean:
	rm *.o

