
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

VTKLIB62 = /usr/lib/x86_64-linux-gnu  -lvtkCommonCore-6.2  -lvtkIOCore-6.2  -lvtkIOGeometry-6.2  -lvtkIOXML-6.2 -lvtkIOImage-6.2 -lvtkIOExport-6.2 \
         -lvtkRenderingCore-6.2   -lvtkFiltersCore-6.2  -lvtkFiltersGeneral-6.2    -lvtkFiltersModeling-6.2  -lvtkFiltersSources-6.2  \
         -lvtkRenderingAnnotation-6.2  -lvtkRenderingVolume-6.2  -lvtkRenderingVolumeOpenGL-6.2  -lvtkRenderingFreeTypeOpenGL-6.2 \
         -lvtkInteractionStyle-6.2   -lvtkIOParallel-6.2  -lvtkIOParallelXML-6.2  -lvtkCommonExecutionModel-6.2  -lvtkCommonDataModel-6.2 \
		 -lvtkRenderingOpenGL-6.2  -lvtkRenderingFreeType-6.2 -lvtkIOLegacy-6.2

CC = gcc 

COPT  = -std=c++11 -O3 -Wno-uninitialized -Wno-sign-compare -Wno-write-strings  -Wno-deprecated  -Wno-unused  -Wno-format -Wno-narrowing  -Wno-reorder  -Wreturn-type -fpermissive -frounding-math


CINCL = -I /usr/include/vtk-6.2 \
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

OBJECTS = writervtk.o solverpetsc.o elementutilitiesbasisfuncs.o elementutilitiespoisson.o elementutilitieselasticity.o


# --- commands -------------------------------------------------------

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


vtktovtu: vtktovtu.o
	${CC} -o vtktovtu vtktovtu.o ${CLIB}
	rm  vtktovtu.o

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

clean:
	rm *.o

