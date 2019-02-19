cmake_minimum_required (VERSION 3.0)

project (ParallelFortran)

enable_language(Fortran)



##################################
# compiler options
##################################

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-w -O3 -cpp -dM -fimplicit-none -ffixed-line-length-none")
    set(bounds "-fbounds-check")
endif()


set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

#set(COMMON_SRC_FILES, elementutilitiesbasisfuncs.F elementutilitieselasticity2D.F elementutilitieselasticity3D.F writervtk.F solverpetsc.F)

# this is where we will place the Fortran module files
#set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)


##################################
# include directories
##################################
#include_directories(/usr/include/petsc)
include_directories(
/usr/include/vtk-6.2  /usr/include/openmpi  /usr/include/petsc  /usr/include/superlu  /usr/include/boost)


##################################
# lib directories
##################################
#link_directories(/usr/lib/gcc/x86_64-linux-gnu /usr/lib /opt/petsc-3.6.4/arch-linux2-c-debug/lib)
link_directories(
/usr/lib/gcc/x86_64-linux-gnu/5.4.0
/usr/lib/x86_64-linux-gnu
/usr/lib/libblas
/usr/lib/lapack 
/opt/petsc-3.6.4/arch-linux2-c-debug/lib 
/usr/lib 
/home/chenna/Documents/myCode/igampap/lib  
/usr/lib/x86_64-linux-gnu  
/usr/x86_64-linux-gnu/lib
)



#src/elemutilitiesquadrature.F 
#file(GLOB_RECURSE sources src/elemutilitiesquadrature.F src/elementutilitiesbasisfuncs.F src/elementutilitiespoisson.F src/elementutilitieselasticity2D.F elementutilitieselasticity3D.F src/writervtk.F src/solverpetsc.F)

##################################

add_executable(partmeshfort src/meshpartitionfort.F)
#add_executable(partmeshfortparmetis src/meshpartitionfortparmetis.f90)
add_executable(parmetisdualtest src/parmetisdualtest.F)


target_link_libraries(partmeshfort metis)
#target_link_libraries(partmeshfortparmetis mpi_mpifh parmetis)
target_link_libraries(parmetisdualtest mpi_mpifh parmetis)

install(TARGETS partmeshfort RUNTIME DESTINATION /home/chenna/Documents/myCode/Fortran/standardfem/bin)
#install(TARGETS partmeshfortparmetis RUNTIME DESTINATION /home/chenna/Documents/myCode/Fortran/standardfem/bin)
install(TARGETS parmetisdualtest RUNTIME DESTINATION /home/chenna/Documents/myCode/Fortran/standardfem/bin)

#add_library(writer writervtk.F)
#add_executable(writer writervtk.F)

#add_executable(tetrapoisson src/tetrapoissonparallelimpl1.F ${sources})
#add_executable(tetraelasticity src/tetraelasticityparallelimpl1.F ${sources})

#add_executable(triapoisson src/triapoissonparallelimpl1.F ${sources})
#add_executable(elasticityexplicittria src/triaelasticityexplicit.F ${sources})
#add_executable(elasticityexplicitquad src/quadelasticityexplicit.F ${sources})
#add_executable(elasticityexplicitprism src/prismelasticityexplicit.F ${sources})
#add_executable(elasticityexplicittria6 src/elasticityexplicittria6node.F ${sources})

#add_executable(elasticityexplicit3d src/tetraelasticityexplicit.F ${sources})

#target_link_libraries(tetrapoisson metis petsc mpi mpi_mpifh)
#target_link_libraries(tetraelasticity metis petsc mpi mpi_mpifh)

#target_link_libraries(triapoisson metis petsc mpi mpi_mpifh)
#target_link_libraries(elasticityexplicittria metis petsc mpi mpi_mpifh)
#target_link_libraries(elasticityexplicitquad metis petsc mpi mpi_mpifh)
#target_link_libraries(elasticityexplicitprism metis petsc mpi mpi_mpifh)
#target_link_libraries(elasticityexplicittria6 metis petsc mpi mpi_mpifh)

#target_link_libraries(elasticityexplicit3d metis petsc mpi mpi_mpifh)


#install(TARGETS triapoisson RUNTIME DESTINATION /home/chenna/Documents/myCode/Fortran/standardfem/bin)
#install(TARGETS elasticityexplicittria RUNTIME DESTINATION /home/chenna/Documents/myCode/Fortran/standardfem/bin)
#install(TARGETS elasticityexplicitquad RUNTIME DESTINATION /home/chenna/Documents/myCode/Fortran/standardfem/bin)
#install(TARGETS elasticityexplicitprism RUNTIME DESTINATION /home/chenna/Documents/myCode/Fortran/standardfem/bin)
#install(TARGETS elasticityexplicittria6 RUNTIME DESTINATION /home/chenna/Documents/myCode/Fortran/standardfem/bin)




