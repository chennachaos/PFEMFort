! Module for PETSc based solver
! This is a wrapper for using PETSc library for creating sparse matrices
! as well as PETSc's in-built iterative solvers
!   - PCG, BiCGStab, GMRES etc.
! or third-party direct solvers
!   - MUMPS, SUPERLU, Spooles, PARDISO etc.
! or third-party preconditioners
!   - HYPRE

!  The following include statements are required for Fortran programs
!     petscsys.h       - base PETSc routines
!     petscvec.h    - vectors
!     petscmat.h    - matrices
!     petscpc.h     - preconditioners
!     petscksp.h    - Krylov subspace methods
!  Additional include statements may be needed if using additional
!  PETSc routines in a Fortran program, e.g.,
!     petscviewer.h - viewers
!     petscis.h     - index sets

! These .h files are not the same as the ones in the standard
! "include" directory of petsc installation.
! They are located in 'finclude' directory, specific for Fortran language


      MODULE Module_SolverPetsc
      ! PRIVATE

      ! ENUM, BIND(C) !:: SOLVER_STATUS
      !   ENUMERATOR :: SOLVER_EMPTY=1
      !   ENUMERATOR :: PATTERN_OK
      !   ENUMERATOR :: INIT_OK
      !   ENUMERATOR :: ASSEMBLY_OK
      !   ENUMERATOR :: FACTORISE_OK
      ! END ENUM !SOLVER_STATUS

#define PETSC_USE_FORTRAN_MODULES

#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petsckspdef.h>
#include <petsc/finclude/petscpcdef.h>

#if defined(PETSC_USE_FORTRAN_MODULES)
      USE petscvec
      USE petscmat
      USE petscksp
      USE petscpc
#endif
      IMPLICIT NONE

#if !defined(PETSC_USE_FORTRAN_MODULES)
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
!#include "petsc/finclude/petscis.h"
!#include "petsc/finclude/petscis.h90"
#endif

      INTEGER, PARAMETER :: SOLVER_EMPTY=1
      INTEGER, PARAMETER :: PATTERN_OK=2
      INTEGER, PARAMETER :: INIT_OK=3
      INTEGER, PARAMETER :: ASSEMBLY_OK=4
      INTEGER, PARAMETER :: FACTORISE_OK=5



      TYPE PetscSolver

      Vec  rhsVec     ! PETSc vector for the righthand side
      Vec  solnVec       ! PETSc vector for the solution
      Vec  solnPrev   ! PETSc vector for the solution at previous timestep
      Mat  mtx        ! PETSc matrix object for the system matrix
      KSP  ksp        ! PETSc linear solver context
      PC   pc         ! PETSc preconditioner context

      PetscInt nRow   ! Number of rows in the matrix
      PetscInt nCol   ! Number of columns in the matrix
      PetscInt nnz    ! Number of non-zero entries in the matrix

      PetscReal norm          ! Variable to compute and store the norm
    
      !PetscViewer  viewer_matx  ! Viewer for PETSc matrix
      !PetscViewer  viewer_vect  ! Viewer for PETSc vector

      PetscInt currentStatus  ! Variable to check the status of the solver

      CONTAINS

      PROCEDURE :: initialise
      PROCEDURE :: setZero
      PROCEDURE :: free
      PROCEDURE :: printInfo
      PROCEDURE :: assembleMatrix
      PROCEDURE :: assembleVector
      PROCEDURE :: assembleMatrixAndVector
      PROCEDURE :: factorise
      PROCEDURE :: solve
      PROCEDURE :: factoriseAndSolve

      END TYPE PetscSolver

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! \brief Initialise Petsc solver
      ! \param size_local  - number of local rows/columns in the matrix
      ! \param size_global - number of global rows/columns in the matrix
      ! \param diag_nnz - number of nonzeros in the diagonal matrix
      ! \param offdiag_nnz - number of nonzeros in the off-diagonal matrix
      SUBROUTINE initialise(this,
     1  size_local, size_global,
     2  diag_nnz, offdiag_nnz)
      IMPLICIT NONE
      class(PetscSolver) :: this

      INTEGER, INTENT(IN) :: size_global, size_local
      INTEGER, DIMENSION(:) :: diag_nnz, offdiag_nnz
      INTEGER :: dummy=50
      PetscErrorCode errpetsc     ! Variable to check the error message during runtime

      this%nRow = size_global
      this%nCol = size_global

!       call PetscPrintf(PETSC_COMM_WORLD,
!      1  " Creating PETSc vectors \n", errpetsc)
      !Create PETSc vector
      call VecCreate(PETSC_COMM_WORLD, this%solnVec, errpetsc)
      CHKERRQ(errpetsc)
      call VecSetSizes(this%solnVec, size_local, size_global, errpetsc)
      CHKERRQ(errpetsc)
      call VecSetFromOptions(this%solnVec, errpetsc)
      CHKERRQ(errpetsc)
      call VecDuplicate(this%solnVec, this%rhsVec, errpetsc)
      CHKERRQ(errpetsc)

      call VecSetOption(this%rhsVec,
     1  VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, errpetsc)
      
!       call PetscPrintf(PETSC_COMM_WORLD,
!      1  " Creating PETSc matrices \n", errpetsc)
      !Create PETSc matrix
      call MatCreate(PETSC_COMM_WORLD, this%mtx, errpetsc)
      CHKERRQ(errpetsc)

      call MatSetSizes(this%mtx, size_local, size_local,
     1  size_global, size_global, errpetsc)
      CHKERRQ(errpetsc)

      call MatSetFromOptions(this%mtx, errpetsc)
      CHKERRQ(errpetsc)

      call MatMPIAIJSetPreallocation(this%mtx,
     1  dummy, diag_nnz,
     2  dummy, offdiag_nnz, errpetsc)
      CHKERRQ(errpetsc)

      call MatSeqAIJSetPreallocation(this%mtx, 
     1  dummy, diag_nnz, errpetsc)
      CHKERRQ(errpetsc)

      call MatSetOption(this%mtx, MAT_NEW_NONZERO_ALLOCATION_ERR,
     1  PETSC_FALSE, errpetsc)
      CHKERRQ(errpetsc)

      call MatSetOption(this%mtx, MAT_NEW_NONZERO_LOCATIONS,
     1  PETSC_TRUE, errpetsc)
      CHKERRQ(errpetsc)

      call MatSetOption(this%mtx, MAT_KEEP_NONZERO_PATTERN,
     1  PETSC_TRUE, errpetsc)
      CHKERRQ(errpetsc)

      ! write(*,*) "Creating KSP context"
      ! Create the KSP context
      call KSPCreate(PETSC_COMM_WORLD, this%ksp, errpetsc)
      CHKERRQ(errpetsc)
      ! Set the operators for the KSP context
      call KSPSetOperators(this%ksp, this%mtx, this%mtx, errpetsc)
      CHKERRQ(errpetsc)
      ! Set the KSP type
      call KSPSetType(this%ksp, KSPCG, errpetsc)
      ! call KSPSetType(this%ksp, KSPBCGS, errpetsc)
      CHKERRQ(errpetsc)

      ! ! Set whether to use non-zero initial guess or not
      ! !KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
      ! !KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);

      ! Set KSP options from the input file
      ! This is convenient as it allows to choose different options
      ! from the input files instead of recompiling the code
      call KSPSetFromOptions(this%ksp, errpetsc)
      CHKERRQ(errpetsc)

      ! Get the PC context
      call KSPGetPC(this%ksp, this%pc, errpetsc)
      CHKERRQ(errpetsc)
      ! Set the PC context
      ! call PCSetType(this%pc, PCILU, errpetsc)
      call PCSetType(this%pc, PCBJACOBI, errpetsc)
      CHKERRQ(errpetsc)
      ! Set PC options from the input file
      call PCSetFromOptions(this%pc, errpetsc)
      CHKERRQ(errpetsc)

      this%currentStatus = SOLVER_EMPTY

      END SUBROUTINE initialise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! \brief Make all non-zero matrix entries to zero
      ! Sets all non-zero matrix and vector components to zero
      ! before starting element/vector assembly
      SUBROUTINE setZero(this)
      IMPLICIT NONE
      class(PetscSolver) :: this

      PetscErrorCode errpetsc     ! Variable to check the error message during runtime

      call MatAssemblyBegin(this%mtx, MAT_FINAL_ASSEMBLY, errpetsc)
      CHKERRQ(errpetsc)
      call MatAssemblyEnd(this%mtx,MAT_FINAL_ASSEMBLY, errpetsc)
      CHKERRQ(errpetsc)
      call MatZeroEntries(this%mtx, errpetsc)
      CHKERRQ(errpetsc)

      call VecAssemblyBegin(this%rhsVec, errpetsc)
      CHKERRQ(errpetsc)
      call VecAssemblyEnd(this%rhsVec, errpetsc)
      CHKERRQ(errpetsc)
      call VecZeroEntries(this%rhsVec, errpetsc)
      CHKERRQ(errpetsc)

      ! call VecAssemblyBegin(reacVec, errpetsc)
      ! call VecAssemblyEnd(reacVec, errpetsc);
      ! call VecZeroEntries(reacVec, errpetsc)

      END SUBROUTINE setZero


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! \brief Make all non-zero matrix entries to zero
      ! Sets all non-zero matrix and vector components to zero
      ! before starting element/vector assembly
      SUBROUTINE free(this)
      IMPLICIT NONE
      class(PetscSolver) :: this

      PetscErrorCode errpetsc     ! Variable to check the error message during runtime

      call VecDestroy(this%solnVec, errpetsc)
      CHKERRQ(errpetsc)
      ! call  VecDestroy(this%solnPrev, errpetsc)
      CHKERRQ(errpetsc)
      call  VecDestroy(this%rhsVec, errpetsc)
      CHKERRQ(errpetsc)
      call  MatDestroy(this%mtx, errpetsc)
      CHKERRQ(errpetsc)
      !ierr = KSPGetPC(ksp,&pc);
      !ierr = PCDestroy(&pc);
      ! call  PCReset(this%pc, errpetsc)
      call  KSPDestroy(this%ksp, errpetsc)
      ! call  KSPReset(this%ksp, errpetsc)
      CHKERRQ(errpetsc)

      !call PetscViewerDestroy(this%viewer_matx,errpetsc)
      !call PetscViewerDestroy(this%viewer_vect,errpetsc)

      END SUBROUTINE free


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! \brief Print the solver information
      !
      !
      SUBROUTINE printInfo(this)
      IMPLICIT NONE
      class(PetscSolver) :: this

      MatInfo info;
      PetscErrorCode errpetsc     ! Variable to check the error message during runtime

      ! call MatGetInfo(this%mtx, MAT_LOCAL, info, errpetsc);
      ! call PetscPrintf(MPI_COMM_WORLD, "Petsc solver:  nRow = %12d \n", this%nRow, errpetsc)
      ! call PetscPrintf(MPI_COMM_WORLD, "               nnz  = %12d \n\n", info.nz_allocated);

      call MatAssemblyBegin(this%mtx, MAT_FINAL_ASSEMBLY, errpetsc)
      CHKERRQ(errpetsc)
      call MatAssemblyEnd(this%mtx, MAT_FINAL_ASSEMBLY, errpetsc)
      CHKERRQ(errpetsc)

      call VecAssemblyBegin(this%rhsVec, errpetsc); !CHKERRQ(errpetsc);
      call VecAssemblyEnd(this%rhsVec, errpetsc); !CHKERRQ(errpetsc);

      ! call VecAssemblyBegin(reacVec, errpetsc); !CHKERRQ(errpetsc);
      ! call VecAssemblyEnd(reacVec, errpetsc); !CHKERRQ(errpetsc);

      !call MatView(this%mtx, PETSC_VIEWER_STDOUT_WORLD, errpetsc)
      !CHKERRQ(errpetsc)
      !call VecView(this%rhsVec, PETSC_VIEWER_STDOUT_WORLD, errpetsc)
      !CHKERRQ(errpetsc)

    !   ! call MatGetInfo(solverpetsc%mtx, MAT_LOCAL, &info);

    !   ! pcout << " mallocs      = " <<  info.mallocs << endl;
    !   ! pcout << " nz_allocated = " <<  info.nz_allocated << endl;
    !   ! pcout << " nz_used      = " <<  info.nz_used << endl;
    !   ! pcout << " nz_unneeded  = " <<  info.nz_unneeded << endl;

      END SUBROUTINE printInfo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! \brief Assembles the element matrix into global matrix
      !
      !
      SUBROUTINE assembleMatrix(this, RINDICES, CINDICES, KLOCAL)
      IMPLICIT NONE
      class(PetscSolver) :: this

      PetscErrorCode errpetsc     ! Variable to check the error message during runtime
      INTEGER :: ii, jj, nsize
      INTEGER, DIMENSION(:) :: RINDICES, CINDICES
      DOUBLE PRECISION, DIMENSION(:,:) :: KLOCAL

      nsize = size(RINDICES)

      DO ii=1, nsize
        DO jj=1, nsize
          call MatSetValue(this%mtx, RINDICES(ii), CINDICES(jj),
     1                     KLOCAL(ii,jj), ADD_VALUES);
        END DO
      END DO

      END SUBROUTINE assembleMatrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! \brief Assembles the element vector into global vector
      !
      !
      SUBROUTINE assembleVector(this, RINDICES, FLOCAL)
      IMPLICIT NONE
      class(PetscSolver) :: this

      PetscErrorCode errpetsc     ! Variable to check the error message during runtime
      INTEGER :: ii, nsize
      INTEGER, DIMENSION(:) :: RINDICES
      DOUBLE PRECISION, DIMENSION(:) :: FLOCAL

      nsize = size(RINDICES)

      DO ii=1, nsize
        call VecSetValue(this%rhsVec, RINDICES(ii), FLOCAL(ii),
     1                   ADD_VALUES, errpetsc)
      END DO

      END SUBROUTINE assembleVector


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! \brief Assembles the element matrix into global matrix
      !
      !
      SUBROUTINE assembleMatrixAndVector(this,
     1  RINDICES, CINDICES, KLOCAL, FLOCAL)

      IMPLICIT NONE
      class(PetscSolver) :: this

      PetscErrorCode errpetsc     ! Variable to check the error message during runtime
      INTEGER :: ii, jj, nsize
      INTEGER, DIMENSION(:) :: RINDICES, CINDICES
      DOUBLE PRECISION, DIMENSION(:,:) :: KLOCAL
      DOUBLE PRECISION, DIMENSION(:) :: FLOCAL

      nsize = size(RINDICES)

      DO ii=1, nsize
        call VecSetValue(this%rhsVec, RINDICES(ii), FLOCAL(ii),
     1                   ADD_VALUES, errpetsc)
        DO jj=1, nsize
          call MatSetValue(this%mtx, RINDICES(ii), CINDICES(jj),
     1                     KLOCAL(ii,jj), ADD_VALUES, errpetsc)
        END DO
      END DO

      END SUBROUTINE assembleMatrixAndVector


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! \brief Factorise the matrix
      !
      !
      SUBROUTINE factorise(this)
      IMPLICIT NONE
      class(PetscSolver) :: this

      PetscErrorCode errpetsc     ! Variable to check the error message during runtime

      IF(this%currentStatus /= ASSEMBLY_OK) THEN
        call PetscPrintf(PETSC_COMM_WORLD,
     1  "Assemble matrix first before solving it! \n", errpetsc)
        STOP " Aborting... in 'solverpetsc->factorise' "
      END IF

      this%currentStatus = FACTORISE_OK

      END SUBROUTINE factorise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! \brief Solve the matrix system
      !
      !
      SUBROUTINE solve(this)
      IMPLICIT NONE
      class(PetscSolver) :: this

      CHARACTER (LEN=100) :: charTemp

      PetscErrorCode errpetsc     ! Variable to check the error message during runtime
      KSPConvergedReason reason;
      PetscInt its;

      IF(this%currentStatus /= FACTORISE_OK) THEN
        call PetscPrintf(PETSC_COMM_WORLD,
     1  "Factorise matrix first before solving it! \n", errpetsc)
        STOP " Aborting... in 'solverpetsc->solve' "
      END IF

      call MatAssemblyBegin(this%mtx, MAT_FINAL_ASSEMBLY, errpetsc)
      CHKERRQ(errpetsc)
      call MatAssemblyEnd(this%mtx, MAT_FINAL_ASSEMBLY, errpetsc)
      CHKERRQ(errpetsc)

      ! call MatView(this%mtx, PETSC_VIEWER_STDOUT_WORLD, errpetsc)

      ! Assemble the solnVec
      call VecAssemblyBegin(this%solnVec, errpetsc)
      CHKERRQ(errpetsc)
      call VecAssemblyEnd(this%solnVec, errpetsc)
      CHKERRQ(errpetsc)
      call VecZeroEntries(this%solnVec, errpetsc)
      CHKERRQ(errpetsc)

      ! call VecView(this%solnVec, PETSC_VIEWER_STDOUT_WORLD, errpetsc)

      ! Assemble the rhsVec 
      call VecAssemblyBegin(this%rhsVec, errpetsc)
      CHKERRQ(errpetsc)
      call VecAssemblyEnd(this%rhsVec, errpetsc)
      CHKERRQ(errpetsc)

      ! call VecView(this%rhsVec, PETSC_VIEWER_STDOUT_WORLD, errpetsc)

      call PetscPrintf(PETSC_COMM_WORLD, 
     1 " Solving the matrix system \n", errpetsc)

      ! Solve the matrix system
      call KSPSolve(this%ksp, this%rhsVec, this%solnVec, errpetsc)
      CHKERRQ( errpetsc )

      call KSPGetConvergedReason(this%ksp, reason, errpetsc)

      IF(reason < 0) THEN
        call PetscPrintf(PETSC_COMM_WORLD, "Divergence.\n", errpetsc)
      ELSE
        call KSPGetIterationNumber(this%ksp, its, errpetsc)

        WRITE(charTemp,*) "Convergence in", its, " iterations.", "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      END IF

      END SUBROUTINE solve


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! \brief Solve the matrix system
      !
      !
      SUBROUTINE factoriseAndSolve(this)
      IMPLICIT NONE
      class(PetscSolver) :: this

      PetscErrorCode errpetsc     ! Variable to check the error message during runtime

      this%currentStatus = ASSEMBLY_OK

      call this%factorise()
      call this%solve()

      END SUBROUTINE factoriseAndSolve

      END MODULE Module_SolverPetsc
      
