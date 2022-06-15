program e_fortran
#include <petsc/finclude/petscsys.h>
  use petscsys

  implicit none

  PetscErrorCode  :: ierr
  PetscMPIInt :: rank,size
  PetscInt :: i
  PetscReal :: localval,globalsum
  character(len=80) :: outputString

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr); CHKERRQ(ierr)

  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr); CHKERRA(ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRA(ierr)

  localval = 1.0
  do i = 1,rank
    localval = localval/i 
  end do

  call MPI_Reduce(localval, globalsum, 1, MPI_DOUBLE, MPI_SUM, &
                   0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)

  write(outputString,*) 'e is about:',globalsum
  call PetscPrintf(PETSC_COMM_WORLD,outputString,ierr); CHKERRA(ierr)
  
  call PetscFinalize(ierr)
end program