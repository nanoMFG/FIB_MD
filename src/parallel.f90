MODULE parallel

  USE mpi

  IMPLICIT none

  integer  :: Np
  integer  :: myid

CONTAINS

  SUBROUTINE initparallel

    call initmpi

  END SUBROUTINE initparallel

  SUBROUTINE initmpi
    
    integer  ::  ierr

    call MPI_INIT( ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Np, ierr )
    !print *, 'Process ', myid, ' of ', Np, ' is alive'
    if (myid .eq. 0) print*, '----------------Parallel Environment Initialized, '
    if (myid .eq. 0) print*, Np, ' processors are alive and running ....-------------------'

  END SUBROUTINE initmpi

  SUBROUTINE closempi

    integer  :: ierr

    call MPI_FINALIZE(ierr)

  END SUBROUTINE closempi


  SUBROUTINE closeparallel
    
    call closempi

  END SUBROUTINE closeparallel

END MODULE parallel
