!
! Program for paritioning the mesh
!

      PROGRAM MeshPartition

      use partition

      IMPLICIT NONE

      DOUBLE PRECISION :: tstart, tend
      !integer :: this_mpi_proc, n_mpi_proc

      call MPI_Init(errpar)

      call MPI_Comm_size(MPI_COMM_WORLD, n_mpi_procs, errpar);
      call MPI_Comm_rank(MPI_COMM_WORLD, this_mpi_proc, errpar);

      if(n_mpi_procs .EQ. 1) then
        WRITE(*,*) "This program should be called with at least two processors "
        STOP "Aborting now ..."
      end if

      WRITE(*,*) " this_mpi_proc = ", this_mpi_proc, "\n"

      !call CPU_TIME(tstart)MPI_WTIME

      call main()

      !call CPU_TIME(tend)
      !WRITE(*,*) "That took ", (tend-tstart), "seconds"

      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) "Program is successful"
      WRITE(*,*) " "
      WRITE(*,*) " "

      call MPI_Finalize(errpar)

      END PROGRAM MeshPartition
