!> @module Forces_Module
!! This module contains the subroutine `find_force_LJ` which calculates the forces and potential energy
!! between particles using the Lennard-Jones potential.
!! The module provides a public interface for the `find_force_LJ` subroutine.
!! The `find_force_LJ` subroutine takes in the positions of particles, box size, cutoff distance,
!! and returns the forces acting on the particles and the total potential energy.
!! The forces are calculated using periodic boundary conditions (PBC).
!! The module also contains a private section for internal use.
!! @endmodule
Module forces
    use pbc_module
    !include 'mpif.h'
    Private
    Public :: find_force_LJ
 
 Contains
 
    !> Calculates the forces and potential energy between particles using the Lennard-Jones potential.
     !!
     !! This subroutine calculates the forces and potential energy between particles in a system using the Lennard-Jones potential. The Lennard-Jones potential is a pairwise potential that models the interaction between neutral atoms or molecules. It is commonly used to simulate the behavior of noble gases and other non-reactive systems.
     !!
     !! Parameters:
     !!   r: real(8), dimension(N, 3), intent(in)
     !!     - The positions of the particles in the system.
     !!   N: integer, intent(in)
     !!     - The number of particles in the system.
     !!   L: real(8), intent(in)
     !!     - The length of the simulation box.
     !!   cutoff: real(8), intent(in)
     !!     - The cutoff distance for the Lennard-Jones potential.
     !!   F: real(8), dimension(N, 3), intent(out)
     !!     - The forces acting on each particle.
     !!   pot: real(8), intent(out)
     !!     - The potential energy of the system.
     !!
     !! Notes:
     !!   - The positions of the particles are given as a 3D array, where each row represents the position of a particle in Cartesian coordinates.
     !!   - The forces acting on each particle are calculated and stored in the `F` array.
     !!   - The potential energy of the system is calculated and stored in the `pot` variable.
     !!   - The subroutine uses periodic boundary conditions to handle particles that are close to the edges of the simulation box.
     !!   - The subroutine assumes that the `pbc` subroutine is defined elsewhere in the code, which handles the periodic boundary conditions.
     !!   - The subroutine assumes that the `isnan` function is available to check for NaN values.
     !!
    Subroutine find_force_LJ(r, N, L, cutoff, F, pot, Ppot, nprocs, rank)
       Implicit none
       include 'mpif.h'
       real(8), dimension(N, 3), intent(in) :: r
       real(8), intent(in) :: L, cutoff
       real(8) :: d, f_ij
       real(8), dimension(3) :: rij
       integer :: i, j, ierror
       integer, intent(in) :: nprocs, rank
       integer, intent(in) :: N
       real(8), dimension(N, 3), intent(out) :: F
       integer, allocatable :: displs_recv(:),counts_recv(:)
       real(8), dimension(N/nprocs, 3) :: F_cut
       real(8), intent(out) :: pot, Ppot
       real(8), dimension(:), allocatable :: pot_list, Ppot_list
       integer :: imin, imax
       !call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
       !call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
       allocate(pot_list(nprocs))
       allocate(Ppot_list(nprocs))
       pot = 0.d0
       Ppot = 0.d0
       F = 0.d0
       pot_list = 0.d0
       Ppot_list = 0.d0

       imin = 1 + (rank * N/nprocs)
       imax = (N/nprocs) + (rank * N/nprocs)

       do i = imin, imax
          do j = i + 1, N
             rij(1) = r(i, 1) - r(j, 1)
             rij(2) = r(i, 2) - r(j, 2)
             rij(3) = r(i, 3) - r(j, 3)
 
             do while (any(rij(:) .gt. L/2.) .or. (any(rij(:) .lt. (-L/2.))))
                call pbc_mic(rij, L, size(rij))
             end do
 
             d = (rij(1)**2 + rij(2)**2 + rij(3)**2)**(1.d0/2.d0)
             if (d .le. cutoff) then
                f_ij = 48.d0/d**14 - 24.d0/d**8
                F(i, :) = F(i, :) + f_ij*rij(:)
                F(j, :) = F(j, :) - f_ij*rij(:)
 
                if (isnan(F(i, 1))) then
                   print*, "ERROR: Force is not a number"
                   print *, i, j
                   stop
                end if
 
                pot = pot + 4.d0*(1.d0/d**12 - 1.d0/d**6) - 4.d0*(1/cutoff**12 - 1.d0/cutoff**6)
                Ppot = Ppot + f_ij*d
             end if
 
          end do
       end do
       F_cut = F(imin:imax,:)
       
       allocate(counts_recv(nprocs))
       counts_recv(1) = 25
       counts_recv(2) = 25
       counts_recv(3) = 25
       counts_recv(4) = 25
       counts_recv(5) = 25
       allocate(displs_recv(nprocs))
       displs_recv(1) = 0
       displs_recv(2) = 25
       displs_recv(3) = 50
       displs_recv(4) = 75
       displs_recv(5) = 100
       call MPI_BARRIER(MPI_COMM_WORLD, ierror)
       !call MPI_ALLGATHERV(pot, 1, MPI_DOUBLE_PRECISION, pot_list, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
       !call MPI_ALLGATHERV(Ppot, 1, MPI_DOUBLE_PRECISION, Ppot_list, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)      
       call MPI_ALLGATHERV(F(imin:imax,1), int(N/nprocs), MPI_DOUBLE_PRECISION, F(:,1), counts_recv, displs_recv, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
       call MPI_ALLGATHERV(F(imin:imax,2), int(N/nprocs), MPI_DOUBLE_PRECISION, F(:,2), counts_recv, displs_recv, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
       call MPI_ALLGATHERV(F(imin:imax,3), int(N/nprocs), MPI_DOUBLE_PRECISION, F(:,3), counts_recv, displs_recv, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
       pot = sum(pot_list)
       Ppot = sum(Ppot_list)


       deallocate(counts_recv, displs_recv)
      
    End Subroutine find_force_LJ
 
 End Module forces
 