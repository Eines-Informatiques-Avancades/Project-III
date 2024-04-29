!> This module contains the subroutine `find_force_LJ` which calculates the forces and potential energy
!> between particles using the Lennard-Jones potential.
!> The module provides a public interface for the `find_force_LJ` subroutine.
!> The `find_force_LJ` subroutine takes in the positions of particles, box size, cutoff distance,
!> and returns the forces acting on the particles and the total potential energy.
!> The forces are calculated using periodic boundary conditions (PBC).
!> The module also contains a private section for internal use.
Module forces
   use pbc_module
   include 'mpif.h'
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
   Subroutine find_force_LJ(r, N, L, cutoff, F, pot, Ppot, nprocs, rank, counts_recv, displs_recv, imin, imax)
      Implicit none
      !include 'mpif.h'
      real(8), dimension(N, 3), intent(in) :: r
      real(8), intent(in) :: L, cutoff
      real(8) :: d, f_ij, pot_rank
      real(8), dimension(3) :: rij
      integer :: i, j, ierror
      integer, intent(in) :: nprocs, rank
      integer, intent(in) :: N
      real(8), dimension(N, 3), intent(out) :: F
      real(8), dimension(N, 3) :: F_new
      integer :: displs_recv(0:nprocs-1), counts_recv(0:nprocs)
      real(8), dimension(N/nprocs, 3) :: F_cut
      real(8), intent(out) :: pot, Ppot
      integer :: imin, imax, k, ii

      !     print*, "Starting forces rutine, rank", rank

!      call MPI_BCAST(r,N*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)

      pot = 0.d0
      Ppot = 0.d0
      F = 0.d0

      do i = imin, imax
         do j = 1, N
         if (i .ne. j) then
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
               !F(j, :) = F(j, :) - f_ij*rij(:)

               if (isnan(F(i, 1))) then
                  print *, "ERROR: Force is not a number"
                  print *, i, j
                  stop
               end if

               pot = pot + 4.d0*(1.d0/d**12 - 1.d0/d**6) - 4.d0*(1/cutoff**12 - 1.d0/cutoff**6)
               Ppot = Ppot + f_ij*d
            end if
         end if
         end do
      end do

      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      if (ierror.ne.0) then
         print*, "Error in MPI_BARRIER"
         stop
      end if

      pot_rank = pot
      pot = 0
      !call MPI_GATHER(pot, 1, MPI_DOUBLE_PRECISION, pot_list, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      !call MPI_GATHER(Ppot, 1, MPI_DOUBLE_PRECISION, Ppot_list, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_REDUCE(pot_rank, pot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      if (ierror.ne.0) then
         print*, "Error in MPI_REDUCE"
         stop
      end if

!      print*, "Ending forces rutine, rank", rank

      pot = pot/2

   End Subroutine find_force_LJ

End Module forces
