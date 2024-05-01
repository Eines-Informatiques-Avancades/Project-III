!> This module contains the subroutine `find_force_LJ` which calculates the forces and potential energy
!> between particles using the Lennard-Jones potential.
!> Principal contributor: Abert Plazas 
Module forces

   use pbc_module

   Implicit none

   Private
   Public :: find_force_LJ

Contains

   !> Calculates the forces and potential energy between particles using the Lennard-Jones potential (serial implementation).
    !!
    !! This subroutine calculates the forces and potential energy between particles in a system using the Lennard-Jones potential. The Lennard-Jones potential is a pairwise potential that models the interaction between neutral atoms or molecules. It is commonly used to simulate the behavior of noble gases and other non-reactive systems.
    !!
    !! Notes:
    !!   - The positions of the particles are given as a 3D array, where each row represents the position of a particle in Cartesian coordinates.
    !!   - The forces acting on each particle are calculated and stored in the `F` array.
    !!   - The potential energy of the system is calculated and stored in the `pot` variable.
    !!   - The subroutine uses periodic boundary conditions to handle particles that are close to the edges of the simulation box.
    !!   - The subroutine assumes that the `pbc` subroutine is defined elsewhere in the code, which handles the periodic boundary conditions.
    !!   - The subroutine assumes that the `isnan` function is available to check for NaN values.
    !!
   Subroutine find_force_LJ(r, N, L, cutoff, F, pot, Ppot)
      Implicit none
      real(8), dimension(N, 3), intent(in) :: r !> Array with the positions of the particles in the system.
      real(8), intent(in) :: L, cutoff !> Length of the simulation box and the cutoff distance for the Lennard-Jones potential.
      real(8) :: d, f_ij !> Variables to store the distance between particles and the force between particles.
      real(8), dimension(3) :: rij !> Vector between particles.
      integer :: i, j !> Loop indices.
      integer, intent(in) :: N !> Number of particles in the system.
      real(8), dimension(N, 3), intent(out) :: F !> Array with the forces acting on each particle.
      real(8), intent(out) :: pot, Ppot !> Variables to store the potential energy of the system and the pressure due to the potential.

      pot = 0.d0
      Ppot = 0.d0
      F = 0.d0

      do i = 1, N
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
                  print *, i, j
                  stop
               end if

               pot = pot + 4.d0*(1.d0/d**12 - 1.d0/d**6) - 4.d0*(1/cutoff**12 - 1.d0/cutoff**6)
               Ppot = Ppot + f_ij*d
            end if

         end do
      end do

   End Subroutine find_force_LJ

End Module forces
