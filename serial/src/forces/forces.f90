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

    Implicit none
    
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
    Subroutine find_force_LJ(r, N, L, cutoff, F, pot)
        Implicit none
        real(8), dimension(N, 3), intent(in) :: r
        real(8), intent(in) :: L, cutoff
        real(8) :: d, f_ij
        real(8), dimension(3) :: d_r
        integer :: i, j
        integer, intent(in) :: N
        real(8), dimension(N, 3), intent(out) :: F
        real(8), intent(out) :: pot

        pot = 0.d0

        F = 0.d0

        do i = 1, N
            do j = i + 1, N
                d_r(:) = r(i, :)

                do while (any(d_r(:) .gt. L/2.) .or. (any(d_r(:) .lt. (-L/2.))))
                    call pbc(d_r, L, size(d_r))
                end do

                d = (d_r(1)**2 + d_r(2)**2 + d_r(3)**2)**(1.d0/2.d0)
                if (d .le. cutoff) then
                    f_ij = 48.d0/d**14 - 24.d0/d**8
                    F(i, :) = F(i, :) + f_ij*d_r(:)
                    F(j, :) = F(j, :) - f_ij*d_r(:)

                    if (isnan(F(i, 1))) then
                        print *, i, j
                        stop
                    end if

                    pot = pot + 4.d0*(1.d0/d**12 - 1.d0/d**6) - 4.d0*(1/cutoff**12 - 1.d0/cutoff**6)

                end if

            end do
        end do

    End Subroutine find_force_LJ
    
End Module forces
