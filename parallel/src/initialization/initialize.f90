!> Module containing initialization rutines for molecular dynamics simulations.
!> Principal contributor: Anna Monclús 
module initialization

   implicit none

contains

   !########################################################################################################


!> Initializes the positions of N particles in a cubic box with density rho.
!> The box is centered at the origin.
!> The positions are stored in the array r.
!> The length of the box is calculated as L = (N/rho)**(1/3).
!> The box is divided in M = N**(1/3) cells.
!> The length of each cell is a = L/M.
!! @param N Number of particles.
!! @param rho Density of the system.
!! @param r Array containing the positions of the particles.
!! @param L Length of the box.
!! @param a Length of each cell.
!! @param x, y, z Coordinates of the particle.
!! @param ini Initial position of the box.
!! @param M Number of cells in each direction.
!! @param i, j, k Loop variables.
!! @param particle Index of the particle.

   Subroutine initialize_positions(N, rho, r)
      Implicit none
      integer, intent(in) :: N
      real(8), intent(in) :: rho
      real(8), dimension(N, 3), intent(out) :: r
      real(8) :: L, a, x, y, z, ini
      integer :: M, i, j, k, particle

      ! Length of the box
      L = (N/rho)**(1./3.)

      M = N**(1./3.)

      a = L/(M)

      ! Set the position of every particle
      particle = 1
      ini = -L/2.d0 ! to center the box at the origin (0, 0, 0)
      ! ini = 0.d0
      do i = 0, M - 1
         do j = 0, M - 1
            do k = 0, M - 1
               x = ini + i*a
               y = ini + j*a
               z = ini + k*a
               r(particle, :) = (/x, y, z/)
               !print *, particle, r(particle, :)
               particle = particle + 1
            end do
         end do
      end do
   End Subroutine

!> Initializes the velocities of N particles with a given velocity v_ini.
!> The velocities are stored in the array v.
!> The velocities are chosen randomly from a bimodal distribution.
!! @param N Number of particles.
!! @param v_ini Initial velocity of the particles.
!! @param v Array containing the velocities of the particles.
!! @param v_i Velocity of the particle.
!! @param rand Random number.
!! @param particle Index of the particle.
!! @param i Loop variable.

   Subroutine initialize_velocities(N, v_ini, v)
      Implicit none
      integer, intent(in) :: N
      real(8), dimension(N, 3), intent(out) :: v
      real(8) :: v_ini, v_i, rand
      integer :: particle, i

      ! Set the velocity of every particle using a bimodal distribution
      do particle = 1, N
         do i = 1, 3
            call random_number(rand)
            if (rand < 0.5d0) then
               v(particle, i) = +v_ini
            else
               v(particle, i) = -v_ini
            end if

         end do
      end do
   End Subroutine

   !########################################################################################################

!> Distributes the particles among the processors.
!! @param N Number of particles.
!! @param rank Rank of the processor.
!! @param nprocs Number of processors.
!! @param imin, imax Indexes of the particles to be handled by the processor.
!! @param chunklength Number of particles to be handled by each processor.
!! @param uneven_parts Number of processors that will handle an extra particle.
!! @param i Loop variable.

   subroutine distribute_particles(N, rank, nprocs, imin, imax)
      implicit none
      integer, intent(in) :: N
      integer, intent(in) :: rank, nprocs
      integer :: chunklength, uneven_parts, imin, imax

      !
      chunklength = N/nprocs
      uneven_parts = mod(N, nprocs) ! 5

      print *, "uneven_parts = ", uneven_parts

      imin = rank*chunklength + 1
      imax = imin + chunklength - 1

      if (rank >= (nprocs - uneven_parts)) then
         print *, "Im rank", rank, "and . They already reparted", rank - (nprocs - uneven_parts)

         imin = imin + rank - (nprocs - uneven_parts)
         imax = imin + chunklength
      end if

   end subroutine

end module initialization
