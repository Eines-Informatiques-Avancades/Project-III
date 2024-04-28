module initialization

   implicit none

contains

   !########################################################################################################

   Subroutine initialize_positions(N, rho, r)
      ! """"
      ! Calculates the positions r of N particles in a SC structure
      ! INPUTS: N, rho
      ! OUTPUT: r(N, 3)
      ! """"
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

   Subroutine initialize_velocities(N, v_ini, v)
      ! """"
      ! Initializes the velocities v of N particles, all with v_ini velocity
      ! INPUTS: N. v_ini
      ! OUTPUT: v(N,3)
      ! """"
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

end module initialization
