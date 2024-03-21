!> Module containing various integrators for molecular dynamics simulations.
module integrators

   use Forces_Module
   use pbc_module

   implicit none

contains

!> Perform a time step using the velocity Verlet integration method.
!! Calculates new positions and velocities for particles based on forces and previous positions/velocities.
!! @param r Input/Output array containing particle positions.
!! @param vel Input/Output array containing particle velocities.
!! @param pot Output potential energy.
!! @param N Number of particles.
!! @param L Box size.
!! @param cutoff Cutoff distance for LJ potential.
!! @param dt Time step size.
   subroutine time_step_vVerlet(r, vel, pot, N, L, cutoff, dt)
      implicit none
      integer, intent(in) :: N                      !< Number of particles
      real(8), dimension(N, 3), intent(inout) :: r  !< Particle positions
      real(8), dimension(N, 3), intent(inout) :: vel  !< Particle velocities
      real(8), intent(out) :: pot                  !< Potential energy
      real(8), intent(in) :: dt, L, cutoff          !< Time step size, box size, cutoff distance
      real(8), dimension(N, 3) :: F                 !< Forces
      integer :: i                                  !< Loop variable

      ! Calculate forces and potential energy using LJ potential
      call find_force_LJ(r, N, L, cutoff, F, pot)

      ! Update positions and velocities using velocity Verlet integration
      do i = 1, N
         r(i, :) = r(i, :) + vel(i, :)*dt + 0.5*F(i, :)*dt*dt

         ! Apply periodic boundary conditions
         do while (any(r(i, :) > L/2.) .or. any(r(i, :) < -L/2.))
            call pbc(r(i, :), L, size(r(i, :)))   !< Apply periodic boundary conditions using the pbc subroutine
         end do

         vel(i, :) = vel(i, :) + F(i, :)*0.5*dt
      end do

      ! Recalculate forces after updating positions
      call find_force_LJ(r, N, L, cutoff, F, pot)

      ! Update velocities using the updated forces
      do i = 1, N
         vel(i, :) = vel(i, :) + F(i, :)*0.5*dt
      end do

   end subroutine time_step_vVerlet

!> Perform a time step using the Euler method with periodic boundary conditions.
!! Calculates new positions and velocities for particles based on forces and previous positions/velocities.
!! @param r_in Input array containing initial particle positions.
!! @param r_out Output array containing updated particle positions.
!! @param vel Input/Output array containing particle velocities.
!! @param N Number of particles.
!! @param L Box size.
!! @param cutoff Cutoff distance for LJ potential.
!! @param dt Time step size.
!! @param pot Output potential energy.
   subroutine time_step_Euler_pbc(r_in, r_out, vel, N, L, cutoff, dt, pot)
      implicit none
      integer, intent(in) :: N                         !< Number of particles
      real(8), dimension(N, 3), intent(in) :: r_in     !< Initial particle positions
      real(8), dimension(N, 3), intent(out) :: r_out   !< Updated particle positions
      real(8), dimension(N, 3) :: vel                   !< Particle velocities
      real(8), dimension(N, 3) :: F                     !< Forces
      real(8), intent(in) :: L, dt, cutoff              !< Box size, time step size, cutoff distance
      integer :: i, counter                             !< Loop variables
      real(8) :: pot                                    !< Potential energy

      ! Calculate forces and potential energy using LJ potential
      call find_force_LJ(r_in, N, L, cutoff, F, pot)

      ! Update positions and velocities using Euler integration and apply periodic boundary conditions
      do i = 1, N
         r_out(i, :) = r_in(i, :) + vel(i, :)*dt + 0.5*F(i, :)*dt*dt
         vel(i, :) = vel(i, :) + F(i, :)*dt

         ! Apply periodic boundary conditions
         do while (any(r_out(i, :) .gt. L/2.) .or. any(r_out(i, :) .lt. (-L/2.)))
            call pbc(r_out, L, size(r_out))
         end do
      end do

   end subroutine time_step_Euler_pbc

!> Generate random numbers following a Box-Muller transformation.
!! Generates normally distributed random numbers using the Box-Muller transformation.
!! @param ndat Number of data points (must be even).
!! @param xnums Output array containing the generated random numbers.
!! @param sigma Standard deviation of the normal distribution.
   subroutine BM(ndat, xnums, sigma)
      implicit none
      integer, intent(in) :: ndat                     !< Number of data points (must be even)
      integer :: i                                     !< Loop variable
      real(8), dimension(ndat), intent(out) :: xnums  !< Output array containing generated random numbers
      real(8) :: r, phi, x1, x2                        !< Variables for the Box-Muller transformation
      real(8), intent(in) :: sigma                     !< Standard deviation of the normal distribution
      real(8), parameter :: pi = 4.d0*atan(1.d0)     !< Mathematical constant pi

      ! Generate random numbers using Box-Muller transformation
      do i = 1, ndat, 2
         r = sqrt(-2.d0*log(1.d0 - rand()))
         phi = 2.d0*pi*rand()
         x1 = r*cos(phi)
         x2 = r*sin(phi)

         ! Store the generated random numbers in the output array
         if (i .ne. ndat) then
            xnums(i) = x1*sigma
            xnums(i + 1) = x2*sigma
         end if
      end do

   end subroutine BM

!> Perform a time step using the Verlet integration method.
!! Calculates new positions and velocities for particles based on forces and previous positions/velocities.
!! @param r Input/Output array containing current particle positions.
!! @param rold Input/Output array containing previous particle positions.
!! @param vel Input/Output array containing particle velocities.
!! @param N Number of particles.
!! @param L Box size.
!! @param cutoff Cutoff distance for LJ potential.
!! @param dt Time step size.
!! @param pot Output potential energy.
   subroutine time_step_Verlet(r, rold, vel, N, L, cutoff, dt, pot)
      implicit none
      integer, intent(in) :: N                           !< Number of particles
      real(8), dimension(N, 3), intent(inout) :: r       !< Current particle positions
      real(8), dimension(N, 3), intent(inout) :: rold    !< Previous particle positions
      real(8), dimension(N, 3) :: F                       !< Forces
      real(8), dimension(N, 3), intent(inout) :: vel     !< Particle velocities
      real(8), intent(in) :: dt, cutoff, L               !< Time step size, cutoff distance, box size
      real(8) :: pot                                      !< Potential energy
      integer :: i, j                                     !< Loop variables
      real(8), dimension(N, 3) :: raux, roldaux           !< Temporary variable for position update

      ! Calculate forces and potential energy using LJ potential
      call find_force_LJ(r, N, L, cutoff, F, pot)

      ! Save current positions to temporary array
      roldaux = r

      ! Update positions and velocities using Verlet integration
      do i = 1, N
         do j = 1, 3
            r(i, j) = 2*roldaux(i, j) - rold(i, j) + F(i, j)*dt*dt

            ! Apply periodic boundary conditions
            raux(i, j) = r(i, j) - rold(i, j)
            call pbc(raux(i, j), L, size(raux))

            ! Update velocities
            vel(i, j) = raux(i, j)/(2*dt)
         end do
      end do

      ! Update previous positions for the next time step
      rold = roldaux

   end subroutine time_step_Verlet

!########################################################################################################

   Subroutine kinetic_energy(vel, K_energy, N)
      Implicit none
      integer, intent(in) :: N
      real(8), dimension(N, 3) :: vel
      integer :: i, k
      real(8) :: K_energy

      K_energy = 0
      ! for each particle
      do i = 1, N
         ! loop over coordinates
         do k = 1, 3
            K_energy = K_energy + 0.5*vel(i, k)*vel(i, k)
         end do
      end do
   End Subroutine

! #########################################################################################################

   Function inst_temp(N, K_energy)
      Implicit none
      integer :: N, N_f
!        real(8), parameter :: k_b = 1.380649e-23
      real(8) :: K_energy, inst_temp

      N_f = 3*N - 3
      ! inst_temp = 2.d0/(N_f * k_b)*K_energy
      inst_temp = 2.d0/(N_f)*K_energy
      Return
   End Function

end module

!#############################################################

Subroutine momentum(vel, p, N)
   Implicit none
   real(8), dimension(N, 3) :: vel
   real(8), dimension(3) :: total_p
   integer :: N, i
   real(8), intent(out) :: p

   total_p(:) = 0

   ! Accumulate p
   do i = 1, N
      total_p(:) = total_p(:) + vel(i, :)
   end do

   ! Produce the module
   p = (total_p(1)**2 + total_p(2)**2 + total_p(3)**2)**(1./2.)

End Subroutine
