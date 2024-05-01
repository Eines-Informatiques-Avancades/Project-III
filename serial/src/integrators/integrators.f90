!> Module containing various integrators for molecular dynamics simulations.
!> Principal contirbutor: Aina Gaya
module integrators

   use forces
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
   subroutine time_step_vVerlet(r, vel, pot, N, L, cutoff, dt, Ppot)
      implicit none
      integer, intent(in) :: N                      !< Number of particles
      real(8), dimension(N, 3), intent(inout) :: r  !< Particle positions
      real(8), dimension(N, 3), intent(inout) :: vel  !< Particle velocities
      real(8), intent(out) :: pot, Ppot                 !< Potential energy
      real(8), intent(in) :: dt, L, cutoff          !< Time step size, box size, cutoff distance
      real(8), dimension(N, 3) :: F                 !< Forces
      integer :: i                                  !< Loop variable
      ! Calculate forces and potential energy using LJ potential
      call find_force_LJ(r, N, L, cutoff, F, pot, Ppot)

      ! Update positions and velocities using velocity Verlet integration
      do i = 1, N
         r(i, :) = r(i, :) + vel(i, :)*dt + 0.5*F(i, :)*dt*dt

         ! Apply periodic boundary conditions
         do while (any(r(i, :) > L/2.) .or. any(r(i, :) < -L/2.))
            call pbc_mic(r(i, :), L, size(r(i, :))) !< Apply periodic boundary conditions using the pbc subroutine
         end do

         vel(i, :) = vel(i, :) + F(i, :)*0.5*dt
      end do

      ! Recalculate forces after updating positions
      call find_force_LJ(r, N, L, cutoff, F, pot, Ppot)

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
      real(8) :: pot, Ppot                                    !< Potential energy

      ! Calculate forces and potential energy using LJ potential
      call find_force_LJ(r_in, N, L, cutoff, F, pot, Ppot)

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
      real(8) :: pot, Ppot                                      !< Potential energy
      integer :: i, j                                     !< Loop variables
      real(8), dimension(N, 3) :: raux, roldaux           !< Temporary variable for position update

      ! Calculate forces and potential energy using LJ potential
      call find_force_LJ(r, N, L, cutoff, F, pot, Ppot)

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

!! Calculate the kinetic energy of particles.
!! @param vel Array containing particle velocities.
!! @param K_energy Output variable for the kinetic energy.
!! @param N Number of particles.
   subroutine kinetic_energy(vel, K_energy, N)
      implicit none
      integer, intent(in) :: N                           !< Number of particles
      real(8), dimension(N, 3), intent(in) :: vel        !< Array containing particle velocities
      real(8), intent(out) :: K_energy                   !< Output variable for the kinetic energy
      integer :: i, k                                    !< Loop variables

      K_energy = 0

      ! Calculate kinetic energy for each particle
      do i = 1, N
         ! Loop over coordinates
         do k = 1, 3
            K_energy = K_energy + 0.5*vel(i, k)*vel(i, k)
         end do
      end do

   end subroutine kinetic_energy

   !> Calculate the instantaneous temperature of the system.
      !! @param N Number of particles.
      !! @param K_energy Kinetic energy of the particles.
      !! @return Instantaneous temperature of the system.
   function inst_temp(N, K_energy) result(temp)
      implicit none
      integer, intent(in) :: N                           !< Number of particles
      real(8), intent(in) :: K_energy                    !< Kinetic energy of the particles
      real(8) :: temp                                    !< Instantaneous temperature of the system

      temp = 2.0d0*K_energy/(3*N - 3)

   end function inst_temp

! > Calculate the total momentum of particles.
!    ! @param vel Array containing particle velocities.
!    ! @param p Output variable for the total momentum.
!    ! @param N Number of particles.
   subroutine momentum(vel, p, N)
      implicit none
      integer, intent(in) :: N                           !< Number of particles
      real(8), dimension(N, 3), intent(in) :: vel        !< Array containing particle velocities
      real(8), dimension(3) :: total_p                    !< Total momentum
      integer :: i                                        !< Loop variable
      real(8), intent(out) :: p                           !< Output variable for the total momentum

      total_p(:) = 0

      ! Accumulate momentum
      do i = 1, N
         total_p(:) = total_p(:) + vel(i, :)
      end do

      ! Calculate the magnitude of the total momentum
      p = sqrt(total_p(1)**2 + total_p(2)**2 + total_p(3)**2)

   end subroutine momentum

!#################################################################

   Subroutine therm_Andersen(vel, nu, sigma_gaussian, N)
      Implicit none
      integer :: i, N
      real(8) :: rand, nu, sigma_gaussian
      real(8), dimension(N, 3) :: vel
      real(8), dimension(2) :: xnums

      do i = 1, N
         call random_number(rand)
         if (rand .lt. nu) then
            call BM(2, xnums, sigma_gaussian)
            !print*, "xnums: ", xnums
            vel(i, 1) = xnums(1)
            vel(i, 2) = xnums(2)
            call BM(2, xnums, sigma_gaussian)
            vel(i, 3) = xnums(1)
         end if
      end do
      !        print*, vel
   End Subroutine
end module integrators
!#################################################################

! Subroutine BM(ndat, xnums, sigma)
!    Implicit none
!    Integer ::  ndat, i
!    real(8), dimension(ndat) :: xnums
!    real(8) :: r, phi, x1, x2, sigma
!    real(8), parameter :: pi = 4.d0*atan(1.d0)
! !     ATENCIÓ! Es generen 2ndat numeros
!    Do i = 1, ndat, 2
!       r = sqrt(-2.d0*log(1.d0 - rand()))
!       phi = 2.d0*pi*rand()
!       x1 = r*cos(phi)
!       x2 = r*sin(phi)
!       if (i .ne. ndat) then ! Ens assegurem que no haguem acabat la llista
!          xnums(i) = x1*sigma
!          xnums(i + 1) = x2*sigma
!       end if
!    end do
!    return
! end Subroutine

