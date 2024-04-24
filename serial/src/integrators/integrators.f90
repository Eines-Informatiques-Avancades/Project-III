!> Module containing various integrators for molecular dynamics simulations.
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
   subroutine time_step_vVerlet(r, vel, pot, N, L, cutoff, dt, Ppot, nprocs, rank)
      implicit none
      integer, intent(in) :: N                      !< Number of particles
      real(8), dimension(N, 3), intent(inout) :: r  !< Particle positions
      real(8), dimension(N, 3), intent(inout) :: vel  !< Particle velocities
      real(8), intent(out) :: pot, Ppot                 !< Potential energy
      real(8), intent(in) :: dt, L, cutoff          !< Time step size, box size, cutoff distance
      real(8), dimension(N, 3) :: F                 !< Forces
      integer :: i, nprocs, rank                                  !< Loop variable
      ! Calculate forces and potential energy using LJ potential
      call find_force_LJ(r, N, L, cutoff, F, pot, Ppot, nprocs, rank)

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
      call find_force_LJ(r, N, L, cutoff, F, pot, Ppot, nprocs, rank)

      ! Update velocities using the updated forces
      do i = 1, N
         vel(i, :) = vel(i, :) + F(i, :)*0.5*dt
      end do

   end subroutine time_step_vVerlet



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

end module integrators

!> Calculate the total momentum of particles.
   !! @param vel Array containing particle velocities.
   !! @param p Output variable for the total momentum.
   !! @param N Number of particles.
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

! Subroutine therm_Andersen(vel, nu, sigma_gaussian, N)
!    Implicit none
!    integer :: i, N
!    real(8) :: rand, nu, sigma_gaussian
!    real(8), dimension(N, 3) :: vel
!    real(8), dimension(2) :: xnums

!    do i = 1, N
!       call random_number(rand)
!       if (rand .lt. nu) then
!          call BM(2, xnums, sigma_gaussian)
!          !print*, "xnums: ", xnums
!          vel(i, 1) = xnums(1)
!          vel(i, 2) = xnums(2)
!          call BM(2, xnums, sigma_gaussian)
!          vel(i, 3) = xnums(1)
!       end if
!    end do
! !        print*, vel
! End Subroutine

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