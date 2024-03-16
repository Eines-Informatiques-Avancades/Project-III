
module subroutines_md

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

      L = (N/rho)**(1./3.)

      M = N**(1./3.)

      a = L/(M)

      ! Set the position of every particle
      particle = 1
      ini = -L/2.d0
      ! ini = 0.d0
      do i = 0, M - 1
         do j = 0, M - 1
            do k = 0, M - 1
               x = ini + i*a
               y = ini + j*a
               z = ini + k*a
               r(particle, :) = (/x, y, z/)
               print *, particle, r(particle, :)
               particle = particle + 1
            end do
         end do
      end do
   End Subroutine

!########################################################################################################

   Subroutine time_step_vVerlet(r, vel, pot, N, L, cutoff, dt)
! """"
! Calculates a timestep using the velicity verlet algorithm
! INPUTS: r, vel, N, L, cutoff, dt
! OUTPUT: r, vel, pot
! """"
      Implicit none
      integer, intent(in) :: N
      real(8), dimension(N, 3), intent(inout) :: r, vel
      real(8), intent(in) :: dt, L, cutoff
      real(8), dimension(N, 3) :: F
      real(8), intent(out) ::  pot
      integer :: i

      Call find_force_LJ(r, N, L, cutoff, F, pot)

      do i = 1, N
         r(i, :) = r(i, :) + vel(i, :)*dt + 0.5*F(i, :)*dt*dt

         do while (any(r(i, :) > L/2.) .or. any(r(i, :) < -L/2.))
            ! Apply periodic boundary conditions using the pbc subroutine
            call pbc(r(i, :), L, size(r(i, :)))
         end do

         vel(i, :) = vel(i, :) + F(i, :)*0.5*dt
      end do

      Call find_force_LJ(r, N, L, cutoff, F, pot)

      do i = 1, N
         vel(i, :) = vel(i, :) + F(i, :)*0.5*dt
      end do

   End Subroutine

!########################################################################################################

   Subroutine pbc(vector, L, D)
! """"
! Applies periodic boundary conditions to a N-dimensional system
! INPUTS: vector, L, D (dimension of the system)
! OUTPUT: vector
! """"
      Implicit none
      integer :: i
      integer, intent(in) :: D
      real(8), dimension(D), intent(inout) :: vector
      real(8), intent(in) :: L

      do i = 1, D
         if (vector(i) .gt. L/2.) then
            vector(i) = vector(i) - L
            if (abs(vector(i)) .gt. 1000) then
               print *, "YOU HAVE INESTABILITIES!"
               stop
            end if
         else if (vector(i) .lt. (-L/2.)) then
            vector(i) = vector(i) + L
            if (abs(vector(i)) .gt. 1000) then
               print *, "YOU HAVE INESTABILITIES!"
               stop
            end if
         end if
      end do

      Return
   End Subroutine

!########################################################################################################

   Subroutine find_force_LJ(r, N, L, cutoff, F, pot)
! """"
! Calculates the forces applied to each particle of the system
! INPUTS: r, Nm L, cutoff, F
! OUTPUT: pot, F
! """"
      Implicit none
      real(8), dimension(N, 3), intent(in) :: r
      real(8), intent(in) :: L, cutoff
      real(8) :: d, f_ij
      real(8), dimension(3) :: d_r
      integer :: i, j, k
      integer, intent(in) :: N
      real(8), dimension(N, 3), intent(out) :: F
      real(8), intent(out) :: pot

      pot = 0.d0

      F = 0.d0

      do i = 1, N
         do j = i + 1, N
            d_r(:) = r(i, :) - r(j, :)

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

   End Subroutine

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

!###########################################################################################################

   Subroutine initialize_velocities(N, absV, vel)
      Implicit none
      integer, intent(in) :: N
      real(8) :: absV, rnd
      integer :: i, j
      real(8), dimension(N, 3) :: vel

      do i = 1, N
         do j = 1, 3
            call random_number(rnd)
            if (rnd .ge. 0.5) then
               vel(i, j) = +absV
            else if (rnd .lt. 0.5) then
               vel(i, j) = -absV
            end if
         end do
      end do

   End Subroutine

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

! ###############################################################

   Subroutine time_step_Euler_pbc(r_in, r_out, vel, N, L, cutoff, dt, pot)
      Implicit none
      integer, intent(in) :: N
      real(8), dimension(N, 3), intent(in) :: r_in
      real(8), dimension(N, 3), intent(out) :: r_out
      real(8), dimension(N, 3) :: vel, F
      real(8), intent(in) :: L, dt
      integer :: i, counter
      real(8) :: cutoff, pot

      Call find_force_LJ(r_in, N, L, cutoff, F, pot)

      do i = 1, N
         r_out(i, :) = r_in(i, :) + vel(i, :)*dt + 0.5*F(i, :)*dt*dt
         vel(i, :) = vel(i, :) + F(i, :)*dt

         do while (any(r_out(i, :) .gt. L/2.) .or. (any(r_out(i, :) .lt. (-L/2.))))
            call pbc(r_out, L, size(r_out))
         end do

      end do

   End Subroutine

! ##################################################################

   Function W(N, r, rho, T, L, cutoff)
      Implicit none
      integer :: N, i, j, step
      real(8) :: suma, W, T, L, rho, d, cutoff, f_ij
      real(8), dimension(N, 3) :: r
      real(8), dimension(3) :: r_ij

      suma = 0

      do i = 1, N
         do j = i + 1, N
            r_ij(:) = r(i, :) - r(j, :)

            do while (any(r_ij(:) .gt. L/2.) .or. (any(r_ij(:) .lt. (-L/2.))))
               call pbc(r_ij, L, size(r_ij))
            end do

            d = (r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2)**(1.d0/2.d0)
            if (d .le. cutoff) then
               f_ij = 48.d0/d**13 - 24.d0/d**7
               !        print*, "f", f_ij
               suma = suma + f_ij*d
            end if

         end do
      end do

      ! average = sum / (N*N)

      W = (1./3.)*suma

   end Function

! ###################################################################

   Subroutine mean_sq_distance(r, r_0, N, MSD)
      Implicit none
      integer, intent(in) :: N
      real(8), dimension(N, 3), intent(in) :: r, r_0
      real(8), intent(out) :: MSD
      real(8), dimension(3) :: resta
      real(8) :: modul_2
      integer :: i

      MSD = 0.d0

      do i = 1, N
         resta = r(i, :) - r_0(i, :)
         modul_2 = resta(1)**2 + resta(2)**2 + resta(3)**2
         MSD = MSD + modul_2
      end do

      MSD = MSD/N

   End Subroutine

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

!#################################################################

   Subroutine BM(ndat, xnums, sigma)
      Implicit none
      Integer ::  ndat, i
      real(8), dimension(ndat) :: xnums
      real(8) :: r, phi, x1, x2, sigma
      real(8), parameter :: pi = 4.d0*atan(1.d0)
!     ATENCIÓ! Es generen 2ndat numeros
      Do i = 1, ndat, 2
         r = sqrt(-2.d0*log(1.d0 - rand()))
         phi = 2.d0*pi*rand()
         x1 = r*cos(phi)
         x2 = r*sin(phi)
         if (i .ne. ndat) then ! Ens assegurem que no haguem acabat la llista
            xnums(i) = x1*sigma
            xnums(i + 1) = x2*sigma
         end if
      end do
      return
   end Subroutine

end module
