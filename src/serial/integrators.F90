
module integrators

   implicit none

contains


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
