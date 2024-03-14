
Program P_final_a

   use subroutines_md

   Implicit none
   real(8), parameter :: mass = 1, rho = 0.7, epsilon = 1, sigma = 1
!        real(8), parameter :: k_b = 1.380649e-23
   integer, parameter :: N = 125
   real(8), dimension(N, 3) :: r, r_ini, vel, vel_ini, r_out, v_fin
   integer :: step, i, dt_index, Nsteps
   real(8) :: pot, K_energy, L, cutoff, M, a, Temp, inst_temp, dt, absV, p, tini, tfin
   real(8), dimension(3) :: dt_list
   integer, allocatable :: seed(:)
   integer :: nn

   dt_list = (/1e-3, 1e-4, 1e-5/)

   call random_seed(size=nn)
   allocate (seed(nn))
   seed = 123456789    ! putting arbitrary seed to all elements
   call random_seed(put=seed)
   deallocate (seed)

   L = (N/rho)**(1./3.)
   Temp = 100.d0

   M = N**(1./3.)
   a = L/(M)

   print *, L, M, a

   cutoff = 2.5

   ! """"
   ! ii) Initialize system and run simulation using velocity Verlet
   !
   ! """"

   ! Initialize bimodal distrubution: v_i = +- sqrt(T' / m)
   absV = (Temp/mass)**(1./2.)
   print *, absV

   open (77, file="Temperatures_euler.dat")

   call initialize_positions(N, rho, r_ini)

   call initialize_velocities(N, absV, vel_ini)

   tini = 0
   tfin = 0.1

   ! """"
   ! iii) Initialize system and run simulation using Euler
   !
   ! """"

   ! Initialize again, now to apply Euler method
   ! Apply Euler algorithm

   open (45, file="energy_euler.dat")

   do dt_index = 1, 3
      dt = dt_list(dt_index)
      write (45, *) ""
      write (45, *) ""
      write (45, *) "dt = ", dt

      r = r_ini
      vel = vel_ini

      Nsteps = int((tfin - tini)/dt)

      ! Initialize velocities and positions

      do step = 1, Nsteps
         call time_step_Euler_pbc(r, r_out, vel, N, L, cutoff, dt, pot)
         call kinetic_energy(vel, K_energy, N)
         ! Momentum p = m*v
         !p = (2*mass*K_energy)**(1./2.)
         call momentum(vel, p, N)

         write (45, *) step*dt, pot, K_energy, pot + K_energy, p
         if (mod(step, 1000) .eq. 0) then
            print *, real(step)/Nsteps
         end if

         if (real(step)/Nsteps .gt. 0.9) then
            v_fin = v_fin + vel
         end if

         r = r_out
      end do

      Temp = inst_temp(N, K_energy)
      write (77, *) "Euler", dt, Temp

   end do

   close (45)
   close (77)

   open (24, file="vel_fin_Euler.dat")
   do i = 1, N
      write (24, *) vel(i, :)/(Nsteps*0.1), (v_fin(i, 1)**2 + v_fin(i, 2)**2 + v_fin(i, 3)**2)**(1./2.)/(Nsteps*0.1)
   end do
   close (24)

End Program

