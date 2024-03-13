
Program P_final_a

   use subroutines_md
   Implicit none
   real(8), parameter :: mass = 1, rho = 0.7, epsilon = 1, sigma = 1
!        real(8), parameter :: k_b = 1.380649e-23
   integer, parameter :: N = 125
   real(8), dimension(N, 3) :: r, r_ini, vel, vel_ini, r_out, v_fin
   integer :: step, i, dt_index, Nsteps
   real(8) :: pot, K_energy, L, cutoff, M, a, Temp, dt, absV, p, tini, tfin
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

   open (22, file="vel_ini.dat")

   call initialize_positions(N, rho, r_ini)

   call initialize_velocities(N, absV, vel_ini)

   do i = 1, N
      write (22, *) vel_ini(i, :)
   end do

   close (22)

   open (44, file="energy_verlet.dat")
   open (77, file="Temperatures_verlet.dat")
   open (23, file="vel_fin_Verlet.dat")

   tini = 0
   tfin = 1

   ! Apply Verlet algorithm
   do dt_index = 1, 3
      dt = dt_list(dt_index)
      write (44, *) ""
      write (44, *) ""
      write (44, *) "dt = ", dt
      print *, "dt = ", dt

      Nsteps = int((tfin - tini)/dt)

      ! We roll back to the initial positions and velocities to initialize
      r = r_ini
      vel = vel_ini

      do step = 1, Nsteps

         call time_step_vVerlet(r, vel, pot, N, L, cutoff, dt)
         call kinetic_energy(vel, K_energy, N)
         call momentum(vel, p, N)
         write (44, *) step*dt, pot, K_energy, pot + K_energy, p
         !        print*, K_energy
         if (mod(step, 1000) .eq. 0) then
            print *, real(step)/Nsteps
         end if

         if (real(step)/Nsteps .gt. 0.9) then
            v_fin = v_fin + vel
         end if

      end do

      Temp = inst_temp(N, K_energy)
      write (77, *) "Verlet", dt, Temp

      do i = 1, N
         write (23, *) v_fin(i, :)/(Nsteps*0.1), (v_fin(i, 1)**2 + v_fin(i, 2)**2 + v_fin(i, 3)**2)**(1./2.)/(Nsteps*0.1)
      end do

      write (23, *) ""
      write (23, *) ""

   end do

   close (44)
   close (23)
   close (77)

End Program
