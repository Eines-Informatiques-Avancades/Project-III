Program main

   use :: integrators
   use :: forces
   use :: initialization
   use :: pbc_module

   Implicit none
   real(8) :: mass, rho, epsilon, sigma, Temp ! (LJ units, input file)
!        real(8), parameter :: k_b = 1.380649e-23
   integer, parameter :: N = 125
   real(8), dimension(N, 3) :: r, r_ini, vel, vel_ini, r_out, v_fin
   integer :: step, i, dt_index, Nsteps
   real(8) :: pot, K_energy, L, cutoff, M, a, dt, absV, p, tini, tfin, Ppot, Pressure
   real(8), dimension(3) :: dt_list
   integer, allocatable :: seed(:)
   integer :: nn, rc

   namelist /md_params/ mass, rho, epsilon, sigma, Temp, tfin

   ! Read parameters from namMD.nml

   inquire (file='namMD.nml', iostat=rc)

   !print *, "rc = ", rc

   open(unit=99, file='namMD.nml', status='old')
   read(99, nml=md_params)
   close(99)

   dt = 1e-4

   call random_seed(size=nn)
   allocate (seed(nn))
   seed = 123456789    ! putting arbitrary seed to all elements
   call random_seed(put=seed)
   deallocate (seed)

   L = (N/rho)**(1./3.)

   M = N**(1./3.)
   a = L/(M)

   print *, "L =", L, "M =", M, "a=", a

   cutoff = L/2.d0 - 0.1d0

   ! """"
   ! ii) Initialize system and run simulation using velocity Verlet
   !
   ! """"

   ! Initialize bimodal distrubution: v_i = +- sqrt(T' / m)
   absV = (Temp/mass)**(1./2.)
   print *, absV

   open (22, file="vel_ini.dat")
   open (33, file="pos_ini.dat")
   open (55, file="pos_out.dat")

   call initialize_positions(N, rho, r_ini)

   ! Write initial positions to file
   do i = 1, N
      write (33, *) r_ini(i, :)
   end do

   close (33)

   call initialize_velocities(N, absV, vel_ini)

   ! Write initial velocities to file
   do i = 1, N
      write (22, *) vel_ini(i, :)
   end do

   close (22)

   open (44, file="energy_verlet.dat")
   open (77, file="Temperatures_verlet.dat")
   open (23, file="vel_fin_verlet.dat")
   open (96, file="pressure_verlet.dat")
   ! Time parameters, initial and final time (input.txt)
   tini = 0
   

   ! Apply Verlet algorithm
   write (44, *) ""
   write (44, *) ""
   write (44, *) "dt = ", dt
   write (44, *) "#  time , pot, kin , total , momentum"
   write (77, *) "#  time , Temp"
   write (96, *) "#  time , pressure"
   print *, "dt = ", dt
   print *, "# time , pot, kin , total , momentum"

   Nsteps = int((tfin - tini)/dt)

   ! We roll back to the initial positions and velocities to initialize
   r = r_ini
   vel = vel_ini

   do step = 1, Nsteps

      call time_step_vVerlet(r, vel, pot, N, L, cutoff, dt, Ppot)
      call kinetic_energy(vel, K_energy, N)
      call momentum(vel, p, N)
      ! Calculate temperature
      Temp = inst_temp(N, K_energy)
      ! Calculate pressure
      Pressure = (2*K_energy + Ppot)/(3*L**3)
      write (96, *) step*dt, Pressure
      write (77, *) step*dt, Temp
      write (44, *) step*dt, pot, K_energy, pot + K_energy, p
      !        print*, K_energy
      if (mod(step, 1000) .eq. 0) then
         print *, int(real(step)/Nsteps*100), "%"
      end if

      ! We save the last 10% positions and velocity components of the simulation
      if (real(step)/Nsteps .gt. 0.9) then
         v_fin = v_fin + vel
         r_out = r_out + r
      end if

   end do

   ! Write final positions to file to plot the distribution of positions
   write (55, *) "#  Positions components (x, y, z) for the last 10% of the simulation"
   write (55, *) "#  x, y, z"
   do i = 1, N
      write (55, *) r_out(i, :)/(Nsteps*0.1)
   end do

   ! Write final velocities to file to plot the distribution of velocities
   write (23, *) "#  Velocities components (x, y, z) and modulus (v) for the last 10% of the simulation"
   write (23, *) "#  v_x, v_y, v_z, v"
   do i = 1, N
      write (23, *) v_fin(i, :)/(Nsteps*0.1), (v_fin(i, 1)**2 + v_fin(i, 2)**2 + v_fin(i, 3)**2)**(1./2.)/(Nsteps*0.1)
   end do

   write (23, *) ""
   write (23, *) ""

   close (55)
   close (44)
   close (23)
   close (77)
   close (96)
End Program
