! Analysis of properties of a Lennard-Jones liquid. Considering that the interaction between argon atoms
! can be described using a Lennard-Jones potential with epsilon = 0.998kJ/mol and σ = 3.4Å (and the atomic mass of
! argon is m = 40g/mol):
! (a) Calculate the kinetic, potential and total energy per particle of a system of (N > 100) argon atoms in a
! fluid state with densities ρ = 0.05, 0.1, 0.2, 0.4, 0.6, 0.8m/σ 3 and a fixed temperature (using a thermal bath)
! kB T = 1.2epsilon
! • Plot a graph with the kinetic, potential and total energy as a function of density. Show your results in
! units of kJ/mol for the energies and g/cm3 for the density. How do you equilibrate your system? How
! many timesteps do you use to calculate your results?
! • Visualize (part of) the equilibrium trajectory and interpret your observations in terms of the phase
! diagram of the Lennard-Jones fluid shown in Fig. 1.
! (b) Calculate the pressure of a system of (N > 100) argon atoms at kB T = 1.2 for densities ρ =
! 0.05, 0.1, 0.2, 0.4, 0.6, 0.8m/σ 3 . Plot your results in Pascal.
!(c) (extra 1) Calculate the mean square displacement over time of argon atoms in a system of N > 100 particles
! with ρ1 = 0.05m/σ 3 and ρ1 = 0.8m/σ 3 at kB T = 1.2. Extract the corresponding diffusion coefficient
! of argon atoms. Show the plot of the mean square displacement (in units of Å2 ) versus time (in units of
! picoseconds). Interpret your results in terms of the phase diagram of the Lennard-Jones fluid shown in
! Fig. 1.
! (d) (extra 2) Obtain the radial distribution function of argon at kB T = 2.0 and ρ = 0.8m/σ 3 . Show the plot
! of the radial distribution function versus distance in units of Ångstrom.

Program P_final_b

   use subroutines_md

   Implicit none
   real(8), parameter :: nu = 0.1, mass = 40, epsilon = 0.998, sigma = 3.4
!        real(8), parameter :: k_b = 1.380649e-23
   integer, parameter :: N = 125, Nsteps_ini = 10000, Nsteps_prod = 150000
   real(8), dimension(N, 3) :: r, r_ini, vel, r_out, F, r_0
   integer :: step, i, rho_index, Nsteps, particle
   real(8) :: pot, K_energy, L, cutoff, M, a, Temp, dt, absV, p, tini, tfin, rho, sigma_gaussian, MSD
   real(8) :: T_inst, P_inst, W_inst
   real(8), dimension(6) :: rho_list
   integer, allocatable :: seed(:)
   integer :: nn
   real(8) :: acc_kin, acc_pot, acc_total, P_total, MSD_total
   character(len=4) :: rho_str

   rho_list = (/0.05, 0.1, 0.2, 0.4, 0.6, 0.8/)

   dt = 1e-5

   call random_seed(size=nn)
   allocate (seed(nn))
   seed = 123456789    ! putting arbitrary seed to all elements
   call random_seed(put=seed)
   deallocate (seed)

   open (44, file="Inicilaization-Verlet.dat")
   open (55, file="thermodynamics.dat")

   do rho_index = 1, 6

      ! system inicialization
      rho = rho_list(rho_index)

      write (rho_str, '(F4.2)') rho

      open (66, file="thermodynamics_"//rho_str//".dat")

      Temp = 100
      sigma_gaussian = Temp**(1.d0/2.d0)

      write (44, *) ""
      write (44, *) ""
      write (44, *) "rho = ", rho

      write (55, *) ""
      write (55, *) ""
      write (55, *) "rho = ", rho

      L = (N/rho)**(1./3.)

      M = N**(1./3.)
      a = L/M

      print *, L, M, a

      cutoff = L/2.

      call initialize_positions(N, rho, r)

      ! Equilibration of the system
      do step = 1, Nsteps_ini
         !        do i = 1, N
         !                print*, i ,vel(i, :)
         !        end do
         call time_step_vVerlet(r, vel, pot, N, L, cutoff, dt)
         call therm_Andersen(vel, nu, sigma_gaussian, N)
         call kinetic_energy(vel, K_energy, N)
         !print*, real(step)/Nsteps
         if (mod(step, 1000) .eq. 0) then
            print *, real(step)/Nsteps_ini
            write (44, *) step, pot, K_energy, pot + K_energy
         end if
      end do

      ! SAve positions to calculate the MSD
      r_0 = r

      print *, "PRODUCTION STARTS"
      print *, "-------------------"

      Temp = 1.2
      sigma_gaussian = Temp**(1.d0/2.d0)
      print *, sigma_gaussian

      acc_kin = 0.d0
      acc_pot = 0.d0
      acc_total = 0.d0
      P_total = 0.d0

      do step = 1, Nsteps_prod
         call time_step_vVerlet(r, vel, pot, N, L, cutoff, dt)
         call therm_Andersen(vel, nu, sigma_gaussian, N)
         call kinetic_energy(vel, K_energy, N)
         call mean_sq_distance(r, r_0, N, MSD)
         T_inst = inst_temp(N, K_energy)
         W_inst = W(N, r, rho, T_inst, L, cutoff)
         P_inst = rho*T_inst + W_inst/L**3

         acc_kin = acc_kin + K_energy
         acc_pot = acc_pot + pot
         acc_total = acc_total + K_energy + pot
         P_total = P_total + P_inst
         MSD_total = MSD_total + MSD
!                        print*, MSD_total

         !print*, real(step)/Nsteps
         if (mod(step, 1001) .eq. 0) then
            ! block average
            print *, real(step)/Nsteps_prod, P_total
            ! average values
            write (66, *) rho, acc_kin/1000, acc_pot/1000, acc_total/1000, P_total/1000, MSD_total/1000
            ! inst values
            write (55, *) step, pot, K_energy, pot + K_energy, T_inst

            acc_kin = 0.d0
            acc_pot = 0.d0
            acc_total = 0.d0
            P_total = 0.d0
            MSD_total = 0.d0

         end if

      end do

      close (66)

      open (33, file="posis_"//rho_str//".xyz")

      write (33, *) N
      write (33, *) ""
      do particle = 1, N
         write (33, *) "A", r(particle, :)
      end do

      close (33)

   end do

   close (44)
   close (55)

End Program

