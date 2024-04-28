Program main

   use :: integrators
   use :: forces
   use :: initialization
   use :: pbc_module

   Implicit none
   real(8) :: mass, rho, epsilon, sigma, Temp ! (LJ units, input file)
!        real(8), parameter :: k_b = 1.380649e-23
   integer, parameter :: N = 125
   real(8), dimension(N, 3) :: r, r_ini, vel, vel_ini, r_out, v_fin, r_new, vel_new
   integer :: step, i, dt_index, Nsteps
   real(8) :: pot, K_energy, L, cutoff, M, a, dt, absV, p, tini, tfin, Ppot, Pressure
   real(8), dimension(3) :: dt_list
   real(8) :: nu, sigma_gaussian
   integer, allocatable :: seed(:)
   integer :: nn, rc
   real(8) :: t1, t2

   ! MPI
   integer :: ierror, rank, nprocs, imin, imax
   integer, allocatable :: particle_distrib(:), counts_recv(:), displs_recv(:)

   include 'mpif.h'

   namelist /md_params/ mass, rho, epsilon, sigma, Temp, tfin

   ! Read parameters from namMD.nml

   inquire (file='namMD.nml', iostat=rc)
   open(unit=99, file='namMD.nml', status='old')
   if (rc /= 0) then
      print*, "Error opening namMD.nml"
      stop
   end if
   read(99, nml=md_params)
   close(99)

   L = (N/rho)**(1./3.)
   M = N**(1./3.)
   a = L/(M)
   nu = 1

   dt = 1e-4

   call random_seed(size=nn)
   allocate (seed(nn))
   seed = 123456789    ! putting arbitrary seed to all elements
   ! TO-DO: aixo cal que ho fagin tots els processadors? O nomes un i que envii la info?
   call random_seed(put=seed)
   deallocate (seed)

   call MPI_INIT(ierror)
   t1 = MPI_Wtime()


   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)

   allocate(particle_distrib(N))
   allocate(counts_recv(nprocs))
   allocate(displs_recv(nprocs))

   ! Send info to the other processors
   call MPI_Bcast(mass,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
   call MPI_Bcast(rho,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
   call MPI_Bcast(epsilon,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
   call MPI_Bcast(sigma,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
   call MPI_Bcast(Temp,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
   call MPI_Bcast(tfin,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
   call MPI_Bcast(L,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
   call MPI_Bcast(M,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
   call MPI_Bcast(a,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

   ! wait until the reading has finished
!   call MPI_Barrier(MPI_COMM_WORLD, ierror)

   print *, "L =", L, "M =", M, "a=", a

   cutoff = L/2.d0 - 0.1d0

   call distribute_particles(N, rank, nprocs, imin, imax)
   ! imin i imax tenen les particules limit de cada processador
   print*, "rank: ", rank, "imin: ", imin, "imax", imax, "particles", imax-imin + 1
   ! build counts_recv: non-negative integer array (of length group size) containing the number of elements that are received from each process (non-negative integer)
   counts_recv(rank) = imax - imin + 1

   call MPI_ALLGATHER(counts_recv(rank), 1, MPI_INTEGER, counts_recv, 1, MPI_INTEGER, MPI_COMM_WORLD, ierror)

   print*, counts_recv

   ! build displs_recv: integer array (of length group size). Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i (integer)
   if (rank > 0) then
      displs_recv(rank) = sum(counts_recv(1:rank))
   end if

   call MPI_ALLGATHER(displs_recv(rank), 1, MPI_INTEGER, displs_recv, 1, MPI_INTEGER, MPI_COMM_WORLD, ierror)
   print*, "displs_recv", displs_recv

   ! """"
   ! ii) Initialize system and run simulation using velocity Verlet
   !
   ! """"

   ! Initialize bimodal distrubution: v_i = +- sqrt(T' / m)
   absV = (Temp/mass)**(1./2.)
   print *, absV

   if ( rank == 0 ) then
      open (22, file="vel_ini.dat")
      open (33, file="pos_ini.dat")
   end if

   call initialize_positions(N, rho, r_ini)

   ! Write initial positions to file

   if ( rank == 0 ) then
      do i = 1, N
         write (33, *) r_ini(i, :)
      end do

      close (33)
   end if

   call initialize_velocities(N, absV, vel_ini)

   ! Write initial velocities to file

   if ( rank == 0 ) then
      do i = 1, N
         write (22, *) vel_ini(i, :)
      end do
   
      close (22)

      open (44, file="energy_verlet.dat")
      open (77, file="Temperatures_verlet.dat")
      open (96, file="pressure_verlet.dat")
   
   end if
   ! Time parameters, initial and final time (input.txt)
   tini = 0
   
   ! Apply Verlet algorithm

   if ( rank == 0 ) then
      write (44, *) ""
      write (44, *) ""
      write (44, *) "dt = ", dt
      write (44, *) "#  time , pot, kin , total , momentum"
      write (77, *) "#  time , Temp"
      write (96, *) "#  time , pressure"
      print *, "dt = ", dt
      print *, "# time , pot, kin , total , momentum"
   end if 

   Nsteps = int((tfin - tini)/dt)

   ! We roll back to the initial positions and velocities to initialize
   r = r_ini
   vel = vel_ini

   print*, "Loop starts"
   do step = 1, Nsteps

      call time_step_vVerlet(r, vel, pot, N, L, cutoff, dt, Ppot, imin, imax, counts_recv, displs_recv, rank, nprocs)

      call MPI_Barrier(MPI_COMM_WORLD, ierror)

      ! sincronitzar els processadors (allgather?)
      ! int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int counts_recv[], const int displs_recv[], MPI_Datatype recvtype, MPI_Comm comm)
      do i = 1, 3
         ! loop over dimensions
    !     print*, "i", i, "step", step, "rank", rank
    !     print*, imin, imax

         call MPI_ALLGATHERV(r(imin:imax, i), int(imax - imin + 1), MPI_DOUBLE_PRECISION, r_new(:,i), counts_recv, &
            displs_recv, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

         call MPI_ALLGATHERV(vel(imin:imax, i), int(imax - imin + 1), MPI_DOUBLE_PRECISION, vel_new(:,i), counts_recv,  &
            displs_recv, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)                

      end do

      call MPI_Barrier(MPI_COMM_WORLD, ierror)    
      
     ! print*, "r_new", size(r_new, 1)

      r = r_new

      vel = vel_new
!     print*, "vel_new(6)_1", vel_new(6,:)
      
      if (rank .eq. 0) then
!         print*, "vel", vel
!         call therm_Andersen(vel, nu, sigma_gaussian, N)
!         print*, "vel", vel
         call kinetic_energy(vel, K_energy, N)
!         print*, "K_energy", vel, K_energy
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

      end if 
      
      call MPI_Bcast(vel, N*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
!      print*, "step", step

   end do

   if (rank .eq. 0) then

      close (77)
      close (96)
      close (44)

      open (23, file="vel_fin_verlet.dat")
      open (55, file="pos_out.dat")
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
      close (23)

   end if

      t2 = MPI_Wtime()
      print *, "Time elapsed: ", t2 - t1
   
      call MPI_Barrier(MPI_COMM_WORLD, ierror)
   
      print*, "Barrier passed"
      call MPI_FINALIZE(ierror)
   
      print*, "MPI finalized"


End Program
