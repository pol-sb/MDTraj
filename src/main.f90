program main

    use initialization
    use boundary
    use integrators
    use mpi

    implicit none
    include "declaration_variables/parallel_variables.h"
    include "../input/parameter.h"
    include "/modules/constants.h"

    integer::natoms
    double precision::L, rc
    integer::tt, gg, si, sj
    double precision::ti
    double precision, allocatable::r(:, :), v(:, :), F(:, :), F_root(:, :)
    double precision::ngr, pressp
    double precision::epot, ekin, ekin_paralel, temperature, deltag
    double precision::rpos, vb, nid
    double precision, allocatable, dimension(:) ::  gr, gr_main
    integer::nhis
    integer :: ii, jj, kk, M, count, seed(33)
    integer, allocatable :: interact_list(:, :), sizes(:), displs(:)
    integer :: particle_range(2), interact_range(2)
    double precision :: start_time, finish_time

    ! Unit conversion variables
    double precision :: density_au, time, time_fact, epsLJ, temp_fact, press_fact

    ! Begin parallel execution code
    call MPI_INIT(ierror)

    ! Find out which is the current process from the set of processes defined by
    ! the communicator MPI_COMM_WORLD (MPI shorthand for all the processors
    ! running this program). Value stored in 'rank' variable.
    call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierror)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)

    ! Starting time measurament
    if (taskid .eq. 0) then
        start_time = MPI_Wtime()
    end if

    ! Random seed initializtaion
    seed(1:33) = rng_seed + taskid
    call random_seed(put=seed)

    ! Initialization of the system structure
    if (structure .eq. 1) then
        natoms = Nc*Nc*Nc
				call reduced(taskid,epsilon,sigma,temp,density,natoms,dt,ntimes,thermo,&
		 		 									epsLJ,time_fact,press_fact,temp_fact)

        L = (float(natoms)/density)**(1.0/3.0)
        allocate (r(natoms, 3))
        if (taskid .eq. 0) then
            call initial_configuration_SC(Nc, L, r, sigma)
        end if

    elseif (structure .eq. 2) then
        natoms = Nc*Nc*Nc*4
				call reduced(taskid,epsilon,sigma,temp,density,natoms,dt,ntimes,thermo,&
		 		 									epsLJ,time_fact,press_fact,temp_fact)
        L = (float(natoms)/density)**(1.0/3.0)
        allocate (r(natoms, 3))
        if (taskid .eq. 0) then
            call initial_configuration_fcc(Nc, L, r, sigma)
        end if
    elseif (structure .eq. 3) then
        natoms = Nc*Nc*Nc*8
				call reduced(taskid,epsilon,sigma,temp,density,natoms,dt,ntimes,thermo,&
		 		 									epsLJ,time_fact,press_fact,temp_fact)
        L = (float(natoms)/density)**(1.0/3.0)
        allocate (r(natoms, 3))
        if (taskid .eq. 0) then
            call initial_configuration_diamond(Nc, L, r, sigma)
        end if
    else
        write (*, *) "Input Error: no structure found. Please input a valid structure."
        stop
    end if

    ! -------------------------------------------------------------------------- !
    !  Selecting range of particles for each processor
    ! -------------------------------------------------------------------------- !
    blocksize = natoms/numproc
    residu = mod(natoms, numproc)
    allocate (sizes(numproc), displs(numproc))

    if (taskid .lt. residu) then
        first_particle = taskid*(blocksize + 1) + 1
        last_particle = blocksize + first_particle
    elseif (taskid .ge. residu) then
        first_particle = taskid*blocksize + 1 + residu
        last_particle = (blocksize - 1) + first_particle
    end if

    count = 0
    do ii = 1, numproc
        if (ii - 1 .lt. residu) then
            sizes(ii) = blocksize + 1
        else
            sizes(ii) = blocksize
        end if
        displs(ii) = count; count = count + sizes(ii)
    end do
    particle_range(1) = first_particle; particle_range(2) = last_particle;

    ! -------------------------------------------------------------------------- !
    !    Selecting range of interactions for each processor                      !
    ! -------------------------------------------------------------------------- !
    num_interacts = natoms*(natoms - 1)/2
    allocate (interact_list(num_interacts, 2))
    kk = 1
    do ii = 1, natoms - 1
        do jj = ii + 1, natoms
            interact_list(kk, :) = (/ii, jj/); kk = kk + 1
        end do
    end do

    inter_blocksize = num_interacts/numproc
    inter_residu = mod(num_interacts, numproc)
    if (taskid .lt. inter_residu) then
        first_inter = taskid*(inter_blocksize + 1) + 1
        last_inter = inter_blocksize + first_inter
    elseif (taskid .ge. inter_residu) then
        first_inter = taskid*inter_blocksize + 1 + inter_residu
        last_inter = (inter_blocksize - 1) + first_inter
    end if
    interact_range(1) = first_inter; interact_range(2) = last_inter;
    ! -------------------------------------------------------------------------!

    nhis = 200; deltag = L/(2.d0*dble(nhis)); rc = L/2.d0
    allocate (gr(nhis), gr_main(nhis)); gr = 0.d0; gr_main = 0.d0
    allocate (v(last_particle - first_particle + 1, 3), F(natoms, 3), F_root(natoms, 3))

    !initialization of velocity
    if (vel_opt .eq. 1) then
        call bimodal(Temp, v)
    else
        v(:, :) = 0.0d0
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    ! Sending positions to all other processors
    call MPI_Bcast(r, natoms*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

    ! Computing the forces
    call force(natoms, r, L, rc, F, epot, pressp, gr, deltag, interact_range, &
               interact_list)

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    ! Adding all of the forces
    call MPI_ALLREDUCE(F, F_root, natoms*3, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierror)

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    F = F_root

    ! Opening the files where results will be written.
    if (taskid .eq. 0) then
        open (11, file='output/temp.dat', status='unknown')
        open (12, file='output/energy.dat', status='unknown')
        open (13, file='output/pressure.dat', status='unknown')
        open (14, file='output/trajectory.xyz', status='unknown')

        print *, ''
        print *, '[!] Starting simulation.'
        print *, ''
        print *, '· Melting the system from the initial configuration...'
    end if

    do tt = 1, ntimes, 1
        ! Updating the instant time
        ti = ti + dt

        ! set g(r) = 0.d0 while the initial structure is melting (equilibrating)
        ! in order to obtain a clean plot of the rdf
        if (tt .lt. tmelt) then
            gr = 0.d0; ngr = 0
        elseif (tt .eq. tmelt) then
            if (taskid .eq. 0) then
                print *, '· Computing the dynamics...'
            end if
        end if

        if (thermo .eq. 0) then
            call vel_verlet(natoms, r, v, F, epot, dt, rc, L, pressp, &
                            gr, deltag, particle_range, interact_range, interact_list, &
                            sizes, displs)
        elseif (thermo .eq. 1) then
            call vel_verlet_with_thermo(natoms, r, v, F, epot, dt, rc, L, temp, pressp, &
                                        gr, deltag, particle_range, interact_range, interact_list, &
                                        sizes, displs)
        else
            write (*, *) "Error, no thermostat status found. Please input a thermostat."
            stop
        end if

        ekin_paralel = kinetic(v, natoms, particle_range)
        ! Adding the kinetic energies
        call MPI_REDUCE(ekin_paralel, ekin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                        MPI_COMM_WORLD, ierror)

        ! Adding the gr values
        call MPI_REDUCE(gr, gr_main, nhis, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                        MPI_COMM_WORLD, ierror)

        ! Computing the temperature from the kinetic energy
        if (taskid .eq. 0) then
            temperature = 2.d0*ekin/(3.d0*dble(natoms) - 3.d0)
            ngr = ngr + 1
        end if

        ! Applying periodic boundary conditions.
        do si = particle_range(1), particle_range(2)
            do sj = 1, 3
                call pbc(r(si, sj), L, L/2.d0)
            end do
        end do

        if (taskid .eq. 0) then
            if (mod(tt, everyt) .eq. 0) then

                ! Unit conversion
                ! Converting time from reduced units to ps
                time = ti*time_fact
                ! Conversion factors to get kcal/mol
                ekin = ekin*epsilon*boltzmann_constant*avogadro_number/1e3
                epot = epot*epsilon*boltzmann_constant*avogadro_number/1e3
                ! Conversion factors to get kg/m^3
                density_au = density*atomic_mass*dble(natoms)/(avogadro_number*1e-4*(sigma**3.d0))
                ! Conversion to MPa
								!print*, 'press_fact', press_fact
								!print*, density*temperature, pressp
                pressp = pressp*press_fact

                ! Saving the results in the files.
                write (11, *) time, temperature*temp_fact
                write (12, *) time, epot/dble(natoms), ekin/dble(natoms), &
                    (epot + ekin)/dble(natoms)
                write (13, *) time, pressp/dble(natoms), density*temperature*press_fact/dble(natoms), &
                    (pressp + density*temperature*press_fact)/dble(natoms)
                write (14, *) natoms
                write (14, *)
                do si = 1, natoms
                    write (14, *) 'He', (r(si, sj), sj=1, 3)
                end do
            end if
        end if

    end do

    ! Finishing and storing the rdf results
    if (taskid .eq. 0) then
        open (15, file='output/rdf.dat', status='unknown')
        do ii = 1, nhis
            rpos = deltag*(ii + 0.5) ! Distance r
            vb = ((ii + 1)**3 - ii**3)*(deltag**3)
            ! Volume between bin i+1 and i
            nid = (4.d0*PI/3.d0)*vb*density
            ! Number of ideal gas part in vb
            gr_main(ii) = gr_main(ii)/(dble(ngr)*dble(natoms)*nid) ! Normalize g(r)
            write (15, *) rpos*sigma, gr_main(ii)
        end do

        ! Closing files
        close (11); close (12); close (13); close (14); close (15)

        print *, ''
        print *, '[!] Simulation Done!'
        print *, ''

    end if

    deallocate (r, v, F, F_root, gr, interact_list)

    ! Saving computation time results
    if (taskid .eq. 0) then
        finish_time = MPI_Wtime()
        open (16, file='output/performance.dat', access="append", status='unknown')
        write (16, *) natoms, numproc, finish_time - start_time
        close (16)
    end if

    ! Ending parallel execution
    call MPI_FINALIZE(ierror)

end program main
