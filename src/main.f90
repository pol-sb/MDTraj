program main

    use initialization
    use boundary
    use integrators
    use mpi

    implicit none
    include "declaration_variables/parallel_variables.h"
    include "../input/parameter.h"

    integer::natoms
    double precision::L,rc
    integer::tt,gg,si,sj
    double precision::ti
    double precision,allocatable::r(:,:),v(:,:),F(:,:),F_root(:,:)
    double precision::ngr,pressp
    double precision::epot,ekin,ekin_paralel,temperature,deltag
    double precision::rpos,vb,nid
    double precision,allocatable,dimension(:) ::  gr, gr_main
    integer::nhis
    integer :: ii,jj,kk,M,count,seed(33)
    integer,allocatable :: interact_list(:,:),sizes(:),displs(:)
    integer :: particle_range(2),interact_range(2)
    double precision :: start, finish
    !Change units
    !double precision::mass
    double precision::density_au,time

    call MPI_INIT(ierror) ! Begin parallel execution code

    call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror) ! Find out which is the current
    ! process from the set of processes defined by the communicator MPI_COMM_WORLD
    ! (MPI shorthand for all the processors running this program). Value stored in
    ! rank variable.

    call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
    if (taskid.eq.0) then
        start_time = MPI_Wtime()
    end if

    ! -------------------------------------------------------------------------- !

    ! Random seed initializtaion
    seed(1:33) = rng_seed + taskid
    call random_seed(put=seed)

	!Change units
	!temp=temp/epsilon
    	!epsilon=epsilon*kb
    	!sigma=sigma*10.0d0**(-10.0d0)
    	!mass=molar_mass/Na/10.0d0/10.0d0/10.0d0
    	!density = density/mass*sigma*sigma*sigma
    
    
    
    !Initialization of the structure
    if (structure .eq. 1) then
        natoms = Nc*Nc*Nc
        L = (float(natoms)/density)**(1.0/3.0)
        !write(*,*) L
        allocate (r(natoms,3))
        if (taskid .eq. 0) then
            call initial_configuration_SC(Nc,L,r,sigma)
        end if
    elseif (structure .eq. 2) then
        natoms = Nc*Nc*Nc*4
        L = (float(natoms)/density)**(1.0/3.0)
        allocate (r(natoms,3))
        if (taskid .eq. 0) then
            call initial_configuration_fcc(Nc,L,r,sigma)
        end if
    elseif (structure .eq. 3) then
        natoms = Nc*Nc*Nc*8
        L = (float(natoms)/density)**(1.0/3.0)
        allocate (r(natoms,3))
        if (taskid .eq. 0) then
            call initial_configuration_diamond(Nc,L,r,sigma)
        end if
    else
        write (*,*) "Error,no structure found"
        stop
    end if
    
    !Change units
    !mass = mass*natoms
    !dt = dt*10.0d0**(-12.0d0)/(sigma*dsqrt(mass/epsilon))

    ! -------------------------------------------------------------------------- !
    !                                                         Select range of particles for each processor
    ! -------------------------------------------------------------------------- !
    blocksize = natoms/numproc
    residu = mod(natoms,numproc)
    allocate (sizes(numproc),displs(numproc))

    if (taskid .lt. residu) then
        first_particle = taskid*(blocksize + 1) + 1
        last_particle = blocksize + first_particle
    elseif (taskid .ge. residu) then
        first_particle = taskid*blocksize + 1 + residu
        last_particle = (blocksize - 1) + first_particle
    end if

    count = 0
    do ii = 1,numproc
        if (ii - 1 .lt. residu) then
            sizes(ii) = blocksize + 1
        else
            sizes(ii) = blocksize
        end if
        displs(ii) = count; count = count + sizes(ii)
    end do
    particle_range(1) = first_particle; particle_range(2) = last_particle;
    ! -------------------------------------------------------------------------- !
    !                                                         Select range of interactions for each processor
    ! -------------------------------------------------------------------------- !
    num_interacts = natoms*(natoms - 1)/2
    allocate (interact_list(num_interacts,2))
    kk = 1
    do ii = 1,natoms - 1
        do jj = ii + 1,natoms
            interact_list(kk,:) = (/ii,jj/); kk = kk + 1
        end do
    end do

    inter_blocksize = num_interacts/numproc
    inter_residu = mod(num_interacts,numproc)
    if (taskid .lt. inter_residu) then
        first_inter = taskid*(inter_blocksize + 1) + 1
        last_inter = inter_blocksize + first_inter
    elseif (taskid .ge. inter_residu) then
        first_inter = taskid*inter_blocksize + 1 + inter_residu
        last_inter = (inter_blocksize - 1) + first_inter
    end if
    interact_range(1) = first_inter; interact_range(2) = last_inter;
    ! -------------------------------------------------------------------------- !
    ! -------------------------------------------------------------------------- !

    nhis = 250; deltag = L/(2.d0*dble(nhis)); rc = L/2.d0
    allocate(gr(nhis),gr_main(nhis)); gr = 0.d0; gr_main = 0.d0
    allocate(v(last_particle-first_particle+1,3),F(natoms,3),F_root(natoms,3))

    !initialization of velocity
    if (vel_opt .eq. 1) then
        call bimodal(Temp,v)
    else
        v(:,:) = 0.0d0
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_Bcast(r,natoms*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call force(natoms,r,L,rc,F,epot,pressp,gr,deltag,interact_range,&
               interact_list)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    !call MPI_ALLGATHER(r,natoms*3,MPI_DOUBLE_PRECISION,r,natoms*3,&
    !      MPI_DOUBLE,MPI_COMM_WORLD,ierror)
    call MPI_ALLREDUCE(F,F_root,natoms*3,MPI_DOUBLE_PRECISION,MPI_SUM,&
                       MPI_COMM_WORLD,ierror)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    F = F_root

    if (taskid .eq. 0) then
        open (11,file='output/temp.dat',status='unknown')
        open (12,file='output/energy.dat',status='unknown')
        open (13,file='output/pressure.dat',status='unknown')
        open (14,file='output/trajectory.xyz',status='unknown')
    end if

    do tt = 1,ntimes,1
        ti = ti + dt ! Actualizing the instant time

        ! set g(r) = 0.d0 while the initial structure is melting (equilibrating)
        ! in order to obtain a clean plot of the rdf
        if (tt .lt. tmelt) then
            gr = 0.d0; ngr = 0
        end if

        if (thermo .eq. 0) then
            call vel_verlet(natoms,r,v,F,epot,dt,rc,L,pressp,&
                                        gr,deltag,particle_range,interact_range,interact_list,&
                                        sizes,displs)
        elseif (thermo .eq. 1) then
            call vel_verlet_with_thermo(natoms,r,v,F,epot,dt,rc,L,temp,pressp,&
                                        gr,deltag,particle_range,interact_range,interact_list,&
                                        sizes,displs)
        else
            write (*,*) "Error,no thermostat status found"
            stop
        end if

        ekin_paralel = kinetic(v,natoms,particle_range)
        call MPI_REDUCE(ekin_paralel,ekin,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
                        MPI_COMM_WORLD,ierror)

        call MPI_REDUCE(gr,gr_main,nhis,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
                        MPI_COMM_WORLD,ierror)

        if (taskid .eq. 0) then
          temperature = 2.d0*ekin/(3.d0*dble(natoms) - 3.d0)
		      ngr = ngr + 1
        end if

        do si = particle_range(1),particle_range(2)
            do sj = 1,3
                call pbc(r(si,sj),L,L/2.d0)
            end do
        end do

        if (taskid .eq. 0) then
          if (mod(tt,everyt) .eq. 0) then
            !Change units
            time=ti!*(sigma*dsqrt(mass/epsilon))
            !temperature=temperature*epsilon/kb
            !ekin=ekin*epsilon
            !epot=epot*epsilon
            density_au=density!*mass/sigma/sigma/sigma
            !pressp=pressp*epsilon/sigma/sigma/sigma
            
            write(11,*) time, temperature
            write(12,*) time, epot/dble(natoms),ekin/dble(natoms),&
                (epot + ekin)/dble(natoms)
            write(13,*) time, pressp/dble(natoms),density_au*temperature/dble(natoms),&
                (pressp + density_au*temperature)/dble(natoms)
            write(14,*) natoms
            write(14,*)
            do si = 1,natoms
            	!Change units (add sigma (m))
                write(14,*) 'He',(r(si,sj)*sigma,sj=1,3)
            end do
          end if
        end if

    end do

    if (taskid .eq. 0) then
        open (15,file='output/rdf.dat',status='unknown')
        do ii = 1,nhis
            rpos = deltag*(ii + 0.5) ! Distance r
            vb = ((ii + 1)**3 - ii**3)*(deltag**3)
            ! Volume between bin i+1 and i
            nid = (4.d0*PI/3.d0)*vb*density
            ! Number of ideal gas part . in vb
            gr_main(ii) = gr_main(ii)/(dble(ngr)*dble(natoms)*nid) ! Normalize g(r)
            write (15,*) rpos*sigma,gr_main(ii)
        end do
        close(11); close(12); close(13); close(14); close(15)
    end if

    deallocate(r,v,F,F_root,gr,interact_list)
    if (taskid.eq.0) then
            finish_time = MPI_Wtime()
            open(16,file='output/performance.dat',access="append",status='old')
            write(16,*) natoms, numproc , finish_time-start_time
            close(16)
    end if

	call MPI_FINALIZE(ierror) ! End parallel execution

end program main
