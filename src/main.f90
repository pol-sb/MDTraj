program main

use	initialization
use boundary
use integrators
use mpi

implicit none
include "declaration_variables/parallel_variables.h"
include "../input/parameter.h"

integer::natoms
double precision::L, rc
integer::tt,gg,si,sj
double precision::ti
double precision,allocatable::r(:,:),v(:,:),F(:,:), F_root(:,:)
double precision::ngr, pressp
double precision::epot, ekin, ekin_paralel, temperature, deltag
double precision::rpos, vb, nid
double precision, allocatable, dimension(:) ::  gr
integer::nhis
integer :: ii, jj, kk, M, count, seed(33)
integer, allocatable :: interact_list(:,:), sizes(:), displs(:)
integer :: particle_range(2), interact_range(2)


	call MPI_INIT(ierror) ! Begin parallel execution code

	call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror) ! Find out which is the current
	! process from the set of processes defined by the communicator MPI_COMM_WORLD
	! (MPI shorthand for all the processors running this program). Value stored in
	! rank variable.

	call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

	! -------------------------------------------------------------------------- !

				! Random seed initializtaion
	seed(1:33) = rng_seed+taskid
	call random_seed(put=seed)

	!Initialization of the structure
	if (structure .eq. 1) then
		natoms=Nc*Nc*Nc
		L= (float(natoms)/density)**(1.0/3.0)
		!write(*,*) L
		allocate(r(natoms,3))
		if (taskid.eq.0) then
			call initial_configuration_SC(Nc,L,r,taskid)
		end if
	elseif (structure .eq. 2) then
		natoms=Nc*Nc*Nc*4
		L= (float(natoms)/density)**(1.0/3.0)
		allocate(r(natoms,3))
		if (taskid.eq.0) then
			call initial_configuration_fcc(Nc,L,r,taskid)
		end if
	elseif (structure .eq. 3) then
		natoms=Nc*Nc*Nc*8
		L= (float(natoms)/density)**(1.0/3.0)
		allocate(r(natoms,3))
		if (taskid.eq.0) then
			call initial_configuration_diamond(Nc,L,r,taskid)
		end if
	else
		write(*,*)"Error, no structure found"
		stop
	endif

	! -------------------------------------------------------------------------- !
	! 							Select range of particles for each processor
	! -------------------------------------------------------------------------- !
	blocksize = natoms/numproc
	residu = mod(natoms,numproc)
	allocate(sizes(numproc),displs(numproc))

	if (taskid.lt.residu) then
		first_particle = taskid*(blocksize+1) + 1
		last_particle = blocksize + first_particle
	elseif (taskid.ge.residu) then
		first_particle = taskid*blocksize + 1 + residu
		last_particle = (blocksize-1) + first_particle
	end if

	count = 0
	do ii = 1,numproc
		if (ii-1.lt.residu) then
			sizes(ii) = blocksize+1
		else
			sizes(ii) = blocksize
		end if
		displs(ii) =  count; count = count+sizes(ii)
	end do
	particle_range(1) = first_particle; particle_range(2) = last_particle

  ! -------------------------------------------------------------------------- !
  ! -------------------------------------------------------------------------- !

	nhis = 250; deltag = L/(2.d0*dble(nhis)); rc = L/2.d0
  allocate(gr(nhis)); gr = 0.d0; ngr = 0
	allocate(v(last_particle-first_particle+1,3))
	allocate(F(last_particle-first_particle+1,3))

	!initialization of velocity
	if (vel_opt .eq. 1) then
		call bimodal(Temp,v)
	else
		v(:,:)=0.0d0
	endif

	call MPI_BARRIER(MPI_COMM_WORLD, ierror)

	call MPI_Bcast(r,natoms*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)

	call force(natoms,r,L,rc,F,particle_range)

	call MPI_BARRIER(MPI_COMM_WORLD, ierror)

	if (taskid.eq.0) then
		open(11,file='output/temp.dat',status='unknown')
		open(12,file='output/energy.dat',status='unknown')
		open(13,file='output/pressure.dat',status='unknown')
		open(14,file='output/trajectory.xyz',status='unknown')
	end if

	do tt = 1,ntimes,1
		ti = ti+dt ! Actualizing the instant time

		! set g(r) = 0.d0 while the initial structure is melting (equilibrating)
		! in order to obtain a clean plot of the rdf
		if (tt.lt.tmelt) then
			gr = 0.d0; ngr = 0
		end if

		if (thermo.eq.0) then
			call vel_verlet(natoms,r,v,F,epot,dt,rc,L,pressp,gr,deltag,&
	      particle_range,taskid)
		elseif (thermo.eq.1) then
			call vel_verlet_with_thermo(natoms,r,v,F,epot,dt,rc,L,Temp,pressp,&
		    gr,deltag,particle_range,sizes,displs,taskid)
		else
			write(*,*) "Error, no thermostat status found"
			stop
		endif

		ekin_paralel = kinetic(v,natoms,particle_range)
		call MPI_REDUCE(ekin_paralel,ekin,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
                  									MPI_COMM_WORLD,ierror)
		if (taskid.eq.0) then
			temperature = 2.d0*ekin/(3.d0*dble(natoms)-3.d0)
			ngr = ngr+1
		end if

		do si = first_particle,last_particle
			do sj = 1,3
				call pbc(r(si,sj),L,L/2.d0)
			end do
		end do

		if (taskid.eq.0) then
			if (mod(tt,everyt).eq.0) then
				write(11,*) ti, temperature
				write(12,*) ti, epot/dble(natoms), ekin/dble(natoms), &
				(epot+ekin)/dble(natoms)
				write(13,*) ti, pressp/dble(natoms), density*temperature/dble(natoms), &
			  	(pressp+density*temperature)/dble(natoms)
				write(14,*) natoms
				write(14,*)
				do si = 1,natoms
					write(14,*) 'He', (r(si,sj), sj=1,3)
				end do
			end if
		end if

  end do

	if (taskid.eq.0) then
		close(11); close(12); close(13); close(14)
	end if

	if (taskid.eq.0) then
		open(15,file='output/rdf.dat',status='unknown')
		do ii= 1,nhis
			rpos = deltag*(ii + 0.5) ! Distance r
			vb = ((ii+1)**3 - ii**3)*(deltag**3)
			! Volume between bin i+1 and i
			nid = (4.d0*PI/3.d0)*vb*density
			! Number of ideal gas part . in vb
			gr(ii) = gr(ii)/(dble(ngr)*dble(natoms)*nid) ! Normalize g(r)
			write(15,*) rpos*sigma, gr(ii)
		enddo
		close(15)
	end if

	deallocate(r,F,v,gr,sizes,displs)
	call MPI_BARRIER(MPI_COMM_WORLD, ierror)

	call MPI_FINALIZE(ierror) ! End parallel execution

endprogram main
