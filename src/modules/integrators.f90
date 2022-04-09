!=====================================================================================!
!                               MODULE INTEGRATORS
! This module includes all the principal integration algorithms used in the
! molecular dynamics code, as the principal algorithm and to check that it
! works properly.
!
! This module needs the parameters module (params.f90) and tools (tools.f90)
! in order to work properly
!=====================================================================================!
! The module contains:
!       -> euler (natoms,r,vel,F,dt,boxlength):
!           Euler Integration Method. It returns the new positions and velocities.
!                   Input:
!                     - natoms (número de átomos)(in) : integer scalar
!	   				  - r (positions)(inout) : double precision array
!      				  - F (force)(inout) : double precision array
! 	   				  -dt (time step)(in) : double precision scalar
!	   				  -boxlength (longitude of one side)(in) double precision scalar
!                   Output:
! 					   - r (positions)(inout) : double precision array
!                      - vel (velocity)(inout) : double precision array
!
!
!       -> vel_verlet (natoms,r,vel,F,Upot,dt,rc,boxlength,pressp,gr,deltag):
!           The verlet algorithm method wihout considering the temperature.
!			It returns the new positions and velocities.
!                   Input:
!                      - natoms (número de átomos)(in) : integer scalar
! 	   				   - r (positions)(inout) : double precision array
!      				   - vel (velocity)(inout) : double precision array
!      				   - F (force)(inout) : double precision array
! 	   				   -dt (time step)(in) : double precision scalar
!	   				   -rc (cut-off distance) (in) : double precision
!	   				   -boxlength (longitude of one side)(in) : double precision scalar
!	   				   -gr ()(in) : double precision array
!      				   -deltag ()(in) : double precision scalar
!                   Output:
!                     - r (positions)(inout) : double precision array
!      				    - vel (velocity)(inout) : double precision array
!	   				    - Upot (potential energy) (out) : double precision scalar
!	   				    - Pressp (pressure) (out) : double precision scalar
!						  Depencency:
!      					  - force() : In module forces (src/modules/forces.f90)
!
!
!       -> vel_verlet_with_thermo (natoms,r,vel,F,Upot,dt,rc,boxlength,Temp,pressp,gr,deltag):
!           The verlet algorithm method considering the temperature.
!			It returns the new positions and velocities.
!                   Input:
!                      - natoms (número de átomos)(in) : integer scalar
! 	   				   - r (positions)(inout) : double precision array
!      				   - vel (velocity)(inout) : double precision array
!      				   - F (force)(inout) : double precision array
! 	   				   -dt (time step)(in) : double precision scalar
!	   				   -rc (cut-off distance) (in) : double precision
!	   				   -boxlength (longitude of one side)(in) : double precision scalar
!	   				   -gr ()(in) : double precision array
!      				   -deltag ()(in) : double precision scalar
!                   Output:
!                      - r (positions)(inout) : double precision array
!      				   - vel (velocity)(inout) : double precision array
!	   				   - Upot (potential energy)(out) : double precision scalar
! 	   				   - temp (temperature)(in) : double precision scalar
!	   				   - Pressp (pressure) (out) : double precision scalar
!						  Depencency:
!      					  - force() : In module forces (src/modules/forces.f90)
! 							  - andersen_thermo() : In module thermostat (src/modules/thermostats.f90)
!=====================================================================================!
module integrators
  use forces
	use thermostat
  use mpi
	implicit none

  contains

!=====================================================================================!
!                   		  EULER INTEGRATION
!=====================================================================================!
!  Input:
!      - natoms (número de átomos)(in) : integer scalar
!	   - r (positions)(inout) : double precision array
!      - F (force)(inout) : double precision array
! 	   -dt (time step)(in) : double precision scalar
!	   -boxlength (longitude of one side)(in) double precision scalar
!  Output:
! 	   - r (positions)(inout) : double precision array
!      - vel (velocity)(inout) : double precision array
!=====================================================================================!
   subroutine euler(natoms,r,vel,F,dt,boxlength)
		integer,intent(in)::natoms
   	double precision, allocatable, intent(in) :: F(:,:)
   	double precision, allocatable, intent(inout) :: r(:,:), vel(:,:)
   	double precision, intent(in) :: dt, boxlength
    integer :: ii, jj

  		do jj = 1,natoms
     		do ii = 1,3
       		r(jj,ii) = r(jj,ii) + vel(jj,ii)*dt + 0.5d0*F(jj,ii)*dt*dt
       		vel(jj,ii) = vel(jj,ii) + F(jj,ii)*dt
     		enddo
   	enddo

   endsubroutine euler

!=====================================================================================!
!                     VELOCITY VERLET INTEGRATION (Without temperature)
!=====================================================================================!
!  Input:
!      - natoms (número de átomos)(in) : integer scalar
! 	   - r (positions)(inout) : double precision array
!      - vel (velocity)(inout) : double precision array
!      - F (force)(inout) : double precision array
! 	   -dt (time step)(in) : double precision scalar
!	   -rc (cut-off distance) (in) : double precision
!	   -boxlength (longitude of one side)(in) : double precision scalar
!	   -gr ()(in) : double precision array
!      -deltag ()(in) : double precision scalar
!  Output:
!	   - r (positions)(inout) : double precision array
!      - vel (velocity)(inout) : double precision array
!	   - Upot (potential energy) (out) : double precision scalar
!	   - Pressp (pressure) (out) : double precision scalar
!  Depencency:
!      - force() : In module forces (src/modules/forces.f90)
!=====================================================================================!
   subroutine vel_verlet(natoms,r,vel,F,Upot,dt,rc,boxlength,pressp,gr,deltag,&
     particle_range,interact_range,interact_list,taskid)
		integer,intent(in)::natoms, particle_range(2), interact_range(2), taskid
    integer, allocatable, intent(in):: interact_list(:,:)
   	double precision, allocatable,intent(inout) :: F(:,:)
   	double precision, allocatable, intent(inout) :: r(:,:), vel(:,:)
   	double precision, allocatable, intent(inout) :: gr(:)
   	double precision, intent(in) :: dt, rc, boxlength, deltag
   	double precision, intent(out) :: Upot, pressp
    integer :: ii, jj


		do jj = particle_range(1),particle_range(2)
		 	do ii = 1,3
				r(jj,ii) = r(jj,ii) + vel(jj,ii)*dt + 0.5d0*F(jj,ii)*dt*dt
				vel(jj,ii) = vel(jj,ii) + F(jj,ii)*0.5d0*dt
			enddo
		enddo

		call force(natoms,r,boxlength,rc,F,Upot,pressp,gr,deltag,interact_range,&
              interact_list,taskid)

		do jj = particle_range(1),particle_range(2)
			do ii = 1,3
				vel(jj,ii) = vel(jj,ii) + F(jj,ii)*0.5d0*dt
			enddo
		enddo

   endsubroutine vel_verlet

!=====================================================================================!
!                     VELOCITY VERLET INTEGRATION (With temperature)
!=====================================================================================!
!  Input:
!      - natoms (número de átomos)(in) : integer scalar
! 	   - r (positions)(inout) : double precision array
!      - vel (velocity)(inout) : double precision array
!      - F (force)(inout) : double precision array
! 	   -dt (time step)(in) : double precision scalar
!	   -rc (cut-off distance) (in) : double precision
!	   -boxlength (longitude of one side)(in) : double precision scalar
!	   -gr ()(in) : double precision array
!      -deltag ()(in) : double precision scalar
!  Output:
!	   - r (positions)(inout) : double precision array
!      - vel (velocity)(inout) : double precision array
!	   - Upot (potential energy)(out) : double precision scalar
! 	   - temp (temperature)(in) : double precision scalar
!	   - Pressp (pressure) (out) : double precision scalar
!  Depencency:
!      - force() : In module forces (src/modules/forces.f90)
!      - andersen_thermo() : In module thermostat (src/modules/thermostats.f90)
!=====================================================================================!
	subroutine vel_verlet_with_thermo(natoms,r,vel,F,epot,dt,rc,boxlength,Temp,pressp,&
    gr,deltag,particle_range,interact_range,interact_list,taskid)
    include "../declaration_variables/parallel_variables.h"
		integer,intent(in)::natoms, particle_range(2), interact_range(2)
    integer, allocatable, intent(in):: interact_list(:,:)
   	double precision, allocatable,  intent(inout) :: F(:,:)
   	double precision, allocatable, intent(inout) :: r(:,:), vel(:,:)
   	double precision, allocatable,  intent(inout) :: gr(:)
   	double precision, intent(in) :: dt, rc, boxlength, Temp, deltag
   	double precision, intent(out) :: pressp, epot
    double precision :: Upot
    double precision :: F_root(size(F,1),size(F,2))
    integer ii, jj

		Upot = 0.d0; pressp = 0.d0
    epot = 0.d0

    ! <------ aqui se necesitan las fuerzas repartidas entre todos los workers
		do jj = particle_range(1),particle_range(2)
			do ii = 1,3
        r(jj,ii) = r(jj,ii) + vel(jj,ii)*dt + 0.5d0*F(jj,ii)*dt*dt
				vel(jj,ii) = vel(jj,ii) + F(jj,ii)*0.5d0*dt
			enddo
		enddo
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    print*, "Before allgather"
    call MPI_ALLGATHER(r, natoms*3, MPI_DOUBLE_PRECISION, r, natoms*3,&
          MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

    print*, "after allgather"
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    print*, "============================================================="
    print*, 'taskid = ', taskid, interact_list(1,:), interact_list(1939,:)
    print*, 'taskid = ', taskid, interact_list(2,:), interact_list(1940,:)
    print*, "============================================================="
    print*, "============================================================="
    print*, "============================================================="
    if (taskid.eq.0) then
      print*, interact_list(1:2*1939,1)
    end if
    ! allgather should be applied into the the r and vel array
    ! <------- aqui se necesita haber repartido todas las posiciones
		call force(natoms,r,boxlength,rc,F,Upot,pressp,gr,deltag,interact_range,&
              interact_list,taskid)
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    print*, "before allreduce"

    ! reduce should be applied into Upot,pressp, and gr
    call MPI_ALLREDUCE(F,F_root,natoms*3,MPI_DOUBLE_PRECISION,MPI_SUM,&
  									MPI_COMM_WORLD,ierror)
    print*, "after allreduce"
    call MPI_REDUCE(Upot,epot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
                  									MPI_COMM_WORLD,ierror)
  	call MPI_BARRIER(MPI_COMM_WORLD, ierror)
  	F = F_root
    ! allgather should be applied into the the F array

		do jj = particle_range(1),particle_range(2)
			do ii = 1,3
        vel(jj,ii) = vel(jj,ii) + F(jj,ii)*0.5d0*dt
      enddo
    enddo
    ! allgather should be applied into the the r and vel array

    call andersen_thermo(Temp,vel,natoms,particle_range)
    ! allgather should be applied into the the vel array

   endsubroutine vel_verlet_with_thermo

endmodule integrators
