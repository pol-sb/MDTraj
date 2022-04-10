!==============================================================================!
!                          MODULE Thermostat
! Module that contains all the tools needed to run a complete and succesfull
! simulation.
!==============================================================================!
! The module contains:
!       -> andersen_thermo (temp,vel,natoms):
!           Thermostat that controls the kinetic energy of the  simulation
!                   Input:
!                       - Temp (temperature)(in) : double precision scalar
!                       - vel (velocity)(inout) : double precision array
!                   Output:
!                       - vel (velocity)(inout) : double precision array
!
!
!     -> normal_rand:
!           Returns a normal distribution
!
!
!      -> kinetic(vel,natoms):
!           Function that given the velocities of the atoms computes the
!           associated kinetic energy.
!                   Input:
!                       - vel(velocity of the atoms)(in) : double precision array
!                       - natoms(number of atoms)(in) : integer scalar
!==============================================================================!
module thermostat
	implicit none


contains
!==============================================================================!
!                   ANDERSEN THERMOSTAT
!==============================================================================!
!  Input:
!      - Temp (temperature)(in) : double precision scalar
!      - vel (velocity)(inout) : double precision array
!  Output:
!      - vel (velocity)(inout) : double precision array
!  Depencency:
!      - random_number() : Intrinsic Fortran 90
!      - normal_rand()  : Tool in this module
!==============================================================================!
   subroutine andersen_thermo(temp,vel,natoms,particle_range)
   	double precision, intent(inout) :: vel(:,:)
   	double precision, intent(in) :: Temp
   	integer, intent(in) :: natoms, particle_range(2)
   	double precision :: nu, sigma
   	double precision :: x_rand(size(vel,2))
		double precision :: vel_normalrand(4)
		integer :: pp, qq
   	nu = 1e-3
   	sigma = dsqrt(Temp)
   	call random_number(x_rand)

   	do pp = particle_range(1)-particle_range(1)+1,particle_range(2)-particle_range(1)+1
     		if (x_rand(pp).lt.nu) then ! choosing if velocity of particle i gets changed
       		call normal_rand(sigma,vel_normalrand(1),vel_normalrand(2))
       		call normal_rand(sigma,vel_normalrand(3),vel_normalrand(4))
      	   ! The subroutine normal_rand returns a random number from a normal
       		! distribution with standard deviation \sigma = sqrt(T)
      	 	do qq = 1,3
      	   	vel(pp,qq) = vel_normalrand(qq)
      	 	enddo
     		endif
   	enddo

   endsubroutine  andersen_thermo


!==============================================================================!
!                        NORMAL RANDOM NUMBER GENERATOR
!==============================================================================!
!  Input:
!      - sigma(standar deviation)(in) : double precision
!
!  Output:
!      - xout1 (random nomal distributed number)(inout) : double precision array
!      - xout2 (random nomal distributed number)(inout) : double precision array
!  Depencency:
!      - random_number() : Intrinsic Fortran 90
!
!==============================================================================!
   subroutine normal_rand(sigma, xout1, xout2)
     double precision :: sigma
     double precision xout1, xout2
     double precision :: x(2)
     double precision, parameter :: PI = 4.d0*datan(1.d0)

     call random_number(x)

     xout1 = sigma*dsqrt(-2d0*(dlog(1d0-x(1))))*dcos(2d0*PI*x(2))
     xout2 = sigma*dsqrt(-2d0*(dlog(1d0-x(1))))*dsin(2d0*PI*x(2))

	endsubroutine normal_rand

!==============================================================================!
!                        kinetic function
!==============================================================================!
!  Input:
!		 - vel(velocity of the atoms)(in) : double precision array
!      - natoms(number of atoms)(in) : integer scalar
!
!==============================================================================!
   double precision function kinetic(vel,natoms,particle_range) result(ekin)
   	integer,intent(in)::natoms, particle_range(2)
   	double precision, allocatable, intent(in) :: vel(:,:)
   	integer::pp

			ekin = 0.d0
			do pp = particle_range(1)-particle_range(1)+1,particle_range(2)-particle_range(1)+1
		   	ekin = ekin + 0.5d0*(vel(pp,1)**2 + vel(pp,2)**2 + vel(pp,3)**2)
		  enddo

   endfunction kinetic


endmodule thermostat
