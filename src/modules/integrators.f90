!==============================================================================!
!                               MODULE INTEGRATORS
! This module includes all the principal integration algorithms used in the
! molecular dynamics code, as the principal algorithm and to check that it
! works properly.
!
! This module needs the parameters module (params.f90) and tools (tools.f90)
! in order to work properly
!==============================================================================!
module integrators
  use forces
	use thermostat
	 implicit none
	  
   contains

    !=========================================================================!
	!                       				EULER ALGORITHM													 !
	!=========================================================================!
   subroutine euler(natoms,r,vel,F,dt,boxlength)
	integer,intent(in)::natoms
   double precision, allocatable, dimension(:,:), intent(in) :: F
   double precision, allocatable, dimension(:,:), intent(inout) :: r, vel
   double precision, intent(in) :: dt, boxlength
   
   		 ! Input variables: r - array which contains the positions of the atoms
		 !											in the lattice
		 !					vel- array which contains the velocity of the atoms
		 !											in the lattice.
		 !					boxlength - length of the box of simulation
		 !             
		 ! 		             F - array of the interacting forces on each atom of the
		 !											simulation box.
		 !					

   do jj = 1,natoms
     do ii = 1,3
       r(jj,ii) = r(jj,ii) + vel(jj,ii)*dt + 0.5d0*F(jj,ii)*dt*dt
       vel(jj,ii) = vel(jj,ii) + F(jj,ii)*dt
     end do
   end do

   end subroutine euler

    !=========================================================================!
	!                       VELOCITY VERLET ALGORITHM													 !
	!=========================================================================!
   subroutine vel_verlet(natoms,r,vel,F,Upot,dt,rc,boxlength,pressp,gr,deltag)
     implicit none
     integer,intent(in)::natoms
     double precision, allocatable, dimension(:,:), intent(inout) :: F
     double precision, allocatable, dimension(:,:), intent(inout) :: r, vel
     double precision, allocatable, dimension(:), intent(inout) :: gr
     double precision, intent(in) :: dt, rc, boxlength, deltag
     double precision, intent(out) :: Upot, pressp
      	 
      	 ! Input variables: r - array which contains the positions of the atoms
		 !											in the lattice
		 !									boxlength - length of the box of simulation
		 !                  rc - cutoff radius from which interactions are neglected
		 !					dt - time step
		 ! Output:          F - array of the interacting forces on each atom of the
		 !											simulation box.
		 !									Upot - potential energy term at time ti
		 ! 									gr - vector which contains the radial distribution function.
		 !                  deltag - width of the bin of the rdf
		 !									pressp - potential pressure term at time ti

     do ii = 1,natoms
       do jj = 1,3
         r(jj,ii) = r(jj,ii) + vel(jj,ii)*dt + 0.5d0*F(jj,ii)*dt*dt
         vel(jj,ii) = vel(jj,ii) + F(jj,ii)*0.5d0*dt
       end do
     end do

     call force(natoms,r,boxlength,rc,F,Upot,pressp,gr,deltag)
     do ii = 1,natoms
       do jj = 1,3
         vel(jj,ii) = vel(jj,ii) + F(jj,ii)*0.5d0*dt
       end do
     end do
   end subroutine vel_verlet

    !=========================================================================!
	!                       VELOCITY VERLET ALGORITHM													 !
	!=========================================================================!
   subroutine vel_verlet_with_thermo(natoms,r,vel,F,Upot,dt,rc,boxlength,Temp,pressp,gr,deltag)
     implicit none
     integer,intent(in)::natoms
     double precision, allocatable, dimension(:,:), intent(inout) :: F
     double precision, allocatable, dimension(:,:), intent(inout) :: r, vel
     double precision, allocatable, dimension(:), intent(inout) :: gr
     double precision, intent(in) :: dt, rc, boxlength, Temp, deltag
     double precision, intent(out) :: Upot, pressp

	! Input variables: r - array which contains the positions of the atoms
		 !											in the lattice
		 !									boxlength - length of the box of simulation
		 !                  rc - cutoff radius from which interactions are neglected
		 !					dt - time step
		 ! Output:          F - array of the interacting forces on each atom of the
		 !											simulation box.
		 !									Upot - potential energy term at time ti
		 ! 									gr - vector which contains the radial distribution function.
		 !                  deltag - width of the bin of the rdf
		 !									pressp - potential pressure term at time ti
		 
	  Upot = 0.d0; pressp = 0.d0

     do jj = 1,natoms
       do ii = 1,3
         r(jj,ii) = r(jj,ii) + vel(jj,ii)*dt + 0.5d0*F(jj,ii)*dt*dt
         vel(jj,ii) = vel(jj,ii) + F(jj,ii)*0.5d0*dt
       end do
     end do

     call force(natoms,r,boxlength,rc,F,Upot,pressp,gr,deltag)

     do jj = 1,natoms
       do ii = 1,3
         vel(jj,ii) = vel(jj,ii) + F(jj,ii)*0.5d0*dt
       end do
     end do

     call andersen_thermo(Temp,vel,natoms)
   end subroutine vel_verlet_with_thermo

endmodule integrators
	
