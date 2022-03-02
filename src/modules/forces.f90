!==============================================================================!
!                             MODULE FORCES
! This module contains all the subroutines regarding the forces that drive the
! trajectories of the molecular dynamics.
!
! The module needs the module $(params.f90) and $(boundary.f90) to work properly
! which incloude the periodic boundary conditions subroutines.
!
!==============================================================================!

module forces
!   use params
	 use boundary
   implicit none
   integer :: ii, jj, kk, M

   contains
		 !=========================================================================!
		 !                       				FORCES SUBROUTINE													 !
		 !=========================================================================!
		 subroutine forces(r,boxlength,rc,F,epot,press,gr,deltag)
		   double precision, allocatable, dimension(:,:), intent(in) :: r
		   double precision, allocatable, dimension(:,:), intent(inout) :: F
		   double precision, allocatable, dimension(:), intent(inout) :: gr
		   double precision, intent(in) :: boxlength, rc, deltag
		   double precision, intent(out) :: epot, press
		   double precision :: vol, rho, factp, facte
			 double precision :: cutoff_press, cutoff_pot, pot, piter
		   integer :: ig
		 ! Input variables: r - array which contains the positions of the atoms
		 !											in the lattice
		 !									boxlength - length of the box of simulation
		 !                  rc - cutoff radius from which interactions are neglected
		 ! Output:          F - array of the interacting forces on each atom of the
		 !											simulation box.
		 !									pot - potential energy term at time ti
		 ! 									gr - vector which contains the radial distribution function.
		 !                  deltag - width of the bin of the rdf
		 !									press - potential pressure term at time ti

			 vol = boxlength**3.; rho = dble(natoms)/vol
			 facte = (8.d0/3.d0)*pi*dfloat(natoms)*rho
			 factp = (16.d0/3.d0)*pi*(rho**2)

			 cutoff_pot = 4.d0*(1.d0/(rc**12) - 1.d0/(rc**6))
			 cutoff_press = factp*((2.d0/3.d0)/(rc**9.) - 1.d0/(rc**3.))
			 !cutoff_pot = facte*((1.d0/3.d0)/(rc**9.) - 1.d0/(rc**3.))

			 press = 0.d0; epot = 0.d0
			 F = 0.d0

			 do ii = 1,natoms-1
				 do jj = ii+1,natoms
					 call lj(r,boxlength,rc,ii,jj,F,pot,piter,d)
					 ! calling function that computes the Lennard-Jones interaction between
					 ! pair of particles i and j
					 press = press + piter; epot = epot + pot

					 if (d.lt.rc) then ! computation of the radial distribution function
						 ! adding the each pair of interaction into the corresponding bin
		         ig = int(d/deltag)
		         gr(ig) = gr(ig) + 2
		       end if
				 end do
			 end do
			 pot = pot - cutoff_pot
			 press = (1.d0/(3.d0*vol))*press
			 press = press + cutoff_press
			 !epot = epot + etail; 			pressp = pressp + ptail
		 end subroutine forces

		 !=========================================================================!
		 !                  LENNARD-JONES INTERACTION COMPUTATION									 !
		 !=========================================================================!
		 subroutine lj(r,boxlength,rc,ii,jj,F,pot,piter,d)
			 double precision, allocatable, dimension(:,:), intent(in) :: r
		   double precision, allocatable, dimension(:,:), intent(inout) :: F
		   double precision, intent(in) :: boxlength, rc
			 integer, intent(in) :: ii, jj
		   double precision, intent(out) :: pot, piter
			 double precision :: dx, dy, dz, d, dU
			 ! Input variables: r - array which contains the positions of the atoms
			 !											in the lattice
			 !									boxlength - length of the box of simulation
			 !                  rc - cutoff radius from which interactions are neglected
			 !									ii - index of the first particle which is looped around
			 !									jj - index of pair particle of i at the actual
			 !											 interaction computation
			 ! Output:          F - array of the interacting forces on each atom of the
			 !											simulation box.
			 !									pot - potential energy interaction between i and j
			 ! 									piter - potential pressure term from interaction
			 !													between particles i and j
			 !                  d - distance between particles i and j

			 pot = 0.d0; piter = 0.d0
			 dx = r(ii,1)-r(jj,1); dy = r(ii,2)-r(jj,2); dz = r(ii,3)-r(jj,3)
			 ! Apply the boundary conditions to the particles distance
			 call pbc(dx,boxlength,0.d0)
			 call pbc(dy,boxlength,0.d0)
			 call pbc(dz,boxlength,0.d0)
			 d = (dx**2. + dy**2. + dz**2.)**0.5 ! Distance between particles i and j

			 if (d.lt.rc) then
				 dU = (48.d0/(d**14.0) - 24.d0/(d**8.0))
				 F(ii,1) = F(ii,1) + dU*dx
				 F(ii,2) = F(ii,2) + dU*dy
				 F(ii,3) = F(ii,3) + dU*dz

				 F(jj,1) = F(jj,1) - dU*dx
				 F(jj,2) = F(jj,2) - dU*dy
				 F(jj,3) = F(jj,3) - dU*dz

				 pot = pot + 4.d0*(1.d0/(d**12.) - 1.d0/(d**6.))
				 piter = piter + dU*dx; piter = piter + dU*dy; piter = piter + dU*dz
			 end if
	   end subroutine

end module forces
