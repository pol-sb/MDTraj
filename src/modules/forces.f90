!==============================================================================!
!                             MODULE FORCES
! This module contains all the subroutines regarding the forces that drive the
! trajectories of the molecular dynamics.
!==============================================================================!
! The module contains:
!                        -> force()
!            Calculates the forces between neighboor particles (using a cutoff)
!             and returns the forces, the potential energy, the pressure
!                                         Input:
!                                                        - natom (number of atoms) (in): integer scalar
!                                                        - r(coordinater)(inout): double precision array
!                                                        - boxlength (in) : double precision scalar
!                                             - rc(cutoff radius)(in):  double precision scalar
!                                                        - F (forces) (inout): double precision array
!                                               - deltag (in): double precision scalar
!                                                        - gr (g(r)) (out): double precision array
!                                         Output:
!                                                        - epot (E potential) (out): double precision scalar
!                                                        - press (pressure) (out): double precision scalar
!                                         Depencency:
!                                                        - lj() : Tool in thi module
!
!                        ->lj()
!                                         Input variables:
!                                                        -r(coordinater)(inout): double precision array
!                                                        -boxlength (in) : double precision scalar
!                                              -rc(cutoff radius)(in):  double precision scalar
!                                                        -ii (in): integer scalar
!                                                        -jj (in): integer scalar
!                                                        -F (forces) (inout): double precision array
!                                         Output:
!                                                        -pot (potential energy) (out): double precision array
!                                                         -piter (potential pressure) (out):
!                                             -d  (distance) (in)
!                                         Depencency:
!                                             -pbc(): Tool in boundary module
!
! Dependency:
!                        -> boundary.f90 module
!                        -> constants.h file
!==============================================================================!

module forces
    use boundary
    implicit none

    include "constants.h"
    !include "../../input/parameter"

contains
!==============================================================================!
!                                                       FORCES SUBROUTINE
!==============================================================================!
! Input:
!                - natom (number of atoms) (in): integer scalar
!                                        Total numbe of atoms in the system.
!                - r(coordinater)(inout): double precision array
!                              Array which contains the positions of the atoms in the lattice
!
!                - boxlength (in) : double precision scalar
!
!     - rc(cutoff radius)(in):  double precision scalar
!                                        Cutoff radius from which interactions are neglected
!
!                - F (forces) (inout): double precision array
!                                        Array of the interacting forces on each atom of the simulation box.
!
!     - deltag (in): double precision scalar
!                                        Width of the bin of the rdf
!
!                - gr (g(r)) (out): double precision array
!                                        Radial pair distribution of particles
! Output:
!                - epot (E potential) (out): double precision scalar
!                                        Potential energy term  calculated at each time step
!
!                - press (pressure) (out): double precision scalar
!                                        Potential pressure term  calculated at each time step
! Depencency:
!                - lj() : Tool in thi module
!
!==============================================================================!
    subroutine force(natoms, r, boxlength, rc, F, epot, press, gr, deltag)
        integer, intent(in)::natoms
        double precision, allocatable, intent(in) :: r(:, :)
        double precision, allocatable, intent(inout) :: F(:, :)
        double precision, allocatable, intent(inout) :: gr(:)
        double precision:: d
        double precision, intent(in) :: boxlength, rc, deltag
        double precision, intent(out) :: epot, press
        double precision :: vol, rho, factp, facte
        double precision :: cutoff_press, cutoff_pot, pot, piter
        integer :: ig

        vol = boxlength**3.; rho = dble(natoms)/vol
        facte = (8.d0/3.d0)*pi*dfloat(natoms)*rho
        factp = (16.d0/3.d0)*pi*(rho**2)

        cutoff_pot = 4.d0*(1.d0/(rc**12) - 1.d0/(rc**6))
        cutoff_press = factp*((2.d0/3.d0)/(rc**9.) - 1.d0/(rc**3.))
        !cutoff_pot = facte*((1.d0/3.d0)/(rc**9.) - 1.d0/(rc**3.))

        press = 0.d0; epot = 0.d0
        F = 0.d0

        do ii = 1, natoms - 1
            do jj = ii + 1, natoms
                call lj(r, boxlength, rc, ii, jj, F, pot, piter, d)
                ! calling function that computes the Lennard-Jones interaction between
                ! pair of particles i and j
                press = press + piter; epot = epot + pot

                if (d .lt. rc) then ! computation of the radial distribution function
                    ! adding the each pair of interaction into the corresponding bin
                    ig = int(d/deltag)
                    gr(ig) = gr(ig) + 2
                end if
            end do
        end do

        !pot = pot - cutoff_pot
        press = (1.d0/(3.d0*vol))*press
        !press = press + cutoff_press
        !epot = epot + etail;                         pressp = pressp + ptail

    end subroutine force

    subroutine force_dpd(natoms, r, boxlength, rc, F, epot, press, dt, v)
        integer, intent(in)::natoms
        double precision, allocatable, intent(in) :: r(:, :), v(:, :)
        double precision, allocatable, intent(inout) :: F(:, :)
        ! double precision, allocatable, intent(inout) :: gr(:)
        double precision:: d
        double precision, intent(in) :: boxlength, rc, dt !deltag
        double precision, intent(out) :: epot, press
        double precision :: vol, rho, factp, facte, dx, dy, dz, U
        double precision :: cutoff_press, cutoff_pot, pot, piter
        double precision :: dv, dvx, dvy, dvz, r1, r2, r3, r4
        double precision :: ig
		
        vol = boxlength**3.; rho = dble(natoms)/vol
        facte = (8.d0/3.d0)*pi*dfloat(natoms)*rho
        factp = (16.d0/3.d0)*pi*(rho**2)

        cutoff_pot = 4.d0*(1.d0/(rc**12) - 1.d0/(rc**6))
        cutoff_press = factp*((2.d0/3.d0)/(rc**9.) - 1.d0/(rc**3.))
        !cutoff_pot = facte*((1.d0/3.d0)/(rc**9.) - 1.d0/(rc**3.))

        press = 0.d0
		epot = 0.d0
        F = 0.d0

        do ii = 1, natoms - 1
            do jj = ii + 1, natoms

                !Distance between two particles
                dx = r(ii, 1) - r(jj, 1)
                dy = r(ii, 2) - r(jj, 2)
                dz = r(ii, 3) - r(jj, 3)

                !If the distance is bigger than the box -> minimum image convention:
                !Move the particle inside the box with the particle coming from the other side
                call pbc(dx, boxlength, 0.d0)
                call pbc(dy, boxlength, 0.d0)
                call pbc(dz, boxlength, 0.d0)

                ! Calculate distance module
                d = (dx**2 + dy**2 + dz**2)**0.5

                ! Diference in velocities:
                dvx = v(ii, 1) - v(jj, 1)
                dvy = v(ii, 2) - v(jj, 2)
                dvz = v(ii, 3) - v(jj, 3)

                if (d .lt. rc) then

                    ! adding the each pair of interaction into the corresponding bin
                    dv = (dx*dvx + dy*dvy + dz*dvz)/d

                    ! Random number generation.
                    call normal_rand(1d0, r1, r2)
                    call normal_rand(1d0, r3, r4)

                    !F(ii, 1) = F(ii, 1)
                    ! Conservative force
                    F(ii, 1) = dpd_aij*(1.0 - d)*dx/d
                    F(ii, 2) = dpd_aij*(1.0 - d)*dy/d
                    F(ii, 3) = dpd_aij*(1.0 - d)*dz/d

                    ! Dissipative forces
                    F(ii, 1) = F(ii, 1) - dpd_gamma*(1.0 - d)*(1.0 - d)*dv*dx/d
                    F(ii, 2) = F(ii, 2) - dpd_gamma*(1.0 - d)*(1.0 - d)*dv*dy/d
                    F(ii, 3) = F(ii, 3) - dpd_gamma*(1.0 - d)*(1.0 - d)*dv*dz/d

                    ! Random noise
                    F(ii, 1) = F(ii, 1) + (1.0 - d)*dpd_sigma*r1/sqrt(dt)*dx/d
                    F(ii, 2) = F(ii, 2) + (1.0 - d)*dpd_sigma*r2/sqrt(dt)*dy/d
                    F(ii, 3) = F(ii, 3) + (1.0 - d)*dpd_sigma*r3/sqrt(dt)*dz/d

					! potential energy computation
					epot = epot + 4.*((1./d**12)-(1./d**6))-4.*((1./(rc)**12)-(1./(rc)**6))

					!piter = piter + dU*dx
					!piter = piter + dU*dy;
					!piter = piter + dU*dz

					! computation of the radial distribution function
                    ! ig = int(d/deltag)
                    ! gr(ig) = gr(ig) + 2
                end if
            end do

		!print*,'epot:', epot

        end do

        epot = epot - cutoff_pot
        press = (1.d0/(3.d0*vol))*press
        !press = press + cutoff_press
        !epot = epot + etail;                         pressp = pressp + ptail

    end subroutine force_dpd

!==============================================================================!
!                  LENNARD-JONES INTERACTION COMPUTATION
!==============================================================================!
! Input variables:
!                -r(coordinater)(inout): double precision array
!                              Array which contains the positions of the atoms in the lattice
!                -boxlength (in) : double precision scalar
!                                        Length of the box of simulation
!     -rc(cutoff radius)(in):  double precision scalar
!                                        Cutoff radius from which interactions are neglected
!                -ii (in): integer scalar
!                                        Index of the first particle which is looped around
!                -jj (in): integer scalar
!                                        Index of pair particle of i at the actual interaction computation
!                -F (forces) (inout): double precision array
!                                        Array of the interacting forces on each atom of the simulation box.
!
! Output:
!                -pot (potential energy) (out): double precision array
!                                        Potential energy interaction between i and j elements
!                 -piter (potential pressure) (out):
!                                        Potential pressure term from interaction between particles i and j
!     -d  (distance) (in)
!                                        Distance between particles i and j
!
! Depencency:
!     -pbc(): Tool in boundary module
!==============================================================================!
    subroutine lj(r, boxlength, rc, ii, jj, F, pot, piter, d)
        double precision, allocatable, intent(in) :: r(:, :)
        double precision, allocatable, intent(inout) :: F(:, :)
        double precision, intent(in) :: boxlength, rc
        integer, intent(in) :: ii, jj
        double precision, intent(out) :: pot, piter
        double precision :: dx, dy, dz, d, dU

        pot = 0.d0; piter = 0.d0
        dx = r(ii, 1) - r(jj, 1); dy = r(ii, 2) - r(jj, 2); dz = r(ii, 3) - r(jj, 3)
        ! Apply the boundary conditions to the particles distance
        call pbc(dx, boxlength, 0.d0)
        call pbc(dy, boxlength, 0.d0)
        call pbc(dz, boxlength, 0.d0)
        d = (dx**2.+dy**2.+dz**2.)**0.5 ! Distance between particles i and j

        if (d .lt. rc) then
            dU = (48.d0/(d**14.0) - 24.d0/(d**8.0))

            F(ii, 1) = F(ii, 1) + dU*dx
            F(ii, 2) = F(ii, 2) + dU*dy
            F(ii, 3) = F(ii, 3) + dU*dz

            F(jj, 1) = F(jj, 1) - dU*dx
            F(jj, 2) = F(jj, 2) - dU*dy
            F(jj, 3) = F(jj, 3) - dU*dz

            pot = pot + 4.d0*(1.d0/(d**12.) - 1.d0/(d**6.))
            piter = piter + dU*dx; piter = piter + dU*dy; piter = piter + dU*dz
        end if
    end subroutine

    subroutine normal_rand(sigma, xout1, xout2)
        double precision,  intent(in) :: sigma
        double precision xout1, xout2
        double precision :: x(2)
        double precision, parameter :: PI = 4.d0*datan(1.d0)

        call random_number(x)

        xout1 = sigma*dsqrt(-2d0*(dlog(1d0 - x(1))))*dcos(2d0*PI*x(2))
        xout2 = sigma*dsqrt(-2d0*(dlog(1d0 - x(1))))*dsin(2d0*PI*x(2))

    end subroutine normal_rand

end module forces
