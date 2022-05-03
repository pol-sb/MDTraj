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
!                                             - r (positions)(inout) : double precision array
!                                        - F (force)(inout) : double precision array
!                                              -dt (time step)(in) : double precision scalar
!                                             -boxlength (longitude of one side)(in) double precision scalar
!                   Output:
!                                            - r (positions)(inout) : double precision array
!                      - vel (velocity)(inout) : double precision array
!
!
!       -> vel_verlet (natoms,r,vel,F,Upot,dt,rc,boxlength,pressp,gr,deltag):
!           The verlet algorithm method wihout considering the temperature.
!                        It returns the new positions and velocities.
!                   Input:
!                      - natoms (número de átomos)(in) : integer scalar
!                                               - r (positions)(inout) : double precision array
!                                         - vel (velocity)(inout) : double precision array
!                                         - F (force)(inout) : double precision array
!                                               -dt (time step)(in) : double precision scalar
!                                              -rc (cut-off distance) (in) : double precision
!                                              -boxlength (longitude of one side)(in) : double precision scalar
!                                              -gr ()(in) : double precision array
!                                         -deltag ()(in) : double precision scalar
!                   Output:
!                     - r (positions)(inout) : double precision array
!                                          - vel (velocity)(inout) : double precision array
!                                               - Upot (potential energy) (out) : double precision scalar
!                                               - Pressp (pressure) (out) : double precision scalar
!                                                  Depencency:
!                                                - force() : In module forces (src/modules/forces.f90)
!
!
!       -> vel_verlet_with_thermo (natoms,r,vel,F,Upot,dt,rc,boxlength,Temp,pressp,gr,deltag):
!           The verlet algorithm method considering the temperature.
!                        It returns the new positions and velocities.
!                   Input:
!                      - natoms (número de átomos)(in) : integer scalar
!                                               - r (positions)(inout) : double precision array
!                                         - vel (velocity)(inout) : double precision array
!                                         - F (force)(inout) : double precision array
!                                               -dt (time step)(in) : double precision scalar
!                                              -rc (cut-off distance) (in) : double precision
!                                              -boxlength (longitude of one side)(in) : double precision scalar
!                                              -gr ()(in) : double precision array
!                                         -deltag ()(in) : double precision scalar
!                   Output:
!                      - r (positions)(inout) : double precision array
!                                         - vel (velocity)(inout) : double precision array
!                                              - Upot (potential energy)(out) : double precision scalar
!                                               - temp (temperature)(in) : double precision scalar
!                                              - Pressp (pressure) (out) : double precision scalar
!                                                  Depencency:
!                                                - force() : In module forces (src/modules/forces.f90)
!                                                           - andersen_thermo() : In module thermostat (src/modules/thermostats.f90)
!=====================================================================================!
module integrators
    use forces
    use thermostat
    implicit none

contains

!=====================================================================================!
!                                     EULER INTEGRATION
!=====================================================================================!
!  Input:
!      - natoms (número de átomos)(in) : integer scalar
!           - r (positions)(inout) : double precision array
!      - F (force)(inout) : double precision array
!            -dt (time step)(in) : double precision scalar
!           -boxlength (longitude of one side)(in) double precision scalar
!  Output:
!            - r (positions)(inout) : double precision array
!      - vel (velocity)(inout) : double precision array
!=====================================================================================!
    subroutine euler(natoms,r,vel,F,Upot,dt,rc,boxlength,pressp)
        integer, intent(in)::natoms
        double precision, allocatable, intent(inout) :: F(:, :)
        double precision, allocatable, intent(inout) :: r(:, :), vel(:, :)
        double precision, intent(in) :: dt, boxlength, rc
		double precision, intent(out) :: Upot, pressp

        do jj = 1, natoms
            do ii = 1, 3
                r(jj, ii) = r(jj, ii) + vel(jj, ii)*dt + 0.5d0*F(jj, ii)*dt*dt
                vel(jj, ii) = vel(jj, ii) + F(jj, ii)*dt
            end do
        end do
		call force_dpd(natoms, r, boxlength, rc, F, Upot, pressp, dt, vel)

    end subroutine euler

!=====================================================================================!
!                     VELOCITY VERLET INTEGRATION (Without temperature)
!=====================================================================================!
!  Input:
!      - natoms (número de átomos)(in) : integer scalar
!            - r (positions)(inout) : double precision array
!      - vel (velocity)(inout) : double precision array
!      - F (force)(inout) : double precision array
!            -dt (time step)(in) : double precision scalar
!           -rc (cut-off distance) (in) : double precision
!           -boxlength (longitude of one side)(in) : double precision scalar
!           -gr ()(in) : double precision array
!      -deltag ()(in) : double precision scalar
!  Output:
!           - r (positions)(inout) : double precision array
!      - vel (velocity)(inout) : double precision array
!           - Upot (potential energy) (out) : double precision scalar
!           - Pressp (pressure) (out) : double precision scalar
!  Depencency:
!      - force() : In module forces (src/modules/forces.f90)
!=====================================================================================!
    subroutine vel_verlet(natoms, r, vel, F, Upot, dt, rc, boxlength, pressp, gr, deltag)
        integer, intent(in)::natoms
        double precision, allocatable, intent(inout) :: F(:, :)
        double precision, allocatable, intent(inout) :: r(:, :), vel(:, :)
        double precision, allocatable, intent(inout) :: gr(:)
        double precision, intent(in) :: dt, rc, boxlength, deltag
        double precision, intent(out) :: Upot, pressp

        do ii = 1, natoms
            do jj = 1, 3
                r(jj, ii) = r(jj, ii) + vel(jj, ii)*dt + 0.5d0*F(jj, ii)*dt*dt
                vel(jj, ii) = vel(jj, ii) + F(jj, ii)*0.5d0*dt
            end do
        end do

        call force_dpd(natoms, r, boxlength, rc, F, Upot, pressp, dt, vel)

        do ii = 1, natoms
            do jj = 1, 3
                vel(jj, ii) = vel(jj, ii) + F(jj, ii)*0.5d0*dt
            end do
        end do

    end subroutine vel_verlet

!=====================================================================================!
!                     VELOCITY VERLET INTEGRATION (With temperature)
!=====================================================================================!
!  Input:
!      - natoms (número de átomos)(in) : integer scalar
!            - r (positions)(inout) : double precision array
!      - vel (velocity)(inout) : double precision array
!      - F (force)(inout) : double precision array
!            -dt (time step)(in) : double precision scalar
!           -rc (cut-off distance) (in) : double precision
!           -boxlength (longitude of one side)(in) : double precision scalar
!           -gr ()(in) : double precision array
!      -deltag ()(in) : double precision scalar
!  Output:
!           - r (positions)(inout) : double precision array
!      - vel (velocity)(inout) : double precision array
!           - Upot (potential energy)(out) : double precision scalar
!            - temp (temperature)(in) : double precision scalar
!           - Pressp (pressure) (out) : double precision scalar
!  Depencency:
!      - force() : In module forces (src/modules/forces.f90)
!      - andersen_thermo() : In module thermostat (src/modules/thermostats.f90)
!=====================================================================================!
    subroutine vel_verlet_with_thermo(natoms, r, vel, F, Upot, dt, rc, boxlength, Temp, pressp, gr, deltag)
        integer, intent(in)::natoms
        double precision, allocatable, intent(inout) :: F(:, :)
        double precision, allocatable, intent(inout) :: r(:, :), vel(:, :)
        double precision, allocatable, intent(inout) :: gr(:)
        double precision, intent(in) :: dt, rc, boxlength, Temp, deltag
        double precision, intent(out) :: Upot, pressp

        Upot = 0.d0; pressp = 0.d0

        do jj = 1, natoms
            do ii = 1, 3
                r(jj, ii) = r(jj, ii) + vel(jj, ii)*dt + 0.5d0*F(jj, ii)*dt*dt
                vel(jj, ii) = vel(jj, ii) + F(jj, ii)*0.5d0*dt
            end do
        end do

        call force(natoms, r, boxlength, rc, F, Upot, pressp, gr, deltag)

        do jj = 1, natoms
            do ii = 1, 3
                vel(jj, ii) = vel(jj, ii) + F(jj, ii)*0.5d0*dt
            end do
        end do

        call andersen_thermo(Temp, vel, natoms)

    end subroutine vel_verlet_with_thermo

end module integrators

