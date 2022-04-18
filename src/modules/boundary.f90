!==============================================================================!
!                          MODULE BOUNDARY
!==============================================================================!
! Module that contains the procces to create periodic boundary conditions in
! a 3D squared system.
!==============================================================================!
! The module contains:														   !
!       -> pbc():															   !
!          Process that controls the periodic boundary conditions			   !
!                   Input:													   !
!                       - x (distance)(inout) : double precision scalar		   !
!                       - boxlength (in) : double precision scalar			   !
!                       - origin (in): double precision scalar				   !
!                   Output:													   !
!                       - x (distance)(inout) : double precision scalar        !
!==============================================================================!

module boundary
    implicit none

contains

!===========================================================================!
!                    PERIODIC BOUNDARY CONDITIONS							!
!===========================================================================!
! Periodic boundary conditions function for forces.
!===========================================================================!
! Input:																	!
!                - x (distance)(inout) : double precision scalar			!
!                - boxlength (in) : double precision scalar					!
!                - origin (in) : double precision scalar					!
!                                        Origin of the cell of simulation.	!
! Output:																	!
!                - origin (in) : double precision scalar					!
!                                        Origin of the cell of simulation.	!
!===========================================================================!

    subroutine pbc(x, boxlength, origin)
        double precision, intent(in) :: boxlength, origin
        double precision, intent(inout) :: x

        if (x .gt. (origin + boxlength/2.d0)) then
            x = x - boxlength

        elseif (x .lt. (origin - boxlength/2.d0)) then
            x = x + boxlength

        end if

    end subroutine pbc

end module boundary
