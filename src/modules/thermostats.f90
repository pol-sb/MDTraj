module thermostat
	 implicit none
   contains

   !===========================================================================!
   !                              ANDERSEN THERMOSTAT
   !===========================================================================!
   subroutine andersen_thermo(Temp,vel,natoms)
   include "constants.h"
   double precision, dimension(:,:) :: vel
   integer,intent(in)::natoms
   double precision :: nu, Temp, sigma
   double precision, dimension(size(vel,2)) :: x_rand
   double precision, dimension(4) :: vel_normalrand

   nu = 1e-3
   sigma = dsqrt(Temp)
   call random_number(x_rand)

   do ii = 1,natoms
     if (x_rand(ii).lt.nu) then ! choosing if velocity of particle i gets changed
       call normal_rand(sigma,vel_normalrand(1),vel_normalrand(2))
       call normal_rand(sigma,vel_normalrand(3),vel_normalrand(4))
       ! The subroutine normal_rand returns a random number from a normal
       ! distribution with standard deviation \sigma = sqrt(T)
       do jj = 1,3
         vel(jj,ii) = vel_normalrand(jj)
       end do
     end if
   end do

   end subroutine

   !===========================================================================!
   !                        NORMAL RANDOM NUMBER GENERATOR
   !===========================================================================!
   subroutine normal_rand(sigma, xout1, xout2)
   include "constants.h"
   double precision :: sigma
   double precision xout1, xout2
   double precision :: x(2)

   call random_number(x)

   xout1 = sigma*dsqrt(-2d0*(dlog(1d0-x(1))))*dcos(2d0*PI*x(2))
   xout2 = sigma*dsqrt(-2d0*(dlog(1d0-x(1))))*dsin(2d0*PI*x(2))

   end subroutine normal_rand


   !===========================================================================!
   !                        NORMAL RANDOM NUMBER GENERATOR
   !===========================================================================!
   double precision function kinetic(vel,natoms) result(ekin)
   include "constants.h"
   integer,intent(in)::natoms
   double precision, allocatable, dimension(:,:), intent(in) :: vel
     do ii = 1,natoms
       ekin = ekin + 0.5d0*(vel(1,ii)**2 + vel(2,ii)**2 + vel(3,ii)**2)
     end do
   end function kinetic
endmodule thermostat
