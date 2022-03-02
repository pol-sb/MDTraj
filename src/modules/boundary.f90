module boundary
!   use params
   implicit none
   integer :: ii, jj, kk, M

   contains

		 !===========================================================================!
	   !               				PERIODIC BOUNDARY CONDITIONS
	   !===========================================================================!
	   ! Periodic boundary conditions function
	   subroutine pbc(x,boxlength,origin)
			 implicit none
	     double precision, intent(in) :: boxlength
	     double precision, intent(inout) :: x
	     if (x.gt.(boxlength/2.d0 + origin)) then
	       x = x - boxlength
	     elseif (x.lt.-(boxlength/2.d0 - origin)) then
	       x = x + boxlength
	     end if
	   end subroutine pbc

end module boundary
