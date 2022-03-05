module boundary

   implicit none
 
   contains
   
   !===========================================================================!
   !               				PERIODIC BOUNDARY CONDITIONS
   !===========================================================================!
   ! Periodic boundary conditions function
	   subroutine pbc(x,boxlength,origin)
			 implicit none
	     double precision, intent(in) :: boxlength,origin
	     double precision, intent(inout) :: x
	     if (x.gt.(origin + boxlength/2.d0)) then
	       x = x - boxlength
	     elseif (x.lt.(origin - boxlength/2.d0)) then
	       x = x + boxlength
	     end if
	   end subroutine pbc

end module boundary
