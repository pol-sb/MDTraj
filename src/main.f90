program main

use initialization

implicit none
include "../input/parameter.h"
     
double precision,allocatable::r(:,:),v(:,:)

	allocate(r(N,3),v(N,3))
      
	if (structure .eq. 0) then
		!read from file
		
		
	elseif (structure .eq. 1) then
		! Call SC
	
	elseif (structure .eq. 2) then
		! Call FCC
		
	elseif (structure .eq. 3) then
		call diamond(L, N, r)

	else 
		write(*,*)"Error, no structure found"
		stop
		
	endif
	
	
	
	
	
	
	
			
endprogram main
