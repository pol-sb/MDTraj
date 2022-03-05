program main

use initialization
use boundary
use integrators

implicit none
include "../input/parameter.h"
integer::tt,gg
double precision::ti
double precision,allocatable::r(:,:),v(:,:),F(:,:)
double precision::ngr, pressp
double precision::epot,deltag
double precision, allocatable, dimension(:) :: r2t, gr
integer::nhis

	nhis = 250; deltag = L/(2.d0*dble(nhis))
   allocate(gr(nhis)); gr = 0.d0

	allocate(r(natoms,3),v(natoms,3),F(natoms,3))
      
	!Initialization of the structure
	if (structure .eq. 1) then
		call sc(natoms,L,r)
	
	elseif (structure .eq. 2) then
		call fcc(natoms,L,r)
		
	elseif (structure .eq. 3) then
		call diamond( natoms,L,r)

	else 
		write(*,*)"Error, no structure found"
		stop
		
	endif
	
	!initialization of velocity
	if (vel_opt .eq. 1) then
		call bimodal(Temp,v)
	else
		v(:,:)=0.0d0
	endif
	
	call force(natoms,r,L,rc,F,epot,pressp,gr,deltag)
	do tt = 1,ntimes,1
		ti = ti+dt ! Actualizing the instant time
     	call vel_verlet_with_thermo(natoms,r,v,F,epot,dt,rc,L,Temp,pressp,gr,deltag)
     	ngr = ngr+1
     	
     	if (gg.lt.1e4) then
       		gr = 0.d0
     	end if
     	
   end do
	
	
	
	
			
endprogram main
