program main

use initialization
use boundary
use integrators

implicit none
include "../input/parameter.h"
integer::natoms
double precision::L
integer::tt,gg,si,sj
double precision::ti
double precision,allocatable::r(:,:),v(:,:),F(:,:)
double precision::ngr, pressp
double precision::epot,deltag
double precision, allocatable, dimension(:) ::  gr
integer::nhis

      
	!Initialization of the structure
	if (structure .eq. 1) then
		natoms=Nc*Nc*Nc
		L= (float(natoms)/density)**(1.0/3.0)
		write(*,*) L
		allocate(r(natoms,3),v(natoms,3),F(natoms,3))
		call sc(Nc,L,r)
	
	elseif (structure .eq. 2) then
		natoms=Nc*Nc*Nc*4
		L= (float(natoms)/density)**(1.0/3.0)
		allocate(r(natoms,3),v(natoms,3),F(natoms,3))
		call fcc(Nc,L,r)
		
	elseif (structure .eq. 3) then
		natoms=Nc*Nc*Nc*8
		L= (float(natoms)/density)**(1.0/3.0)
		allocate(r(natoms,3),v(natoms,3),F(natoms,3))
		call diamond( Nc,L,r)

	else 
		write(*,*)"Error, no structure found"
		stop
		
	endif
	
	
	
	
	nhis = 250; deltag = L/(2.d0*dble(nhis))
   allocate(gr(nhis)); gr(:) = 0.d0

	
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
     	
     	do si = 1,natoms 
            do sj = 1,3
                call pbc(r(si,sj),L,L/2.d0)
            end do 
        end do
   end do
	
	deallocate(r)
	deallocate(v)
	!deallocate(F)
	deallocate(gr)
	!deallocate(r,v,F, gr)
	
			
endprogram main
