program main

use	initialization
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
double precision::epot, ekin, temperature, deltag
double precision::rpos, vb, nid
double precision, allocatable, dimension(:) ::  gr
integer::nhis


	!Initialization of the structure
	if (structure .eq. 1) then
		natoms=Nc*Nc*Nc
		L= (float(natoms)/density)**(1.0/3.0)
		write(*,*) L
		allocate(r(natoms,3),v(natoms,3),F(natoms,3))
		call initial_configuration_SC(Nc,L,r)

	elseif (structure .eq. 2) then
		natoms=Nc*Nc*Nc*4
		L= (float(natoms)/density)**(1.0/3.0)
		allocate(r(natoms,3),v(natoms,3),F(natoms,3))
		call initial_configuration_fcc(Nc,L,r)

	elseif (structure .eq. 3) then
		natoms=Nc*Nc*Nc*8
		L= (float(natoms)/density)**(1.0/3.0)
		allocate(r(natoms,3),v(natoms,3),F(natoms,3))
		call initial_configuration_diamond( Nc,L,r)

	else
		write(*,*)"Error, no structure found"
		stop

	endif

	nhis = 250; deltag = L/(2.d0*dble(nhis)); rc = L/2.d0
   allocate(gr(nhis)); gr(:) = 0.d0

	!initialization of velocity
	if (vel_opt .eq. 1) then
		call bimodal(Temp,v)
	else
		v(:,:)=0.0d0
	endif

	call force(natoms,r,L,rc,F,epot,pressp,gr,deltag)

	open(11,file='output/temp.dat',status='unknown')
	open(12,file='output/energy.dat',status='unknown')
	open(13,file='output/pressure.dat',status='unknown')

	do tt = 1,ntimes,1
		ti = ti+dt ! Actualizing the instant time
    if (thermo.eq.0) then
      call vel_verlet(natoms,r,v,F,epot,dt,rc,L,pressp,gr,deltag)
    elseif (thermo.eq.1) then
     	call vel_verlet_with_thermo(natoms,r,v,F,epot,dt,rc,L,Temp,pressp,gr,deltag)
    else
      write(*,*)"Error, no thermostat status found"
  		stop
    endif

		ekin = kinetic(v,natoms)
		temperature = 2.d0*ekin/(3.d0*dble(natoms)-3.d0)

    ngr = ngr+1

   	do si = 1,natoms
      do sj = 1,3
        call pbc(r(si,sj),L,L/2.d0)
      end do
    end do

		if (mod(tt,everyt).eq.0) then
			write(11,*) ti, temperature
			write(12,*) ti, epot, ekin, epot+ekin
			write(13,*) ti, pressp, density*temperature, pressp+density*temperature
		end if
   end do

	open(14,file='output/rdf.dat',status='unknown')
	do ii= 1,nhis
		rpos = deltag*(ii + 0.5) ! Distance r
		vb = ((ii+1)**3 - ii**3)*(deltag**3)
		! Volume between bin i+1 and i
		nid = (4.d0*PI/3.d0)*vb*density
		! Number of ideal gas part . in vb
		gr(ii) = gr(ii)/(dble(ngr)*dble(natoms)*nid) ! Normalize g(r)
		write(14,*) rpos*sigma, gr(ii)
	enddo
	close(14)

	deallocate(r,v,F, gr)

endprogram main
