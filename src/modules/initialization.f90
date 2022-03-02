!==============================================================================!
!                             MODULE INITIALIZORS
! This module contains all the subroutines regarding the initialization of
! the spatial coordinates of the system and the velocities of the particles.
!
! The module needs the module $(params.f90) to work properly
! It includes a simple cubic initialization, fcc initialization and a bimodal
! velocity distribution.
!==============================================================================!

module initialization
!   use params
   implicit none
   integer :: ii, jj, kk, M
   double precision :: a
   integer :: nn, nx, ny, nz

   contains

!===========================================================================!
!                       DIAMOND (SC)
!===========================================================================!
subroutine diamond(boxlength,N, r)
double precision, allocatable, intent(inout)::r(:,:)
double precision, allocatable, dimension(:,:) :: r0
integer, intent(in) :: N
double precision, intent(in) :: boxlength

! Input variables: N - number of particles of the system
!                  boxlength - length of the box of simulation
! Output:          file with the configuration of a simple cubic network

! N = total number of nodes
! a = lattice spacing
! L = length of simulation box

M = int((float(N)/8.0d0)**(1.0/3.0)) ! M = units cells in each dimension
allocate(r0(8,3))
	
	
	a = boxlength/dfloat(M)


	r0(1,:) = [0.d0, 0.0d0, 0.0d0]
	r0(2,:)= [0.25d0,0.25d0,0.25d0]
	r0(3,:) = [0.0d0, 0.5d0, 0.5d0]
	r0(5,:) = [0.5d0,0.5d0,0.0d0]
	r0(6,:)=r0(3,:) + 0.25d0
	r0(7,:)=[0.5d0,0.0d0,0.5d0]
	r0(4,:)=r0(7,:) +0.25d0
	r0(8,:)=r0(5,:)+0.25d0


	open(11,file="../output/init_conf_diamond.xyz")
    write(11,*) N
	write(11,*) " "
   	nn = 1


   	do nz = 0, M - 1,1
   		do nx = 0, M - 1,1
   			do ny = 0, M - 1,1
				do ii = 1, 8
					r(nn,1) = ( nx + r0(ii,1) ) * a
					r(nn,2) = ( ny + r0(ii,2) ) * a
					r(nn,3) = ( nz + r0(ii,3) ) * a
					write(11,*) "A", r(nn,1), r(nn,2), r(nn,3)
         			nn = nn + 1
       			end do
   			end do
		end do
   	end do
	!Computing number of atoms.
   	nn = nn - 1

	close(11)
	deallocate(r0)


end subroutine diamond

endmodule initialization
