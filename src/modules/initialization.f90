!==============================================================================!
!                             MODULE INITIALIZATORS
! This module contains all the subroutines regarding the initialization of
! the spatial coordinates of the system and the velocities of the particles.
! 
!It includes a simple cubic initialization, fcc initialization and a bimodal
! velocity distribution.
!==============================================================================!

module initialization

   implicit none
   integer :: ii, jj, kk, M
   double precision :: a
   integer :: nn, nx, ny, nz

   contains

   !===========================================================================!
   !                       SIMPLE CUBIC CONFIGURATION (SC)
   !===========================================================================!
   subroutine sc(N,boxlength,r)
   double precision, allocatable, intent(inout)::r(:,:)
   integer, intent(inout) :: N
   integer::natoms
   double precision, intent(in) :: boxlength
   ! Input variables: N - number of particles of the system
   !                  L - length of the box of simulation
   ! Output:          file with the configuration of a simple cubic network

   ! N = total number of nodes
   ! a = lattice spacing
   ! D = dimensionality (D=3)
   ! L = length of simulation box
   !M = int(N**(1.0/3.0)) ! M = units cells in each dimension
   a = boxlength/dfloat(N)
   natoms=N*N*N
   open(11,file="output/init_conf_sc.xyz")
   write(11,*) natoms
   write(11,*) " "

   nn=0
   do ii = 0,N-1
     do jj = 0,N-1
       do kk = 0,N-1
       	 nn=nn+1
       	 r(nn,:)= [a*ii,a*jj,a*kk]
         write(11,*) "A", a*ii, a*jj, a*kk
       end do
     end do
   end do

   close(11)

   end subroutine sc

   !===========================================================================!
   !                   FACE CUBIC CENTERED CONFIGURATION (FCC)
   !===========================================================================!
   subroutine fcc(N,boxlength,r)
   double precision, allocatable, intent(inout)::r(:,:)
   double precision, allocatable, dimension(:,:) :: r0
   integer, intent(in) :: N
   integer::natoms
   double precision, intent(in) :: boxlength
   ! Input variables: N - number of particles of the system
   !                  L - length of the box of simulation
   ! Output:          file with the configuration of a simple cubic network
   !				  r- matrix with all the positions.
   
   !N = total number of nodes
   ! a = lattice spacing
   ! D = dimensionality (D=3)
   ! L = length of simulation box
   !M = int((float(N)/4.0)**(1.0/3.0)) ! M = units cells in each dimension
   allocate(r0(4,3))
   natoms=N*N*N*4
   a = boxlength/dfloat(N)

   
   r0(1,:) = [0.d0, 0.d0, 0.d0]
   r0(2,:) = [a/2.d0, a/2.d0, 0.d0]
   r0(3,:) = [0.d0, a/2.d0, a/2.d0]
   r0(4,:) = [a/2.d0, 0.d0, a/2.d0]

   ii = 0
   do nx = 0,N-1,1
     do ny = 0,N-1,1
       do nz = 0,N-1,1
         do jj = 1,4,1
           !print*, 4*ii+jj
           write(*,*)r0(jj,:)
           r(4*ii + jj,:) = a*[nx, ny, nz] + r0(jj,:)
         end do
         ii = ii+1
       end do
     end do
   end do

   open(12,file="output/init_conf_fcc.xyz")
   write(12,*) natoms
   write(12,*)
   do nn = 1,natoms,1
      write(12,*) "A", r(nn,1), r(nn,2), r(nn,3)
   end do
   close(12)
	deallocate(r0)
   end subroutine fcc



!===========================================================================!
!                             DIAMOND
!===========================================================================!
subroutine diamond(N,boxlength,r)
double precision, allocatable, intent(inout)::r(:,:)
double precision, allocatable, dimension(:,:) :: r0
integer, intent(in) :: N
double precision, intent(in) :: boxlength
integer::natoms

! Input variables: N - number of particles of the system
!                  boxlength - length of the box of simulation
! Output:          file with the configuration of a simple cubic network
!				   r- matrix with all the positions.

! N = total number of nodes
! a = lattice spacing
! L = length of simulation box

!M = int((float(N)/8.0d0)**(1.0/3.0)) ! M = units cells in each dimension
allocate(r0(8,3))
	
	natoms=N*N*N*8
	a = boxlength/dfloat(N)


	r0(1,:) = [0.d0, 0.0d0, 0.0d0]
	r0(2,:)= [0.25d0,0.25d0,0.25d0]
	r0(3,:) = [0.0d0, 0.5d0, 0.5d0]
	r0(5,:) = [0.5d0,0.5d0,0.0d0]
	r0(6,:)=r0(3,:) + 0.25d0
	r0(7,:)=[0.5d0,0.0d0,0.5d0]
	r0(4,:)=r0(7,:) +0.25d0
	r0(8,:)=r0(5,:)+0.25d0


	open(11,file="output/init_conf_diamond.xyz")
    write(11,*) natoms
	write(11,*) " "
   	nn = 1

   	do nz = 0, N - 1,1
   		do nx = 0, N - 1,1
   			do ny = 0, N - 1,1
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

   !===========================================================================!
   !                      BIMODAL VELOCITY DISTRIBUTION
   !===========================================================================!
   subroutine bimodal(Temp,vel)
   double precision, intent(in) :: Temp
   double precision, allocatable, dimension(:,:), intent(inout) :: vel
   double precision, allocatable, dimension(:) :: vel_decider
   integer :: ii, jj, dims, ind_pos, ind_neg, natoms

   natoms = size(vel,1); dims = size(vel,2)
   ind_pos = 0; ind_neg = 0

   allocate(vel_decider(natoms))
   call random_number(vel_decider)
   do ii = 1,natoms,1
     if ((vel_decider(ii).lt.0.5).and.(ind_pos.lt.int(natoms/2))) then
       ind_pos = ind_pos + 1
       do jj = 1,dims,1
         vel(ii,jj) = dsqrt(Temp)!/dsqrt(3.d0)
       end do
     elseif (ind_neg.lt.int(natoms/2)) then
       ind_neg = ind_neg + 1
       do jj = 1,dims,1
         vel(ii,jj) = -dsqrt(Temp)!/dsqrt(3.d0)
       end do
     else
       ind_pos = ind_pos + 1
       do jj = 1,dims,1
         vel(ii,jj) = dsqrt(Temp)!/dsqrt(3.d0)
       end do
     end if
   end do

   end subroutine bimodal

endmodule initialization
