!==============================================================================!
!                             MODULE INITIALIZORS
! This module contains all the subroutines regarding the initialization of
! the spatial coordinates of the system and the velocities of the particles.
!
! The module needs the module $(params.f90) to work properly
! It includes a simple cubic initialization, fcc initialization and a bimodal
! velocity distribution.
!==============================================================================!
! Contains
!
!
!
!
!
!
!
!
!
!
!
!
!
!===========================================================================!

module initialization
!   use params
   implicit none
   ! ## las variables que definamos aqui es por quermos guardarlas [ el comportamiento default es save]
   integer, save :: M
   double precision, save :: a
   integer, save :: nn

   contains

!===========================================================================!
!                       SIMPLE CUBIC CONFIGURATION (SC)
!===========================================================================!
! Input variables: N - number of particles of the system
!                  boxlength - length of the box of simulation
!                  r - empty array
! Output:          
!                  file with the configuration of a simple cubic network
!                  r - array with the configuration
!
! N = total number of nodes
! a = lattice spacing
! L = length of simulation box
!===========================================================================!
   subroutine initial_configuration_SC(boxlength,N, r)
   implicit none
      double precision, intent(out)::r(:,:)
      integer, intent(in) :: N
      double precision, intent(in) :: boxlength
      logical :: ext
      integer :: nx, ny, nz
      integer :: out_ref
      integer :: ii

      M = int(N**(1.0/3.0)) ! M = units cells in each dimension
      a = boxlength/dfloat(M)

      ! Creating the .xyz file with the FCC structure

      inquire(file="../output/",exist=ext)
      if (.NOT.ext) then
          call execute_command_line("mkdir ../output/")
      end if

      inquire(file="../output/init_conf_sc.xyz",exist=ext)
      if (.NOT.ext) then
          open(newunit=out_ref,file="../output/init_conf_sc.xyz", status="new")
      else 
          open(newunit=out_ref,file="../output/init_conf_sc.xyz", status="replace")
      end if

      nn = 1


outer:do nx = 0,M-1
         do ny = 0,M-1
            do nz = 0,M-1
                  r(nn,:)=(/a*nx, a*ny, a*nz/)

               if (nn.eq.(N)) then  ! Controller of the number of atoms
                  exit outer
               end if
               nn = nn + 1
            end do
         end do
      end do outer
      !Computing number of atoms.
      nn = nn-1

      write(out_ref,*) nn
      write(out_ref,*) " "      
      do ii =1, nn 
         write(out_ref,*) "A", r(ii,1), r(ii,2), r(ii,3)
      end do
      close(out_ref)
    end subroutine initial_configuration_SC

!===========================================================================!
!                   FACE CUBIC CENTERED CONFIGURATION (FCC)
!===========================================================================!
! Input variables: 
!                  N - number of particles of the system
!                  boxlength - length of the box of simulation
!                  r - empty array
! Output:          file with the configuration of a face centered cubic network
!
! N = total number of nodes
! a = lattice spacing
! L = length of simulation box
!===========================================================================!
   subroutine initial_configuration_fcc(boxlength,N,r)
   implicit none
      integer, intent(in) :: N
      double precision, intent(in) :: boxlength
      double precision,allocatable, intent(out) :: r(:,:)
      double precision, allocatable :: r0(:,:)

      logical :: ext
      integer :: nx, ny, nz
      integer :: out_ref, ii, jj

      M = int((float(N)/4.0)**(1.0/3.0)) ! M = units cells in each dimension

      a = boxlength/dfloat(M)

      inquire(file="../output/",exist=ext)
      if (.NOT.ext) then
          call execute_command_line("mkdir ../output/")
      end if

      inquire(file="../output/init_conf_fcc.xyz",exist=ext)
      if (.NOT.ext) then
          open(newunit=out_ref,file="../output/init_conf_fcc.xyz", status="new")
      else 
          open(newunit=out_ref,file="../output/init_conf_fcc.xyz", status="replace")
      end if

      M = int((float(N)/4.0)**(1.0/3.0)) ! M = units cells in each dimension
      allocate(r(3,N),r0(3,4))
      a = boxlength/dfloat(M)
   
      r0(:,1) = [0.d0, 0.d0, 0.d0]
      r0(:,2) = [a/2.d0, a/2.d0, 0.d0]
      r0(:,3) = [0.d0, a/2.d0, a/2.d0]
      r0(:,4) = [a/2.d0, 0.d0, a/2.d0]
   
      ii = 0
      do nx = 0,M-1,1
        do ny = 0,M-1,1
          do nz = 0,M-1,1
            do jj = 1,4,1
              !print*, 4*ii+jj
              r(:,4*ii + jj) = a*[nx, ny, nz] + r0(:,jj)
            end do
            ii = ii+1
          end do
        end do
      end do
   
      write(out_ref,*) N
      write(out_ref,*)
      do nn = 1,N
         write(out_ref,*) "A", r(1,nn), r(2,nn), r(3,nn)
      end do
      close(out_ref)
   end subroutine initial_configuration_fcc


!===========================================================================!
!                              DIAMOND
!===========================================================================!
! Input variables: 
!                  N - number of particles of the system
!                  boxlength - length of the box of simulation
!                  r - empty array
! Output:          
!                  file with the configuration of a diamond network
!                  r - array ith the configuration
!
! N = total number of nodes
! a = lattice spacing
! L = length of simulation box
!===========================================================================!
   subroutine initial_configuration_diamond(boxlength,N,r)
   implicit none
      integer, intent(in) :: N
      double precision, intent(in) :: boxlength
      double precision,intent(out) :: r(:,:)
      double precision, allocatable :: r0(:,:)
      logical :: ext
      integer :: nx, ny, nz
      integer :: out_ref
      integer :: ii




      M = int((float(N)/8.0d0)**(1.0/3.0)) ! M = units cells in each dimension
      print*,M, float(N)/8.0d0
      a = boxlength/dfloat(M)
      allocate(r0(8,3))

      ! Vector of all atoms contained in a unit cell

      r0(1,:) = [0.d0, 0.0d0, 0.0d0]
      r0(2,:)= [0.25d0,0.25d0,0.25d0]
      r0(3,:) = [0.0d0, 0.5d0, 0.5d0]
      r0(5,:) = [0.5d0,0.5d0,0.0d0]
      r0(6,:)=r0(3,:) + 0.25d0
      r0(7,:)=[0.5d0,0.0d0,0.5d0]
      r0(4,:)=r0(7,:) +0.25d0
      r0(8,:)=r0(5,:)+0.25d0


      ! Creating the .xyz file with the diamond structure

      inquire(file="../output/",exist=ext)
      if (.NOT.ext) then
          call execute_command_line("mkdir ../output/")
      end if

      inquire(file="../output/init_conf_diamond.xyz",exist=ext)
      if (.NOT.ext) then
          open(newunit=out_ref,file="../output/init_conf_diamond.xyz", status="new")
      else 
          open(newunit=out_ref,file="../output/init_conf_diamond.xyz", status="replace")
      end if

      nn = 1

outer:do nz = 0, M - 1,1
         do nx = 0, M - 1,1
            do ny = 0, M - 1,1
               do ii = 1, 8
                  
                  r(nn,:) = (/(( nx + r0(ii,1) ) * a), &
                            ((ny + r0(ii,2) ) * a), &
                            ((nz + r0(ii,3) ) * a)/)
                  nn = nn + 1
               end do
            end do
         end do
      end do outer
      !Computing number of atoms.
      nn = nn - 1

      write(out_ref,*) nn
      write(out_ref,*) " "
      do ii =1, nn 
         write(out_ref,*) "A", r(ii,1), r(ii,2), r(ii,3)
      end do
      close(out_ref)
      deallocate(r0)


   end subroutine initial_configuration_diamond

!===========================================================================!
!                       READ FROM FILE
!===========================================================================!
! Input variables: 
!                  N (in)- number of particles of the system
!                  coord_path (in) - path to the xyz fle containing the structure
!                  vel_path(inout) (OPTIONAL) : path to velocities file
!                  initial_velocities(inout)(OPTIONAL) : empty array
!                  initial_position(inout) : empty array
! Output:
!                  initial_velocities(inout)(OPTIONAL) : empty array
!                  initial_position(inout) : empty array
!
!===========================================================================!
subroutine initial_reading(N, coord_path, initial_position, initial_velocities, vel_path)
   implicit none
       integer,intent(in) :: N
       real(8), intent(out) :: initial_position(:,:)
       character(len=*), intent(in) :: coord_path
       real(8), optional, intent(out) :: initial_velocities(:,:)
       character(len=*),optional,intent(in) :: vel_path
   
       ! Internal Parameters declaration
       character(len=4) ::atom_name
       real(8) :: x,y,z
       logical :: ext
       integer :: error
       integer :: file_id
       integer :: i


      !////////////// Read the old coordinates
       open(newunit=file_id,file=coord_path)
       read(file_id,*)
       read(file_id,*)
       do i=1,N
           read(file_id,*,iostat=error)atom_name,x,y,z
           if (error>0) then
               print*, "Error in ",coord_path,"reading"
           else if (error<0) then
               exit
           else
               initial_position(i,:)=(/(x),(y),(z)/)
           end if
       end do
       close(file_id)
   
   
       !////////////// Read the old velocities
       if (present(vel_path)) then 
          open(newunit=file_id,file=vel_path)
          read(file_id,*)
          read(file_id,*)
          do i=1,N
              read(file_id,*,iostat=error) x,y,z
              if (error>0) then
              print*, "Error in ",vel_path,"reading"            
              else if (error<0) then
                  exit
              else
                  initial_velocities(i,:)= (/x,y,z/)
              end if
          end do
         end if
   end subroutine initial_reading
   
end module initialization
   
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

end module initialization
