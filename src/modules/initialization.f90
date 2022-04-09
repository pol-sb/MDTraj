!=====================================================================================!
!                             MODULE INITIALIZORS
! This module contains all the subroutines regarding the initialization of
! the spatial coordinates of the execute_command_line and the velocities of the particles.
!
! It includes a simple cubic (sc), an face centered ccubic (FCC) and diamond
! cristalline initializarion of the structure.
!
! It is included as well a subroutine to read an structure from a file.
!
! It also contain a bimodal velocity distribution.
!=====================================================================================!
! Contains
!
! 	  -> initial_configuration_SC (N, boxlength, r):
! 		  This subroutine creates a SC cristalline structure.
!				Input:
!					- N (number particles of one side)(in): inegter scalar
!					- Boxlength (longitude of one side)(in): double precision scalar
!				Output:
!					- r (positions of the atoms)(out): double precision array
!
!
!     -> initial_configuration_fcc (N, boxlength, r):
! 		  This subroutine creates a FCC cristalline structure.
!				Input:
!					- N (number particles of one side)(in): inegter scalar
!					- Boxlength (longitude of one side)(in): double precision scalar
!				Output:
!					- r (positions of the atoms)(out): double precision array
!
!
!      -> initial_configuration_diamond (N, boxlength, r):
!		   This subroutine creates a SC cristalline structure.
!				Input:
!					- N (number particles of one side)(in): inegter scalar
!					- Boxlength (longitude of one side)(in): double precision scalar
!				Output:
!					- r (positions of the atoms)(out): double precision array
!
!
!      -> initial_reading (N, coord_path, initial_position, initial_velocities, vel_path):
!		   This subroutine creates a SC cristalline structure.
!				Input:
!				   - N (number of sides of the execute_command_line)(in): integer scalar
!    			   - coord_path (path to the xyz fle containing the structure)(in): character
!	 			   - vel_path (path to velocities file)(inout): OPTIONAL, empty array
! 	 			   - initial_velocities (velocidades iniciales)(inout): OPTIONAL, empty array
! 	       		   - initial_position (posiciones iniciales)(inout): OPTIONAL, empty array
!				Output:
!					- initial_velocities (velocidades iniciales)(inout): OPTIONAL, empty array
! 				    - initial_position (posiciones iniciales)(inout): OPTIONAL, empty array
!
!
!      -> bimodal (Temp,vel):
!		   This subroutine generates the velocities of the atoms following a bimodal
!		   distribution of a temperature.
!				Input:
!					- Temp (temperature)(in): double precision scalar
!				Output:
!					- vel (velocities of the atoms)(inout): double precision array
!=====================================================================================!

module initialization
	implicit none

   double precision, save :: a
   integer, save :: nn

contains

!=====================================================================================!
!                       SIMPLE CUBIC CONFIGURATION (SC)
!=====================================================================================!
! Input:
!	 - N (number particles of one side)(in): inegter scalar
!	 - Boxlength (longitude of one side)(in): double precision scalar
!
! Output:
!	 - r (positions of the atoms)(out): double precision array
!=====================================================================================!
	subroutine initial_configuration_SC(N,boxlength, r, rank)
	double precision, intent(out)::r(:,:)
      integer, intent(in) :: N,rank
      double precision, intent(in) :: boxlength
      logical :: ext
      integer :: nx, ny, nz
      integer :: out_ref
      integer :: ii
      integer :: natoms
      double precision :: r_pert(size(r,2))

        
      a = boxlength/dfloat(N)
      natoms = N*N*N

      ! Creating the .xyz file with the FCC structure

      nn = 1

	print*, size(r_pert)

outer:do nx = 0,N-1
         do ny = 0,N-1
            do nz = 0,N-1
	       
	     call random_number(r_pert)
	     r_pert = (r_pert - 0.5d0)*a/4.d0
            r(nn,:)=(/a*nx + r_pert(1), a*ny + r_pert(2), a*nz + r_pert(3)/)
            nn = nn+1
            
            enddo
         enddo
      enddo outer

  
      if (rank.eq.0) then
	      inquire(file="./output/",exist=ext)
	      if (.NOT.ext) then
		   call execute_command_line("mkdir ./output/")
	      endif

	      inquire(file="./output/structure",exist=ext)
	      if (.NOT.ext) then
		   call execute_command_line("mkdir ./output/structure")
	      endif

	      inquire(file="./output/structure/init_conf_sc.xyz",exist=ext)
	      if (.NOT.ext) then
		   open(newunit=out_ref,file="./output/structure/init_conf_sc.xyz", status="new")
	      else
		   open(newunit=out_ref,file="./output/structure/init_conf_sc.xyz", status="replace")
	      endif
	      write(out_ref,*) natoms
	      write(out_ref,*) " "
	      do ii = 1,natoms
		  write(out_ref,*) "A", r(ii,1), r(ii,2), r(ii,3)
	      enddo

	      close(out_ref)
      endif

    endsubroutine initial_configuration_SC

!=====================================================================================!
!                   FACE CUBIC CENTERED CONFIGURATION (FCC)
!=====================================================================================!
! Input:
!	 - N (number particles of one side)(in): inegter scalar
!	 - Boxlength (longitude of one side)(in): double precision scalar
!
! Output:
!	 - r (positions of the atoms)(out): double precision array
!=====================================================================================!
   subroutine initial_configuration_fcc(N,boxlength,r,rank)
   	include "../declaration_variables/parallel_variables.h"
      integer, intent(in) :: N,rank
      double precision, intent(in) :: boxlength
      double precision,allocatable, intent(out) :: r(:,:)
      double precision, allocatable :: r0(:,:)

      logical :: ext
      integer :: nx, ny, nz, natoms
      integer :: out_ref, ii, jj
			double precision :: r_pert(size(r,2))

      a = boxlength/dfloat(N)
      natoms=N*N*N*4

 

      allocate(r0(4,3))

      r0(1,:) = [0.d0, 0.d0, 0.d0]
      r0(2,:) = [a/2.d0, a/2.d0, 0.d0]
      r0(3,:) = [0.d0, a/2.d0, a/2.d0]
      r0(4,:) = [a/2.d0, 0.d0, a/2.d0]

      nn = 0
      ii = 0
      do nx = 0,N-1,1
        do ny = 0,N-1,1
          do nz = 0,N-1,1
            do jj = 1,4,1
              !print*, 4*ii+jj
							call random_number(r_pert)
							r_pert = (r_pert - 0.5d0)*a/4.d0
              r(4*ii + jj,:) = a*[nx, ny, nz] + r0(jj,:) + r_pert(:)
            enddo
            ii = ii+1
          enddo
        enddo
      enddo


      if (rank.eq.0) then
	      inquire(file="./output/",exist=ext)
	      if (.NOT.ext) then
		   call execute_command_line("mkdir ./output/")
	      endif

	      inquire(file="./output/structure",exist=ext)
	      if (.NOT.ext) then
		   call execute_command_line("mkdir ./output/structure/structure")
	      endif

	      inquire(file="./output/init_conf_fcc.xyz",exist=ext)
	      if (.NOT.ext) then
		   open(newunit=out_ref,file="./output/structure/init_conf_fcc.xyz", status="new")
	      else
		   open(newunit=out_ref,file="./output/structure/init_conf_fcc.xyz", status="replace")
	      endif
	      write(out_ref,*) natoms
	      write(out_ref,*)

	      do ii = 1,nn
		  write(out_ref,*) "A", r(nn,1), r(nn,2), r(nn,3)
	      enddo

	      close(out_ref)
      endif
      
      deallocate(r0)

   endsubroutine initial_configuration_fcc

!=====================================================================================!
!                              DIAMOND
!=====================================================================================!
! Input:
!	 - N (number particles of one side)(in): inegter scalar
!	 - Boxlength (longitude of one side)(in): double precision scalar
!
! Output:
!	 - r (positions of the atoms)(out): double precision array
!=====================================================================================!
   subroutine initial_configuration_diamond(N,boxlength,r,rank)
      integer, intent(in) :: N,rank
      double precision, intent(in) :: boxlength
      double precision,intent(out) :: r(:,:)
      double precision, allocatable :: r0(:,:)
      logical :: ext
      integer :: nx, ny, nz
      integer :: out_ref
      integer :: ii, natoms
      double precision :: r_pert(size(r,2))

      a = boxlength/dfloat(N)
      natoms=8*N*N*N
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


      

      nn = 1

outer:do nz = 0, N - 1,1
         do nx = 0, N - 1,1
            do ny = 0, N - 1,1
               do ii = 1, 8
								 call random_number(r_pert)
								 r_pert = (r_pert - 0.5d0)*a/4.d0
                  r(nn,:) = (/(( nx + r0(ii,1) ) * a + r_pert(1)), &
                            ((ny + r0(ii,2) ) * a + r_pert(2)), &
                            ((nz + r0(ii,3) ) * a + r_pert(3))/)
                  nn = nn + 1
               enddo
            enddo
         enddo
      enddo outer

      !Computing number of atoms.
      nn = nn - 1
  
      if (rank.eq.0) then
	      ! Creating the .xyz file with the diamond structure

	      inquire(file="./output/",exist=ext)
	      if (.NOT.ext) then
		   call execute_command_line("mkdir ./output/")
	      endif
	      inquire(file="./output/structure",exist=ext)
	      if (.NOT.ext) then
		   call execute_command_line("mkdir ./output/structure")
	      endif
	      inquire(file="./output/structure/init_conf_diamond.xyz",exist=ext)
	      if (.NOT.ext) then
		   open(newunit=out_ref,file="./output/structure/init_conf_diamond.xyz", status="new")
	      else
		   open(newunit=out_ref,file="./output/structure/init_conf_diamond.xyz", status="replace")
	      endif

	      write(out_ref,*) natoms
	      write(out_ref,*) " "

	      do ii =1, nn !
		  write(out_ref,*) "A", r(ii,1), r(ii,2), r(ii,3)
	      enddo

	      close(out_ref)
      endif
      deallocate(r0)

   endsubroutine initial_configuration_diamond

!=====================================================================================!
!                       READ FROM FILE
!=====================================================================================!
!=====================================================================================!
! Input:
!	 - N (number of sides of the execute_command_line)(in): integer scalar
!    - coord_path (path to the xyz fle containing the structure)(in): character
!	 - vel_path (path to velocities file)(inout): OPTIONAL, empty array
! 	 - initial_velocities (velocidades iniciales)(inout): OPTIONAL, empty array
! 	 - initial_position (posiciones iniciales)(inout): OPTIONAL, empty array
!
! Output:
!	 - initial_velocities (velocidades iniciales)(inout): OPTIONAL, empty array
! 	 - initial_position (posiciones iniciales)(inout): OPTIONAL, empty array
!=====================================================================================!

	subroutine initial_reading(N, coord_path, initial_position, initial_velocities, vel_path)
	
	integer,intent(in) :: N
      double precision, intent(out) :: initial_position(:,:)
      character(len=*), intent(in) :: coord_path
      double precision, optional, intent(out) :: initial_velocities(:,:)
      character(len=*),optional,intent(in) :: vel_path

      ! Internal Parameters declaration
      character(len=4) ::atom_name
      double precision :: x,y,z
      logical :: ext
      integer :: error
		integer :: file_id
      integer :: i


      !////////////// Read the old coordinates
      open(unit=101,file=coord_path)
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
           endif
       enddo

       close(file_id)


       !////////////// Read the old velocities
      if (present(vel_path)) then
			open(unit=102,file=vel_path)
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
				endif
			enddo
		endif

   endsubroutine initial_reading

!=====================================================================================!
!                      BIMODAL VELOCITY DISTRIBUTION
!=====================================================================================!
! Input:
!	 - Temp (temperature)(in): double precision scalar
!
! Output:
!	 - vel (velocities of the atoms)(inout): double precision array
!=====================================================================================!
	subroutine bimodal(Temp,vel)
   	double precision, intent(in) :: Temp
   	double precision, allocatable, intent(inout) :: vel(:,:)
   	double precision, allocatable :: vel_decider(:)
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
       		enddo
     		elseif (ind_neg.lt.int(natoms/2)) then
       		ind_neg = ind_neg + 1
       		do jj = 1,dims,1
         		vel(ii,jj) = -dsqrt(Temp)!/dsqrt(3.d0)
       		enddo
     		else
       		ind_pos = ind_pos + 1
       		do jj = 1,dims,1
         		vel(ii,jj) = dsqrt(Temp)!/dsqrt(3.d0)
       		enddo
     		endif
   	enddo

   endsubroutine bimodal

endmodule initialization
