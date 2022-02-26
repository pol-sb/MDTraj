module function_module
	

	
	implicit none

	! Default seed
    integer, parameter :: defaultsd = 4357
! Period parameters
    integer, parameter :: N = 624, N1 = N + 1

! the array for the state vector
    integer, save, dimension(0:N-1) :: mt
    integer, save                   :: mti = N1
    

	contains
	
!-------------------------------------------------------------------------------
!===============================================================================
!Construction===================================================================
!===============================================================================
!*******************************************************************************
subroutine periodical_init(inc,nbr)
!*******************************************************************************
!subroutine that given a 
!*******************************************************************************
!Input: 
!Output: 
!-------------------------------------------------------------------------------
	implicit none
	integer,intent(out)::inc(:,:),nbr(:,:)
	integer::L,i,x,y
		
	L=size(inc(0,:))
		
	do i=1,L
		inc(0,i)=i-1
		inc(1,i)=i+1
	enddo
	inc(0,1)=L
	inc(1,L)=1
		
	!nbr definition
	!--------------
	i=0
		
	do y=1,L
		do x=1,L
			i=i+1
			nbr(i,1) = inc(1,x) + L*(y-1) 
			nbr(i,2) = inc(0,x) + L*(y-1) 
			nbr(i,3) = x + L*(inc(1,y)-1) 
			nbr(i,4) = x + L*(inc(0,y)-1) 
		enddo
	enddo
		
	return
		
endsubroutine periodical_init


!-----------------------------------------------------------
end module function_module
	
