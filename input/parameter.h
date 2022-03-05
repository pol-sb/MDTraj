
	!	Initialization parameters
	!---------------------------
	integer::nc=2 !Number of particles per side
	!double precision::density=2.5
	double precision::density=10.0d0 !Lomgitude of the size box
	integer::dimension=3 !Dimension of the system.
	
	integer::structure=1 !For SC:1, for fcc:2, for diamond 3.  
	
	double precision::temp= 300.0d0  !Temperature (K)
	integer::vel_opt=1 ! if 1, implemented bimodal. Else, set up to 0
	
	!Simulation set up
	!-----------------
	integer::ntimes=6 !Steps of the simulation
	double precision:: dt= 0.005 !Time step (ps)
	
	
	!Force parameters
	!----------------
	double precision::rc= 0.05d0 !Cut-off
